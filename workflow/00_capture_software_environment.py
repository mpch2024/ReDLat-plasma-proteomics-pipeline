#!/usr/bin/env python3
"""Capture the software environment used by a ReDLat analysis workflow.

The report is written to an untracked internal directory. Absolute user paths
are redacted by default. The script records detected R and Python dependencies,
runtime versions, system tools, environment-variable names and source hashes.
"""

from __future__ import annotations

import argparse
import ast
import csv
import datetime as dt
import hashlib
import importlib.metadata as metadata
import json
import os
from pathlib import Path
import platform
import re
import shutil
import subprocess
import sys
import tempfile
from typing import Iterable

EXCLUDED_DIRS = {
    ".git", "archive", "internal", "result", "results", "publication_candidate",
    "data_private", "data", "__pycache__", ".venv", ".Rproj.user",
}
SOURCE_SUFFIXES = {".R", ".r", ".py", ".yml", ".yaml", ".toml", ".md", ".txt", ".json"}
HASH_SUFFIXES = {".R", ".r", ".py", ".yml", ".yaml", ".toml", ".md", ".txt", ".json", ".csv"}
R_ASSIGNMENT_NAMES = {
    "required_pkgs", "required_packages", "packages", "cran_pkgs",
    "bioc_pkgs", "cran_packages", "bioconductor_packages",
}
DISTRIBUTION_ALIASES = {"skopt": "scikit-optimize", "sklearn": "scikit-learn", "PIL": "Pillow"}

SYSTEM_TOOLS: dict[str, list[list[str]]] = {
    "Rscript": [["Rscript", "--version"]],
    "Python": [[sys.executable, "--version"]],
    "Git": [["git", "--version"]],
    "Conda": [["conda", "--version"]],
    "Pandoc": [["pandoc", "--version"]],
    "Quarto": [["quarto", "--version"]],
    "Java": [["java", "-version"]],
    "Inkscape": [["inkscape", "--version"]],
    "Ghostscript": [["gs", "--version"], ["gswin64c", "--version"], ["gswin32c", "--version"]],
    "ImageMagick": [["magick", "-version"], ["convert", "-version"]],
    "LibreOffice": [["libreoffice", "--version"], ["soffice", "--version"]],
    "GCC": [["gcc", "--version"]],
    "Make": [["make", "--version"]],
}


def locate_root(start: Path) -> Path:
    start = start.resolve()
    for candidate in [start, *start.parents]:
        if (candidate / ".redlat-root").exists() or (candidate / ".here").exists():
            return candidate
    raise FileNotFoundError("Repository root not found. Run the script inside the workflow repository.")


def iter_source_files(root: Path, suffixes: set[str] = SOURCE_SUFFIXES) -> Iterable[Path]:
    for path in root.rglob("*"):
        if not path.is_file() or path.suffix not in suffixes:
            continue
        rel = path.relative_to(root)
        if any(part in EXCLUDED_DIRS for part in rel.parts):
            continue
        yield path


def read_text(path: Path) -> str:
    try:
        return path.read_text(encoding="utf-8")
    except UnicodeDecodeError:
        return path.read_text(encoding="latin-1", errors="replace")


def detect_r_packages(root: Path) -> tuple[set[str], dict[str, list[str]]]:
    packages: set[str] = set()
    evidence: dict[str, list[str]] = {}
    patterns = [
        re.compile(r"(?:library|require)\s*\(\s*['\"]?([A-Za-z][A-Za-z0-9.]*)"),
        re.compile(r"requireNamespace\s*\(\s*['\"]([A-Za-z][A-Za-z0-9.]*)['\"]"),
        re.compile(r"\b([A-Za-z][A-Za-z0-9.]*):::{0,1}[A-Za-z][A-Za-z0-9._]*"),
    ]
    assignment = re.compile(
        r"\b(" + "|".join(map(re.escape, sorted(R_ASSIGNMENT_NAMES))) + r")\s*<-\s*c\((.*?)\)",
        re.S,
    )
    quoted = re.compile(r"['\"]([A-Za-z][A-Za-z0-9.]*)['\"]")

    for path in iter_source_files(root, {".R", ".r"}):
        if path.name == "00_capture_software_environment.py":
            continue
        text = read_text(path)
        rel = str(path.relative_to(root))
        found: set[str] = set()
        for pattern in patterns:
            found.update(pattern.findall(text))
        for match in assignment.finditer(text):
            found.update(quoted.findall(match.group(2)))
        found.discard("base")
        for package in found:
            packages.add(package)
            evidence.setdefault(package, []).append(rel)
    return packages, evidence


def detect_python_imports(root: Path) -> tuple[set[str], dict[str, list[str]]]:
    modules: set[str] = set()
    evidence: dict[str, list[str]] = {}
    for path in iter_source_files(root, {".py"}):
        if path.name == "00_capture_software_environment.py":
            continue
        rel = str(path.relative_to(root))
        try:
            tree = ast.parse(read_text(path), filename=str(path))
        except SyntaxError:
            continue
        found: set[str] = set()
        for node in ast.walk(tree):
            if isinstance(node, ast.Import):
                found.update(alias.name.split(".")[0] for alias in node.names)
            elif isinstance(node, ast.ImportFrom) and node.module:
                found.add(node.module.split(".")[0])
        for module in found:
            modules.add(module)
            evidence.setdefault(module, []).append(rel)
    return modules, evidence


def python_distribution_versions(modules: set[str]) -> list[dict[str, str]]:
    mapping = metadata.packages_distributions()
    stdlib = getattr(sys, "stdlib_module_names", set())
    rows: list[dict[str, str]] = []
    for module in sorted(modules):
        if module in stdlib:
            rows.append({"module": module, "distribution": "Python standard library", "version": platform.python_version(), "status": "available"})
            continue
        if module == "redlat_ml":
            rows.append({"module": module, "distribution": "local project module", "version": "", "status": "available in repository"})
            continue
        distributions = mapping.get(module, []) or [DISTRIBUTION_ALIASES.get(module, module)]
        seen: set[str] = set()
        for distribution in distributions:
            if distribution in seen:
                continue
            seen.add(distribution)
            try:
                version = metadata.version(distribution)
                status = "installed"
            except metadata.PackageNotFoundError:
                version = ""
                status = "not detected"
            rows.append({"module": module, "distribution": distribution, "version": version, "status": status})
    return rows


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def run_command(command: list[str], timeout: int = 15) -> tuple[str, str, int]:
    try:
        completed = subprocess.run(command, capture_output=True, text=True, timeout=timeout, check=False)
        output = (completed.stdout or completed.stderr or "").strip()
        first_line = output.splitlines()[0] if output else ""
        return first_line, output, completed.returncode
    except (FileNotFoundError, PermissionError, subprocess.TimeoutExpired) as exc:
        return "", str(exc), 127


def capture_tools() -> list[dict[str, str]]:
    rows: list[dict[str, str]] = []
    for label, alternatives in SYSTEM_TOOLS.items():
        selected: list[str] | None = None
        version_line = ""
        full_output = ""
        return_code = 127
        executable = ""
        for command in alternatives:
            candidate = command[0]
            resolved = candidate if Path(candidate).exists() else shutil.which(candidate)
            if not resolved:
                continue
            selected = [str(resolved), *command[1:]]
            version_line, full_output, return_code = run_command(selected)
            executable = Path(str(resolved)).name
            break
        rows.append({
            "tool": label,
            "available": str(selected is not None and return_code == 0),
            "executable": executable,
            "version": version_line,
            "return_code": str(return_code),
            "details": full_output[:2000],
        })
    return rows


def redact(value: str, root: Path, include_paths: bool) -> str:
    if include_paths:
        return value
    replacements = {
        str(root.resolve()): "<PROJECT_ROOT>",
        str(Path.home().resolve()): "~",
    }
    result = value
    for source, replacement in sorted(replacements.items(), key=lambda item: len(item[0]), reverse=True):
        result = result.replace(source, replacement)
        result = result.replace(source.replace("\\", "/"), replacement)
    return result


def write_csv(path: Path, rows: list[dict[str, object]], fieldnames: list[str] | None = None) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    if fieldnames is None:
        fieldnames = list(rows[0].keys()) if rows else []
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)


def capture_r(root: Path, output: Path, packages: set[str], rscript_override: str | None, include_all: bool) -> dict[str, object]:
    rscript = rscript_override or shutil.which("Rscript")
    if not rscript:
        return {"available": False, "message": "Rscript was not found on PATH."}

    package_file = output / "r_packages_requested.txt"
    package_file.write_text("\n".join(sorted(packages)), encoding="utf-8")
    r_code = r'''
args <- commandArgs(trailingOnly = TRUE)
out_dir <- args[[1]]
package_file <- args[[2]]
include_all <- identical(args[[3]], "TRUE")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
requested <- if (file.exists(package_file)) unique(readLines(package_file, warn = FALSE)) else character()
requested <- requested[nzchar(requested)]
installed <- installed.packages()
rows <- data.frame(
  package = requested,
  installed = requested %in% rownames(installed),
  version = NA_character_,
  library_path = NA_character_,
  stringsAsFactors = FALSE
)
for (i in seq_len(nrow(rows))) {
  pkg <- rows$package[[i]]
  if (rows$installed[[i]]) {
    rows$version[[i]] <- as.character(packageVersion(pkg))
    rows$library_path[[i]] <- installed[pkg, "LibPath"]
  }
}
write.csv(rows, file.path(out_dir, "r_packages_used.csv"), row.names = FALSE, na = "")
if (include_all) {
  all_rows <- data.frame(
    package = rownames(installed),
    version = installed[, "Version"],
    library_path = installed[, "LibPath"],
    priority = installed[, "Priority"],
    stringsAsFactors = FALSE
  )
  write.csv(all_rows, file.path(out_dir, "r_packages_all_installed.csv"), row.names = FALSE, na = "")
}
runtime <- data.frame(
  key = c("R_version", "platform", "arch", "os", "Bioconductor_version", "renv_version"),
  value = c(
    R.version.string,
    R.version$platform,
    R.version$arch,
    R.version$os,
    if (requireNamespace("BiocManager", quietly = TRUE)) as.character(BiocManager::version()) else NA_character_,
    if (requireNamespace("renv", quietly = TRUE)) as.character(packageVersion("renv")) else NA_character_
  ),
  stringsAsFactors = FALSE
)
write.csv(runtime, file.path(out_dir, "r_runtime.csv"), row.names = FALSE, na = "")
writeLines(capture.output(sessionInfo()), file.path(out_dir, "r_session_info.txt"))
'''
    with tempfile.NamedTemporaryFile("w", suffix=".R", delete=False, encoding="utf-8") as handle:
        handle.write(r_code)
        temp_script = Path(handle.name)
    try:
        completed = subprocess.run(
            [str(rscript), "--vanilla", str(temp_script), str(output), str(package_file), str(include_all).upper()],
            capture_output=True,
            text=True,
            timeout=120,
            check=False,
        )
    finally:
        temp_script.unlink(missing_ok=True)
    (output / "r_capture_console.txt").write_text(
        (completed.stdout or "") + ("\n" + completed.stderr if completed.stderr else ""), encoding="utf-8"
    )
    return {
        "available": completed.returncode == 0,
        "executable": Path(str(rscript)).name,
        "return_code": completed.returncode,
        "message": (completed.stderr or completed.stdout or "").strip()[:2000],
    }


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--project-root", type=Path, default=None)
    parser.add_argument("--output", type=Path, default=None)
    parser.add_argument("--rscript", default=None, help="Optional path to Rscript.")
    parser.add_argument("--include-paths", action="store_true", help="Retain absolute executable and library paths in the internal report.")
    parser.add_argument("--all-installed-r", action="store_true", help="Also record every installed R package.")
    parser.add_argument("--all-installed-python", action="store_true", help="Also record every installed Python distribution.")
    args = parser.parse_args()

    root = locate_root(args.project_root or Path.cwd())
    timestamp = dt.datetime.now(dt.timezone.utc).astimezone().strftime("%Y%m%d_%H%M%S%z")
    output = (args.output or root / "internal" / "software_environment" / timestamp).resolve()
    output.mkdir(parents=True, exist_ok=False)

    r_packages, r_evidence = detect_r_packages(root)
    python_modules, py_evidence = detect_python_imports(root)
    py_rows = python_distribution_versions(python_modules)
    if args.all_installed_python:
        seen = {(row["distribution"], row["version"]) for row in py_rows}
        for dist in metadata.distributions():
            name = dist.metadata.get("Name", "") or ""
            version = dist.version or ""
            if name and (name, version) not in seen:
                py_rows.append({"module": "", "distribution": name, "version": version, "status": "installed (environment)"})
    write_csv(output / "python_packages_used.csv", py_rows, ["module", "distribution", "version", "status"])

    evidence_rows: list[dict[str, str]] = []
    for package, files in sorted(r_evidence.items()):
        evidence_rows.append({"language": "R", "dependency": package, "source_files": " | ".join(sorted(set(files)))})
    for module, files in sorted(py_evidence.items()):
        evidence_rows.append({"language": "Python", "dependency": module, "source_files": " | ".join(sorted(set(files)))})
    write_csv(output / "detected_dependencies.csv", evidence_rows, ["language", "dependency", "source_files"])

    tools = capture_tools()
    for row in tools:
        row["details"] = redact(str(row["details"]), root, args.include_paths)
    write_csv(output / "system_tools.csv", tools, ["tool", "available", "executable", "version", "return_code", "details"])

    r_result = capture_r(root, output, r_packages, args.rscript, args.all_installed_r)
    if not args.include_paths:
        for filename in ("r_packages_used.csv", "r_packages_all_installed.csv", "r_runtime.csv", "r_session_info.txt", "r_capture_console.txt"):
            path = output / filename
            if path.exists():
                path.write_text(redact(read_text(path), root, False), encoding="utf-8")

    hash_rows: list[dict[str, str]] = []
    for path in iter_source_files(root, HASH_SUFFIXES):
        hash_rows.append({
            "relative_path": str(path.relative_to(root)).replace("\\", "/"),
            "sha256": sha256(path),
            "size_bytes": str(path.stat().st_size),
        })
    write_csv(output / "source_file_hashes.csv", hash_rows, ["relative_path", "sha256", "size_bytes"])

    env_names: set[str] = set()
    env_pattern = re.compile(r"\b(REDLAT_[A-Z0-9_]+)\b")
    for path in iter_source_files(root):
        env_names.update(env_pattern.findall(read_text(path)))
    env_rows = [{"variable": name, "set": str(bool(os.getenv(name, ""))), "value_recorded": "False"} for name in sorted(env_names)]
    write_csv(output / "environment_variables.csv", env_rows, ["variable", "set", "value_recorded"])

    lock_files = []
    for name in ("renv.lock", "environment.yml", "requirements.txt", "pyproject.toml", "poetry.lock", "uv.lock"):
        for path in root.rglob(name):
            if any(part in EXCLUDED_DIRS for part in path.relative_to(root).parts):
                continue
            lock_files.append({"relative_path": str(path.relative_to(root)).replace("\\", "/"), "sha256": sha256(path)})
    write_csv(output / "environment_lock_files.csv", lock_files, ["relative_path", "sha256"])

    summary = {
        "captured_at": dt.datetime.now(dt.timezone.utc).isoformat(),
        "project_root": redact(str(root), root, args.include_paths),
        "platform": platform.platform(),
        "machine": platform.machine(),
        "processor": platform.processor(),
        "python": {
            "version": platform.python_version(),
            "implementation": platform.python_implementation(),
            "executable": Path(sys.executable).name if not args.include_paths else sys.executable,
        },
        "R": r_result,
        "detected_R_packages": len(r_packages),
        "detected_Python_modules": len(python_modules),
        "source_files_hashed": len(hash_rows),
        "lock_files_found": lock_files,
        "paths_redacted": not args.include_paths,
    }
    (output / "software_environment_summary.json").write_text(json.dumps(summary, indent=2, ensure_ascii=False), encoding="utf-8")
    (output / "README.txt").write_text(
        "Internal software-environment record\n"
        "====================================\n\n"
        "This directory records the runtime used to execute the workflow. It is ignored by Git.\n"
        "Keep it with the internal analysis archive and provide a reviewed copy only if requested.\n"
        "Absolute user paths are redacted unless --include-paths was explicitly supplied.\n"
        "The presence of a package in the report does not prove that every code branch was executed.\n",
        encoding="utf-8",
    )

    print(f"Software environment captured in: {output}")
    if not r_result.get("available"):
        print("Warning: R runtime details were not captured because Rscript was unavailable or failed.", file=sys.stderr)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
