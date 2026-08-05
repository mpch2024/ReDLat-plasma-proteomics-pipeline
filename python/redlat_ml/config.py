from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
import os
from typing import Iterable


@dataclass(frozen=True)
class MLConfig:
    project_root: Path
    data_dir: Path
    metadata_file: Path
    adat_file: Path
    master_file: Path
    ptau_metadata_file: Path
    matched_ids_file: Path
    matched_full_file: Path
    no_exclusions_file: Path | None
    excluded_ids_file: Path | None
    result_root: Path
    private_root: Path
    public_root: Path
    publication_root: Path
    allow_participant_exports: bool
    matching_focus_id: str


def _flag(name: str, default: bool = False) -> bool:
    value = os.getenv(name, str(default)).strip().lower()
    return value in {"1", "true", "t", "yes", "y"}


def _optional_path(value: str) -> Path | None:
    value = value.strip()
    return Path(value).expanduser().resolve() if value else None


def locate_project_root(script_file: str | Path | None = None) -> Path:
    explicit = os.getenv("REDLAT_PROJECT_ROOT", "").strip()
    if explicit:
        root = Path(explicit).expanduser().resolve()
        if not root.exists():
            raise FileNotFoundError(f"REDLAT_PROJECT_ROOT does not exist: {root}")
        return root

    starts: list[Path] = []
    if script_file is not None:
        starts.append(Path(script_file).resolve().parent)
    starts.append(Path.cwd().resolve())

    for start in starts:
        for candidate in [start, *start.parents]:
            if (candidate / ".redlat-root").exists() or (candidate / ".here").exists():
                return candidate
    raise FileNotFoundError(
        "Project root was not found. Run inside the repository or set REDLAT_PROJECT_ROOT."
    )


def load_config(script_file: str | Path | None = None) -> MLConfig:
    root = locate_project_root(script_file)
    data_dir = Path(os.getenv("REDLAT_ML_DATA_DIR", str(root / "data_private"))).expanduser().resolve()
    metadata = Path(os.getenv("REDLAT_ML_METADATA_FILE", str(data_dir / "clinical_metadata.csv"))).expanduser().resolve()
    adat = Path(os.getenv("REDLAT_ML_ADAT_FILE", str(data_dir / "proteomics.adat"))).expanduser().resolve()
    master = Path(os.getenv("REDLAT_ML_MASTER_FILE", str(data_dir / "ml_master_matrix.csv"))).expanduser().resolve()
    result_root = Path(os.getenv("REDLAT_ML_RESULT_ROOT", str(root / "result" / "ML"))).expanduser().resolve()
    ptau_metadata = Path(os.getenv(
        "REDLAT_ML_PTAU_METADATA_FILE",
        str(metadata),
    )).expanduser().resolve()
    matched_ids = Path(os.getenv(
        "REDLAT_ML_MATCHED_IDS_FILE",
        str(result_root / "private" / "matching" / "matched_ids_SELECTED.csv"),
    )).expanduser().resolve()
    matched_full = Path(os.getenv(
        "REDLAT_ML_MATCHED_FULL_FILE",
        str(result_root / "private" / "matching" / "Matched_Output_SELECTED.csv"),
    )).expanduser().resolve()
    no_exclusions = _optional_path(os.getenv("REDLAT_ML_NO_EXCLUSIONS_FILE", ""))
    excluded_ids = _optional_path(os.getenv("REDLAT_ML_EXCLUDED_IDS_FILE", ""))
    publication_root = Path(os.getenv(
        "REDLAT_ML_PUBLICATION_ROOT",
        str(root / "publication_candidate" / "ML"),
    )).expanduser().resolve()
    private_root = result_root / "private"
    public_root = result_root / "public"
    for directory in (data_dir, result_root, private_root, public_root, publication_root):
        directory.mkdir(parents=True, exist_ok=True)
    return MLConfig(
        project_root=root,
        data_dir=data_dir,
        metadata_file=metadata,
        adat_file=adat,
        master_file=master,
        ptau_metadata_file=ptau_metadata,
        matched_ids_file=matched_ids,
        matched_full_file=matched_full,
        no_exclusions_file=no_exclusions,
        excluded_ids_file=excluded_ids,
        result_root=result_root,
        private_root=private_root,
        public_root=public_root,
        publication_root=publication_root,
        allow_participant_exports=_flag("REDLAT_ALLOW_PARTICIPANT_LEVEL_EXPORTS", False),
        matching_focus_id=os.getenv("REDLAT_ML_MATCHING_FOCUS_ID", "").strip(),
    )


def read_excluded_ids(config: MLConfig) -> set[str]:
    path = config.excluded_ids_file
    if path is None or not path.exists():
        return set()
    values: set[str] = set()
    for line in path.read_text(encoding="utf-8").splitlines():
        value = line.strip().split(",")[0].strip()
        if value and value.lower() not in {"sampleid", "sample_id", "id"}:
            values.add(value)
    return values


def require_files(items: Iterable[tuple[Path, str]]) -> None:
    missing = [f"{label}: {path}" for path, label in items if not path.exists()]
    if missing:
        raise FileNotFoundError("Missing required files:\n" + "\n".join(missing))


def assert_public_table(dataframe, label: str = "table") -> None:
    forbidden = {
        "sampleid", "sample_id", "participantid", "participant_id",
        "subjectid", "subject_id", "record_id", "match_pair_id", "subclass",
    }
    columns = {str(column).strip().lower() for column in dataframe.columns}
    found = sorted(columns & forbidden)
    if found:
        raise ValueError(f"{label} contains direct or stable identifiers: {found}")
