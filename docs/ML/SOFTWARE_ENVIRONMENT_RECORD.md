# Internal software environment record

Run the environment-capture utility from the repository root after restoring the exact analysis environment and immediately before the final production run:

```bash
python workflow/00_capture_software_environment.py
```

The utility creates an untracked timestamped folder under `internal/software_environment/` containing:

- detected R and Python dependencies and their source files;
- installed package versions available in the active environment;
- R session information and Bioconductor version when `Rscript` is available;
- Python, operating-system and system-tool versions;
- hashes of analysis, workflow and configuration source files;
- hashes of available environment lock files;
- the names, but not values, of configured `REDLAT_*` environment variables.

Absolute home and project paths are redacted by default. Use `--include-paths` only for a governed internal record. The generated directory is excluded from Git and must be reviewed before sharing externally.

Useful options:

```bash
python workflow/00_capture_software_environment.py --all-installed-r
python workflow/00_capture_software_environment.py --all-installed-python
python workflow/00_capture_software_environment.py --rscript "<path-to-Rscript>"
```

The runtime record complements, but does not replace, the repository-level `renv.lock` and the locked Python environment used for the final release.
