# ReDLat Nature Aging pipeline update overlay

This package contains the public-code update for the DEP, WGCNA and strict machine-learning workflows.

## Important

- Apply this package to a **new Git branch**, not directly to `main`.
- Do not copy participant-level data, local outputs or a populated `.Renviron` file.
- Keep the repository's existing root `README.md`, `LICENSE`, `CITATION.cff` and manuscript-specific metadata unless they are reviewed separately.
- This overlay intentionally does not replace the repository root README or licensing files.

## Included modules

- `scripts/DEP/`
- `scripts/WGCNA/`
- `scripts/ML/`
- shared workflow, configuration, documentation, tests and static checks

## Local validation before merge

1. Copy `.Renviron.example` to `.Renviron` and populate only local paths.
2. Run `python workflow/00_capture_software_environment.py`.
3. Run the DEP, WGCNA and ML static audits.
4. Run full analyses only on the governed local environment.
5. Review `git status` before every commit and confirm that no restricted file is staged.
