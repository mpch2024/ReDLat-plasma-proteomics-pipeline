#!/usr/bin/env bash
set -euo pipefail
cd "$(dirname "$0")"
python scripts/validation/00_prepare_reviewer_inputs.py
python scripts/validation/01_static_repository_audit.py
python scripts/validation/02_check_submitted_locked_counts.py
python scripts/validation/03_build_dependency_map.py
echo "PRECHECK COMPLETE."
