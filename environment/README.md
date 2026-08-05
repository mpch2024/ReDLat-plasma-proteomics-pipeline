# Python environment

Create the strict rerun environment from `environment.yml`. Python 3.11 is specified to provide a conservative cross-platform environment for scikit-learn and scikit-optimize.

After the first successful local run:

1. capture the installed environment with `workflow/00_capture_software_environment.py`;
2. export an exact lock file from the validated environment;
3. rerun the aggregate checkpoints;
4. retain the lock and environment report in the governed internal archive;
5. add the validated lock to the journal release only after the strict outputs are approved.

The package ranges in this folder are installation constraints, not a claim that the historical submission used those exact versions.
