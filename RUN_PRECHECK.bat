@echo off
setlocal
cd /d "%~dp0"
echo [1/4] Preparing processed reviewer-data inputs...
python scripts\validation\00_prepare_reviewer_inputs.py
if errorlevel 1 exit /b %errorlevel%
echo [2/4] Static repository audit...
python scripts\validation\01_static_repository_audit.py
if errorlevel 1 exit /b %errorlevel%
echo [3/4] Checking locked manuscript counts...
python scripts\validation\02_check_submitted_locked_counts.py
if errorlevel 1 exit /b %errorlevel%
echo [4/4] Building script dependency map...
python scripts\validation\03_build_dependency_map.py
if errorlevel 1 exit /b %errorlevel%
echo.
echo PRECHECK COMPLETE.
pause
