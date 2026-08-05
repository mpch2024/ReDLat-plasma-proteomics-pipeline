from pathlib import Path
import subprocess, sys
ROOT = Path(__file__).resolve().parents[1]
subprocess.run([sys.executable, str(ROOT / "workflow/run_ml_analysis.py")], check=True, cwd=ROOT)
subprocess.run([sys.executable, str(ROOT / "workflow/run_ml_reporting.py")], check=True, cwd=ROOT)
