"""Small complete DQMC -> ED -> analysis runs for installation checks."""
from pathlib import Path
import platform
import subprocess
import sys


def smoke(output: Path, models, python_ed: str):
    if platform.system() != "Linux":
        raise ValueError("simulation smoke checks require Linux/WSL")
    root = Path(__file__).resolve().parents[2]
    for model in models:
        solver = root / "src" / model
        filename = "pipeline_smoke.json" if model == "number_conserving" else "triangle_pairing_pipeline_smoke.json"
        manifest = solver / "benchmarks/campaigns" / filename
        subprocess.run(["make", "-C", str(solver), "build"], check=True)
        for mode in ("dqmc", "ed", "analyze"):
            command = [sys.executable, str(solver / "run_paper.py"), "--manifest", str(manifest), "--output", str(output / "smoke" / model), "--mode", mode, "--python", python_ed]
            print("+", " ".join(command), flush=True)
            subprocess.run(command, check=True)
    print(f"Completed {len(models)} small DQMC/ED pipelines: {output / 'smoke'}")
