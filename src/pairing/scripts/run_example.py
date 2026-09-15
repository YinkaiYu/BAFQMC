#!/usr/bin/env python3
"""Copy a committed example into a fresh run directory and execute it."""
import argparse
from pathlib import Path
import shutil
import subprocess

ROOT = Path(__file__).resolve().parents[1]
parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument("--input", type=Path, default=ROOT.parents[1] / "examples/pairing/examples/triangle_pairing_3x2")
parser.add_argument("--output", type=Path, required=True)
parser.add_argument("--np", type=int, default=1)
args = parser.parse_args()
output = args.output.resolve()
if output.exists() and any(output.iterdir()):
    parser.error(f"output directory is not empty: {output}; choose a fresh output directory")
output.mkdir(parents=True, exist_ok=True)
for filename in ("paramC_sets.txt", "confin.txt", "seeds.txt"):
    shutil.copy2(args.input / filename, output / filename)
subprocess.run(["bash", str(ROOT / "scripts/run_local.sh"), str(output), str(args.np)], check=True)
