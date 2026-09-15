#!/usr/bin/env python3
"""Run copied example inputs without changing the distributed seed/input files."""
import argparse
from pathlib import Path
import shutil
import subprocess

root = Path(__file__).resolve().parents[1]
parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument('--output', type=Path, default=root / 'build' / 'example')
parser.add_argument('--np', type=int, default=1)
args = parser.parse_args()
output = args.output.resolve()
if output.exists() and any(output.iterdir()):
    parser.error(f'output must be new or empty to avoid mixing Monte Carlo samples: {output}')
output.mkdir(parents=True, exist_ok=True)
for name in ('paramC_sets.txt', 'confin.txt', 'seeds.txt'):
    shutil.copy2(root.parents[1] / 'examples' / 'number_conserving' / 'examples' / 'triangle_3x3_free' / name, output / name)
subprocess.run(['bash', str(root / 'scripts' / 'run_local.sh'), str(output), str(args.np)], check=True)
