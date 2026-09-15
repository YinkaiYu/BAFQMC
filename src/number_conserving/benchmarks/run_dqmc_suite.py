#!/usr/bin/env python3
"""Run every live DQMC benchmark case listed in a suite manifest."""

from __future__ import annotations

import argparse
import json
import os
import shutil
import subprocess
import sys
import tempfile
from pathlib import Path
from typing import Any


REPO_ROOT = Path(__file__).resolve().parents[1]
REQUIRED_INPUTS = ("paramC_sets.txt", "confin.txt", "seeds.txt")


def resolve_repo_path(path: str) -> Path:
    candidate = Path(path)
    if candidate.is_absolute():
        return candidate
    return REPO_ROOT / candidate


def load_suite(path: Path) -> list[dict[str, Any]]:
    suite = json.loads(path.read_text(encoding="utf-8"))
    cases = suite.get("cases")
    if not isinstance(cases, list) or not cases:
        raise ValueError(f"{path} must define a non-empty 'cases' list")
    return cases


def copy_inputs(input_dir: Path, run_dir: Path) -> None:
    if not input_dir.is_dir():
        raise FileNotFoundError(f"benchmark input directory does not exist: {input_dir}")
    for filename in REQUIRED_INPUTS:
        source = input_dir / filename
        if not source.is_file():
            raise FileNotFoundError(f"missing required input file: {source}")
        shutil.copy2(source, run_dir / filename)


def run_command(command: list[str], env: dict[str, str]) -> int:
    print("+ " + " ".join(command), flush=True)
    return subprocess.run(command, cwd=REPO_ROOT, env=env, check=False).returncode


def run_case(case: dict[str, Any], np: str, python: str) -> int:
    name = case["name"]
    input_dir = resolve_repo_path(case["input_dir"])
    reference = resolve_repo_path(case["reference"])
    if not reference.is_file():
        raise FileNotFoundError(f"benchmark reference does not exist: {reference}")

    run_dir = Path(tempfile.mkdtemp(prefix=f"bosonDQMC-{name}.", dir="/tmp"))
    copy_inputs(input_dir, run_dir)

    print("=" * 72, flush=True)
    print(f"DQMC benchmark case: {name}", flush=True)
    print(f"Input: {input_dir}", flush=True)
    print(f"Reference: {reference}", flush=True)
    print(f"Run directory: {run_dir}", flush=True)

    env = os.environ.copy()
    run_rc = run_command(["bash", "scripts/run_local.sh", str(run_dir), np], env)
    if run_rc != 0:
        return run_rc

    compare_rc = run_command(
        [
            python,
            "benchmarks/compare.py",
            "--reference",
            str(reference),
            "--run-dir",
            str(run_dir),
        ],
        env,
    )
    return compare_rc


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--suite",
        type=Path,
        default=REPO_ROOT / "benchmarks" / "dqmc_suite.json",
        help="JSON manifest listing live DQMC benchmark cases",
    )
    parser.add_argument("--np", default=os.environ.get("MPI_NP", "1"))
    parser.add_argument("--python", default=sys.executable)
    args = parser.parse_args()

    try:
        suite_path = args.suite if args.suite.is_absolute() else REPO_ROOT / args.suite
        failures: list[str] = []
        for case in load_suite(suite_path):
            rc = run_case(case, str(args.np), str(args.python))
            if rc != 0:
                failures.append(case["name"])
                break
    except Exception as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        return 2

    if failures:
        print("Failed live DQMC benchmark cases:")
        for name in failures:
            print(f"- {name}")
        return 1

    print("All live DQMC benchmark cases passed.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
