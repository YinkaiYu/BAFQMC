#!/usr/bin/env python3
"""Check the local BAFQMC environment without installing packages or running jobs."""
from __future__ import annotations

import argparse
import json
import os
import platform
import shlex
import shutil
import subprocess
import sys


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--python-ed", default=sys.executable, help="interpreter containing QuSpin")
    parser.add_argument("--plot-only", action="store_true", help="check only dependencies for processed data and figures")
    parser.add_argument("--json", action="store_true", help="emit a machine-readable report")
    args = parser.parse_args()
    checks = []

    def record(name, ok, detail):
        checks.append({"name": name, "ok": bool(ok), "detail": detail})

    def packages(label, interpreter, names):
        code = (
            "import importlib, json; "
            f"names={names!r}; "
            "print(json.dumps({n: getattr(importlib.import_module(n), '__version__', 'installed') for n in names}))"
        )
        env = os.environ.copy()
        for name in ("OMP_NUM_THREADS", "MKL_NUM_THREADS", "OPENBLAS_NUM_THREADS", "NUMBA_NUM_THREADS"):
            env.setdefault(name, "1")
        try:
            result = subprocess.run([interpreter, "-c", code], capture_output=True, text=True, timeout=45, env=env)
            messages = (result.stderr or result.stdout).strip().splitlines()
            detail = result.stdout.strip() if result.returncode == 0 else (messages[-1] if messages else f"interpreter exited with status {result.returncode}")
            record(label, result.returncode == 0, detail)
        except (OSError, subprocess.TimeoutExpired) as exc:
            record(label, False, str(exc))

    record("Python >= 3.11", sys.version_info >= (3, 11), f"{sys.executable}: {platform.python_version()}")
    packages("analysis packages", sys.executable, ["numpy", "matplotlib"])
    if not args.plot_only:
        record("Linux or WSL", platform.system() == "Linux", platform.platform())
        record("make", shutil.which("make") is not None, shutil.which("make") or "install make")
        compiler = shlex.split(os.environ.get("FC", ""))
        if not compiler:
            compiler = next(([name] for name in ("mpifort", "mpiifx", "mpiifort") if shutil.which(name)), [])
        if compiler:
            try:
                result = subprocess.run(compiler + ["--version"], capture_output=True, text=True, timeout=15)
                detail = (result.stdout or result.stderr).strip().splitlines()
                record("MPI Fortran compiler", result.returncode == 0, detail[0] if detail else shlex.join(compiler))
            except (OSError, subprocess.TimeoutExpired) as exc:
                record("MPI Fortran compiler", False, str(exc))
        else:
            record("MPI Fortran compiler", False, "install gfortran and an MPI development package, or load Intel oneAPI")
        launcher = os.environ.get("MPIEXEC", "mpirun")
        record("MPI launcher", shutil.which(launcher) is not None, shutil.which(launcher) or f"missing {launcher}")
        packages("ED packages", args.python_ed, ["quspin", "numpy", "scipy"])
    report = {"ready": all(c["ok"] for c in checks), "checks": checks}
    if args.json:
        print(json.dumps(report, indent=2))
    else:
        for check in checks:
            print(f"{'OK' if check['ok'] else 'MISSING'}  {check['name']}: {check['detail']}")
        if report["ready"]:
            print("Next: run the small installation check described in docs/getting-started.md.")
        else:
            print("See docs/getting-started.md for environment setup.")
    return 0 if report["ready"] else 1


if __name__ == "__main__":
    raise SystemExit(main())
