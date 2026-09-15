#!/usr/bin/env python3
"""Reproduce BAFQMC and ED data for the main-text and supplemental benchmarks."""
from __future__ import annotations
import argparse
import json
import os
from pathlib import Path
import platform
import subprocess
import sys
import time

ROOT = Path(__file__).resolve().parent
DATA = ROOT / "benchmarks/paper/data"
MODELS = ("number_conserving", "pairing")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--mode", choices=("archived", "plot", "check", "raw", "full", "dqmc", "ed", "analyze", "smoke"), default="full", help="full (default): rerun BAFQMC and ED; plot: draw the stored processed data; smoke: small installation check")
    parser.add_argument("--model", choices=("both", *MODELS), default="both")
    parser.add_argument("--scope", choices=("all", "main", "supplement"), default="all", help="all (default): main-text U/Delta scans and supplemental U1 scan; intersects --model")
    parser.add_argument("--output", type=Path, help="results directory; default: benchmarks/paper/output/full-TIMESTAMP for production")
    parser.add_argument("--python-ed", default=sys.executable, help="Python interpreter containing QuSpin for fresh ED")
    parser.add_argument("--threads", type=int, default=1, help="threads per numerical-library process (default 1); cases run sequentially, one MPI rank")
    parser.add_argument("--work-dir", type=Path, help="empty Linux scratch directory; default: a new bafqmc directory in the system temp directory")
    parser.add_argument("--resume", action="store_true", help="resume completed per-case stages in the given --output")
    parser.add_argument("--plan", action="store_true", help="print resource estimates and exit without launching computations")
    args = parser.parse_args()
    if args.threads < 1:
        parser.error("--threads must be positive")
    if args.resume and args.output is None:
        parser.error("--resume requires --output pointing to the existing run")
    if args.mode == "smoke" and args.scope != "all":
        parser.error("--scope selects paper benchmarks; use --model to select a smoke pipeline")
    default = ROOT / "benchmarks/paper/output"
    if args.mode == "full":
        default = default / time.strftime("full-%Y%m%d-%H%M%S")
    output = (args.output or default).resolve()
    # Keep generated runs away from tracked source/data and manuscript figures.
    if output == ROOT or output.is_relative_to(DATA) or output.is_relative_to(ROOT / "src") or output.is_relative_to(ROOT / "figures"):
        parser.error("choose a separate generated-output directory")
    models = MODELS if args.model == "both" else (args.model,)
    from benchmarks.paper.manuscript import select_index, selection_summary
    # Planning uses only the standard library, so it also works before installation.
    index = select_index(json.loads((DATA / "index.json").read_text()), models=models, scope=args.scope)
    models = tuple(m for m in models if any(c["model"] == m for c in index["cases"]))
    if args.mode != "smoke":
        print(f"Scope: {args.scope}; {len(index['cases'])} paper points.\n{selection_summary(index['cases'])}", flush=True)
    if args.mode in {"full", "dqmc", "ed"} or args.plan:
        resources = json.loads((ROOT / "benchmarks/paper/resources.json").read_text())
        print(f"Mode: {args.mode}; models: {', '.join(models)}; {args.threads} thread(s), one MPI rank; cases run sequentially.")
        print(resources["planning_summary"])
        print(f"Results: {output}\nSee benchmarks/paper/RESOURCES.md for measurements and stage estimates.", flush=True)
    if args.plan:
        return
    output.mkdir(parents=True, exist_ok=True)
    os.environ.setdefault("MPLCONFIGDIR", str(output / ".matplotlib"))
    os.environ.setdefault("XDG_CACHE_HOME", str(output / ".cache"))
    for name in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS", "NUMBA_NUM_THREADS"):
        os.environ[name] = str(args.threads)
    from benchmarks.paper.analysis import load_index, verify_checksums, processed_cases, recompute, save_tables, plot
    index = select_index(load_index(DATA), models=models, scope=args.scope)
    count = verify_checksums(DATA)
    print(f"Verified {count} archived files.", flush=True)
    if args.mode == "smoke":
        from benchmarks.paper.smoke import smoke
        smoke(output, models, args.python_ed)
        return
    fresh = args.mode in {"full", "dqmc", "ed"}
    if fresh:
        if platform.system() != "Linux":
            parser.error("fresh simulations require Linux/WSL; use scripts/reproduce.ps1 on Windows")
        modes = ("dqmc", "ed") if args.mode == "full" else (args.mode,)
        from benchmarks.paper.production import produce
        produce(ROOT, DATA, index, output, models, modes, args.python_ed, args.threads, resume=args.resume, work_dir=args.work_dir)
        if args.mode != "full":
            print(f"Completed {args.mode}; results: {output / 'runs'}")
            return
    if args.mode in {"archived", "plot", "check"}:
        cases, blocks = processed_cases(index, DATA, models=models)
    else:
        cases, blocks = recompute(index, DATA, fresh_dir=output / "runs" if args.mode in {"full", "analyze"} else None, models=models)
    save_tables(cases, blocks, output)
    if args.mode != "check":
        plot(cases, output / "figures")
    import numpy
    import matplotlib
    provenance = {"mode": args.mode, "scope": args.scope, "models": models, "python": sys.version, "numpy": numpy.__version__, "matplotlib": matplotlib.__version__, "platform": platform.platform(), "cases": len(cases), "case_ids": [f"{c['model']}/{c['id']}" for c in cases], "data_index_sha256": __import__("hashlib").sha256((DATA / "index.json").read_bytes()).hexdigest()}
    (output / "environment.json").write_text(json.dumps(provenance, indent=2) + "\n")
    print(f"Reproduced {len(cases)} benchmark points. Results: {output}")


if __name__ == "__main__":
    try:
        main()
    except (OSError, ValueError, subprocess.CalledProcessError) as exc:
        print(f"error: {exc}", file=sys.stderr)
        raise SystemExit(1)
