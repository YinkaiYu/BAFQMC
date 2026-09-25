#!/usr/bin/env python3
"""Run the paired BAFQMC paper campaign from a portable JSON manifest.

Fresh simulations are written separately from the stored paper data. The
analysis exports observables, standard errors, and differences from ED.
"""
from __future__ import annotations

import argparse
import csv
import json
import os
from pathlib import Path
import shutil
import subprocess
import sys
from typing import Any

ROOT = Path(__file__).resolve().parent
sys.path.insert(0, str(ROOT.parents[1]))
from src.pairing.benchmarks import campaign_pairing_delta as campaign
from src.pairing.benchmarks import pairing_delta_analysis as analysis

MAIN_OBSERVABLES = ("density_total", "energy_density", "S_SF_K", "S_DW_K")


def write_json(path: Path, payload: Any) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8")


def validate_manifest(manifest: dict[str, Any]) -> list[dict[str, Any]]:
    cases = list(campaign.iter_cases(manifest))
    if not cases:
        raise ValueError("manifest must contain at least one case")
    names = [campaign.case_name(case) for case in cases]
    if len(set(names)) != len(names):
        raise ValueError("manifest contains duplicate case names")
    for case in cases:
        params = case["parameters"]
        if float(params.get("t", 1)) != 1.0:
            raise ValueError("the Fortran solver fixes t=1; change RT in src/calc_basic.f90 for other hopping")
        defaults = campaign._merged_defaults(case, "dqmc", manifest.get("dqmc_defaults", {}))
        beta = float(params["beta"])
        dtau = float(defaults.get("dtau", 0.01))
        if beta <= 0 or dtau <= 0 or min(int(params[k]) for k in ("Lx", "Ly")) <= 0:
            raise ValueError("beta, dtau and lattice dimensions must be positive")
        ltrot = int(defaults.get("Ltrot", round(beta / dtau)))
        nwrap = int(defaults["Nwrap"])
        if ltrot <= 0 or nwrap <= 0 or ltrot % nwrap:
            raise ValueError("Ltrot must be a positive multiple of Nwrap")
        if int(defaults["Nbin"]) <= 0 or int(defaults["Nsweep"]) <= 0:
            raise ValueError("Nbin and Nsweep must be positive")
    return cases


def initialize(manifest: dict[str, Any], output: Path, seed: int) -> None:
    snapshot = output / "run_manifest.json"
    payload = {"manifest": manifest, "base_seed": seed, "seed_policy": "base_seed + case_index"}
    if snapshot.exists():
        if json.loads(snapshot.read_text()) != payload:
            raise ValueError(f"{output} belongs to a different manifest or seed; choose a fresh output directory")
        for case in campaign.iter_cases(manifest):
            for name in ("paramC_sets.txt", "confin.txt", "seeds.txt", "params.json"):
                path = output / "inputs" / campaign.case_name(case) / name
                if not path.is_file():
                    raise ValueError(f"incomplete initialization: {path}")
        return
    if (output / "inputs").exists() and any((output / "inputs").iterdir()):
        raise ValueError(f"refusing to replace unowned inputs in {output}; choose a fresh directory")
    campaign.init_local(manifest, output / "inputs", seed=seed)
    write_json(snapshot, payload)


def count_rows(path: Path) -> int:
    with path.open(encoding="utf-8") as handle:
        return sum(bool(line.strip()) for line in handle)


def run_dqmc(cases: list[dict[str, Any]], output: Path, np: int, manifest: dict[str, Any]) -> None:
    # Every run is a fresh chain: the executable appends measurements, so never
    # silently reuse even a partial output directory.
    for case in cases:
        run = output / "inputs" / campaign.case_name(case)
        for name in (*campaign.KNOWN_APPEND_OUTPUTS, "dqmc.log"):
            if (run / name).exists():
                raise ValueError(f"existing simulation output {run / name}; choose a fresh output directory")
    subprocess.run(["make", "build"], cwd=ROOT, check=True)
    launcher = shutil.which(os.environ.get("MPIEXEC", "mpirun"))
    if launcher is None:
        raise FileNotFoundError("MPI launcher not found; activate MPI or set MPIEXEC")
    env = os.environ.copy()
    env.setdefault("OMP_NUM_THREADS", "1")
    env.setdefault("MKL_NUM_THREADS", "1")
    env.setdefault("OPENBLAS_NUM_THREADS", "1")
    for case in cases:
        run = output / "inputs" / campaign.case_name(case)
        print(f"BAFQMC: {run.name}", flush=True)
        with (run / "dqmc.log").open("w") as log:
            subprocess.run([launcher, "-np", str(np), str(ROOT / "build/bosonDQMC.out")],
                           cwd=run, env=env, stdout=log, stderr=subprocess.STDOUT, check=True)
        expected = int(campaign._merged_defaults(case, "dqmc", manifest.get("dqmc_defaults", {}))["Nbin"])
        actual = count_rows(run / "density_total")
        if actual != expected:
            raise ValueError(f"{run.name}: expected {expected} measurements, found {actual}; see dqmc.log")


def run_ed(cases: list[dict[str, Any]], output: Path, python: str, memory_cap: float) -> None:
    env = os.environ.copy()
    env.setdefault("OMP_NUM_THREADS", "1")
    env.setdefault("MKL_NUM_THREADS", "1")
    env.setdefault("OPENBLAS_NUM_THREADS", "1")
    (output / "ed_results").mkdir(parents=True, exist_ok=True)
    for case in cases:
        name = campaign.case_name(case)
        result = output / "ed_results" / f"{name}.json"
        if result.exists():
            raise ValueError(f"existing ED result {result}; choose a fresh output directory")
    # Serial execution preserves the per-process dense-memory guard's meaning.
    for case in cases:
        name = campaign.case_name(case)
        print(f"ED: {name}", flush=True)
        subprocess.run([python, str(ROOT / "benchmarks/ed/ed_pairing_triangle_general.py"),
                        "--params", str(output / "inputs" / name / "params.json"),
                        "--output", str(output / "ed_results" / f"{name}.json"),
                        "--dense-memory-cap-gib", str(memory_cap)], env=env, check=True)


def analyze(manifest: dict[str, Any], cases: list[dict[str, Any]], output: Path) -> dict[str, Any]:
    records = []
    table = []
    for case in cases:
        name = campaign.case_name(case)
        params = case["parameters"]
        defaults = campaign._merged_defaults(case, "dqmc", manifest.get("dqmc_defaults", {}))
        settings = manifest.get("analysis_defaults", {})
        block_size = int(settings.get("block_size", defaults["block_size"]))
        expected = int(defaults["Nbin"])
        run = output / "inputs" / name
        nrows = count_rows(run / "density_total")
        if nrows != expected:
            raise ValueError(f"{name}: expected {expected} rows, got {nrows}")
        ed_path = output / "ed_results" / f"{name}.json"
        ed_payload = json.loads(ed_path.read_text(encoding="utf-8"))
        expected_params = json.loads((run / "params.json").read_text(encoding="utf-8"))
        for key in ("Lx", "Ly", "t", "U1", "U2", "mu", "Delta", "beta"):
            if ed_payload["parameters"][key] != expected_params[key]:
                raise ValueError(f"{name}: ED parameter {key} differs from the manifest")
        for key in ("nmax", "ncut"):
            if ed_payload["cutoffs"][key] != expected_params[key]:
                raise ValueError(f"{name}: ED cutoff {key} differs from the manifest")
        comparison = analysis.compare_dqmc_ed_case(
            run, block_size=block_size, Lx=int(params["Lx"]), Ly=int(params["Ly"]),
            ed_result_path=ed_path, stderr_tolerance=3.0,
            skip_samples=int(settings.get("skip_samples", 0)))
        if comparison["status"] in ("dqmc_unavailable", "ed_unavailable"):
            raise ValueError(f"{name}: {comparison['reason']}")
        comparison.update({"case": name, "Delta": float(params["Delta"]), "row_count": nrows})
        records.append(comparison)
        for observable, row in comparison["observables"].items():
            table.append({"case": name, "Delta": params["Delta"], "observable": observable,
                          "dqmc": row.get("dqmc"), "ed": row.get("ed"),
                          "stderr": row.get("actual_stderr", row.get("stderr")),
                          "reliable": row.get("reliable", False),
                          "pass_uncertainty": row.get("pass_uncertainty", False),
                          "reason": row.get("reason", "")})
    main_pass = all(bool(case["observables"].get(name, {}).get("pass_uncertainty"))
                    and bool(case["observables"].get(name, {}).get("reliable"))
                    for case in records for name in MAIN_OBSERVABLES)
    payload = {"cases": records, "manifest": manifest, "main_observables_pass": main_pass,
               "status": "ok" if all(case["status"] == "ok" for case in records) else "comparison_has_failures"}
    summary = output / "summary"
    if not table:
        raise ValueError("no observables are shared by the DQMC and ED results")
    write_json(summary / "comparison_dqmc_ed_nmax3_ncut4.json", payload)
    with (summary / "comparison_observables.csv").open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(table[0]))
        writer.writeheader()
        writer.writerows(table)
    print(f"Analyzed {len(records)} cases; tables: {summary}", flush=True)
    return payload


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--manifest", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--mode", choices=("init", "dqmc", "ed", "analyze"), required=True)
    parser.add_argument("--case", action="append", dest="case_names", help="run this case name only; repeat to select several cases")
    parser.add_argument("--dry-run", action="store_true", help="print the validated plan without creating or changing files")
    parser.add_argument("--python", default=sys.executable, help="Python interpreter with QuSpin for ED")
    parser.add_argument("--np", type=int, default=1)
    parser.add_argument("--seed", type=int, help="base seed, incremented in manifest case order")
    parser.add_argument("--dense-memory-cap-gib", type=float, default=12.0)
    parser.add_argument("--require-agreement", action="store_true", help="fail analysis if any main observable fails")
    args = parser.parse_args()
    try:
        if args.np < 1 or args.np > 128:
            raise ValueError("--np must be between 1 and 128 (number of initialized rank seeds)")
        if args.dense_memory_cap_gib <= 0:
            raise ValueError("--dense-memory-cap-gib must be positive")
        manifest = campaign.load_manifest(args.manifest)
        cases = validate_manifest(manifest)
        seed = args.seed if args.seed is not None else int(manifest.get("dqmc_defaults", {}).get("seeds", [13579113])[0])
        output = args.output.resolve()
        all_names = {campaign.case_name(case) for case in cases}
        selected_names = set(args.case_names) if args.case_names else all_names
        unknown_names = selected_names - all_names
        if unknown_names:
            raise ValueError(f"unknown case names: {', '.join(sorted(unknown_names))}")
        if args.dry_run:
            print(json.dumps({
                "mode": args.mode, "manifest": str(args.manifest.resolve()),
                "output": str(output), "mpi_ranks": args.np,
                "python": args.python, "dense_memory_cap_gib": args.dense_memory_cap_gib,
                "initializes_all_manifest_inputs": True,
                "cases": [{"case": campaign.case_name(case), "seed": seed + index,
                           "parameters": case["parameters"],
                           "dqmc": campaign._merged_defaults(case, "dqmc", manifest.get("dqmc_defaults", {})),
                           "ed": campaign._merged_defaults(case, "ed", manifest.get("ed_defaults", {}))}
                          for index, case in enumerate(cases) if campaign.case_name(case) in selected_names],
            }, indent=2, sort_keys=True))
            return 0
        # Initialize the full manifest before subsetting so each case keeps its
        # original seed position when jobs are split or resumed individually.
        initialize(manifest, output, seed)
        cases = [case for case in cases if campaign.case_name(case) in selected_names]
        if args.mode == "dqmc":
            run_dqmc(cases, output, args.np, manifest)
        elif args.mode == "ed":
            run_ed(cases, output, args.python, args.dense_memory_cap_gib)
        elif args.mode == "analyze":
            result = analyze(manifest, cases, output)
            if args.require_agreement and not result["main_observables_pass"]:
                return 1
    except (ValueError, OSError, KeyError, subprocess.CalledProcessError) as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        return 2
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
