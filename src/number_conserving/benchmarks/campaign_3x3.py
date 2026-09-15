#!/usr/bin/env python3
"""Import-light helpers for the 3x3 observable benchmark campaign."""

from __future__ import annotations

import argparse
import csv
import html
import json
import math
from pathlib import Path
from typing import Any, Iterable


CAMPAIGN_NAME = "triangle_3x3_observables"
DEFAULT_JOB_PREFIX = "bafqmc/"
DQMC_OUTPUT_FILES = (
    "info.txt",
    "confout.txt",
    "density_up",
    "density_do",
    "density_total",
    "density_site_total",
    "kinetic",
    "doubleOcc",
    "squareOcc",
    "num_up",
    "num_do",
    "numsquare_up",
    "numsquare_do",
    "onsite_n2_up",
    "onsite_n2_do",
    "interaction_energy_density",
    "energy_density",
    "sf_K",
    "psf_Gamma",
    "dw_K",
    "den_updo",
    "pole_z",
    "pole_distance",
    "pole_x",
    "green_spectral_radius",
    "green_smax",
    "log_weight",
)
DQMC_OUTPUT_GLOBS = (
    "den_upup_sub*",
    "den_dodo_sub*",
)
OBSERVABLE_ORDER = (
    "density_total",
    "energy_density",
    "doubleOcc",
    "IPR",
    "S_SF_K",
    "S_PSF_Gamma",
    "S_DW_K",
)
OBSERVABLE_DQMC_FILES = {
    "density_total": ("density_total",),
    "energy_density": ("energy_density",),
    "doubleOcc": ("doubleOcc",),
    "IPR": ("density_site_total",),
    "S_SF_K": ("sf_K",),
    "S_PSF_Gamma": ("psf_Gamma",),
    "S_DW_K": ("dw_K",),
}


def default_manifest() -> dict[str, Any]:
    """Return the reviewed default manifest for the 3x3 observable campaign."""

    return {
        "campaign": CAMPAIGN_NAME,
        "description": "3x3 triangular-lattice BAFQMC vs ED observable benchmark",
        "lattice": {"Lx": 3, "Ly": 3},
        "observables": [
            "density_total",
            "energy_density",
            "doubleOcc",
            "IPR",
            "S_SF_K",
            "S_PSF_Gamma",
            "S_DW_K",
        ],
        "sweeps": [
            {
                "name": "U2_sweep",
                "vary": "U2",
                "mu": -3.5,
                "beta": 4.0,
                "U1": 0.0,
                "U2": [0.0, 0.5, 1.0, 1.5, 2.0],
            },
            {
                "name": "U1_sweep",
                "vary": "U1",
                "mu": -5.0,
                "beta": 1.0,
                "U1": [0.0, -0.05, -0.1, -0.15, -0.2],
                "U2": 1.0,
            },
        ],
        "dqmc_defaults": {
            "dtau": 0.01,
            "Nwrap": 10,
            "Nbin": 100000,
            "Nsweep": 1,
            "shiftLoc": 1.5,
            "is_tau": False,
            "Nthermal": 0,
            "is_warm": True,
            "Nwarm": 500,
            "shiftWarm1": 1.0,
            "shiftWarm2": 1.0,
            "iniType": 2,
            "iniAmpl": 0.1,
            "iniBias1": 0.0,
            "iniBias2": 0.0,
            "block_size": 10000,
            "mpi_np": 1,
            "seeds": [64464988],
        },
        "ed_policy": {
            "tail_tolerance": 1.0e-3,
            "max_total_particles": 6,
            "max_density_total": 1.0,
            "min_completed_shells": 2,
            "statuses": [
                "converged",
                "low_density_cutoff_accepted",
                "incomplete",
            ],
        },
        "analysis_defaults": {
            "stderr_tolerance": 3.0,
            "atol": 0.0,
            "rtol": 0.0,
            "min_ipr_number": 1.0e-8,
        },
        "stability_notes": {
            "u1_zero": "For 3x3 triangular U1=0 runs, use mu < -3.",
            "negative_u1": (
                "Negative U1 is a finite-window/cutoff comparison only; accepted "
                "points use status low_density_cutoff_accepted."
            ),
        },
        "hpc": {
            "host": "cluster.example",
            "ssh": "ssh cluster.example",
            "work_dir": "/path/to/campaign",
            "queues": ["compute", "large-memory"],
            "job_prefix": DEFAULT_JOB_PREFIX,
            "quspin_env": "/path/to/quspin-env",
            "dqmc_env": "/path/to/python-env",
        },
    }


def iter_cases(manifest: dict[str, Any]) -> Iterable[dict[str, Any]]:
    """Yield concrete case dictionaries expanded from manifest sweeps."""

    lattice = manifest["lattice"]
    for sweep in manifest["sweeps"]:
        vary = sweep["vary"]
        values = sweep[vary]
        if not isinstance(values, list):
            raise ValueError(f"sweep {sweep['name']} vary field must be a list")
        for value in values:
            case = {
                "campaign": manifest["campaign"],
                "sweep": sweep["name"],
                "Lx": lattice["Lx"],
                "Ly": lattice["Ly"],
                "beta": float(sweep["beta"]),
                "mu": float(sweep["mu"]),
                "U1": float(value if vary == "U1" else sweep["U1"]),
                "U2": float(value if vary == "U2" else sweep["U2"]),
            }
            case["label"] = _case_label(case)
            yield case


def case_name(case: dict[str, Any]) -> str:
    """Return the stable directory/reference stem for a concrete case."""

    return f"triangle_{int(case['Lx'])}x{int(case['Ly'])}_{_case_label(case)}"


def slurm_job_name(kind: str, case: dict[str, Any]) -> str:
    """Return the compact Slurm job name for a BAFQMC or ED case."""

    return f"{DEFAULT_JOB_PREFIX}{kind}/{_case_label(case)}"


def write_dqmc_case(
    root: Path, case: dict[str, Any], defaults: dict[str, Any], seed: int
) -> Path:
    """Create a DQMC run directory with the fixed runtime input filenames."""

    run_dir = root / case_name(case)
    run_dir.mkdir(parents=True, exist_ok=True)
    clean_dqmc_outputs(run_dir)
    ltrot = _ltrot(case["beta"], defaults["dtau"])
    (run_dir / "paramC_sets.txt").write_text(
        _param_text(case, defaults, ltrot),
        encoding="utf-8",
    )
    (run_dir / "confin.txt").write_text("0\n", encoding="utf-8")
    (run_dir / "seeds.txt").write_text(f"{int(seed)}\n", encoding="utf-8")
    return run_dir


def clean_dqmc_outputs(run_dir: Path) -> None:
    """Remove known append-only DQMC outputs before regenerating a run directory."""

    for name in DQMC_OUTPUT_FILES:
        path = run_dir / name
        if path.exists():
            path.unlink()
    for pattern in DQMC_OUTPUT_GLOBS:
        for path in run_dir.glob(pattern):
            if path.is_file():
                path.unlink()


def write_ed_case(root: Path, case: dict[str, Any], policy: dict[str, Any]) -> Path:
    """Create a 3x3 ED parameter directory compatible with the ED script."""

    run_dir = root / case_name(case)
    run_dir.mkdir(parents=True, exist_ok=True)
    max_total_particles = int(policy["max_total_particles"])
    convergence_policy = {
        key: value for key, value in policy.items() if key != "max_total_particles"
    }
    payload = {
        "parameters": {
            "Lx": int(case["Lx"]),
            "Ly": int(case["Ly"]),
            "beta": float(case["beta"]),
            "mu": float(case["mu"]),
            "U1": float(case["U1"]),
            "U2": float(case["U2"]),
            "t": float(case.get("t", 1.0)),
            "max_total_particles": max_total_particles,
        },
        "convergence_policy": convergence_policy,
    }
    (run_dir / "params.json").write_text(
        json.dumps(payload, ensure_ascii=False, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )
    return run_dir


def write_stage1_report(output_dir: Path, manifest: dict[str, Any]) -> Path:
    """Write the Chinese Stage 1 observable-definition report skeleton."""

    output_dir.mkdir(parents=True, exist_ok=True)
    path = output_dir / "report_stage_1_observable_definitions.html"
    path.write_text(_stage1_html(manifest), encoding="utf-8")
    return path


def write_manifest(path: Path, manifest: dict[str, Any] | None = None) -> Path:
    """Write a campaign manifest JSON file."""

    data = default_manifest() if manifest is None else manifest
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(
        json.dumps(data, ensure_ascii=False, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )
    return path


def write_results_summary(
    input_dir: Path,
    output_dir: Path,
    manifest: dict[str, Any],
    *,
    block_size: int | None = None,
    skip_samples: int = 0,
    min_ipr_number: float = 1.0e-8,
    stderr_tolerance: float = 3.0,
    atol: float = 0.0,
    rtol: float = 0.0,
    report_name: str = "report_stage_results.html",
) -> dict[str, Path]:
    """Write fail-closed JSON/CSV/HTML result summaries for campaign cases."""

    try:
        from . import campaign_analysis as ca
    except ModuleNotFoundError:
        import campaign_analysis as ca

    output_dir.mkdir(parents=True, exist_ok=True)
    lq = int(manifest["lattice"]["Lx"]) * int(manifest["lattice"]["Ly"])
    effective_block_size = int(
        block_size
        if block_size is not None
        else manifest["dqmc_defaults"]["block_size"]
    )

    records: list[dict[str, Any]] = []
    for case in iter_cases(manifest):
        run_dir = input_dir / case_name(case)
        comparison = ca.compare_dqmc_ed_case(
            run_dir,
            block_size=effective_block_size,
            lq=lq,
            expected_parameters={
                "Lx": int(case["Lx"]),
                "Ly": int(case["Ly"]),
                "beta": float(case["beta"]),
                "mu": float(case["mu"]),
                "U1": float(case["U1"]),
                "U2": float(case["U2"]),
            },
            expected_policy=manifest["ed_policy"],
            skip_samples=skip_samples,
            min_ipr_number=min_ipr_number,
            stderr_tolerance=stderr_tolerance,
            atol=atol,
            rtol=rtol,
        )
        records.append(_compact_case_record(case, comparison))

    figures = _write_result_figures(output_dir, records, manifest)
    submitted_jobs = _load_submitted_jobs(output_dir)
    payload = {
        "campaign": manifest["campaign"],
        "input_dir": str(input_dir),
        "block_size": effective_block_size,
        "skip_samples": skip_samples,
        "min_ipr_number": min_ipr_number,
        "stderr_tolerance": stderr_tolerance,
        "atol": atol,
        "rtol": rtol,
        "cases": records,
        "counts": _status_counts(records),
        "figures": figures,
        "submitted_jobs": submitted_jobs,
        "job_counts": _job_counts(submitted_jobs),
    }

    json_path = output_dir / "summary_campaign_results.json"
    csv_path = output_dir / "summary_campaign_results.csv"
    html_path = output_dir / report_name
    json_path.write_text(
        json.dumps(payload, ensure_ascii=False, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )
    _write_results_csv(csv_path, records)
    html_path.write_text(_results_html(manifest, payload), encoding="utf-8")
    return {"json": json_path, "csv": csv_path, "html": html_path}


def write_queue_snapshot_report(
    summary_dir: Path,
    output_dir: Path,
    *,
    checked_at: str,
    hpc_root: str,
    repo_commit: str,
    default_state: str = "unknown",
    default_reason: str = "",
    default_start_time: str = "",
    default_time_limit: str = "",
    dqmc_output_files: int = 0,
    ed_results_json: int = 0,
    report_name: str = "report_stage_4_hpc_queue_snapshot.html",
) -> Path:
    """Write a Chinese queue-state snapshot report from summary artifacts."""

    output_dir.mkdir(parents=True, exist_ok=True)
    submitted_jobs = _load_submitted_jobs(summary_dir)
    summary_payload = _load_optional_json(summary_dir / "summary_campaign_results.json")
    jobs = [
        _queue_job_record(
            job,
            default_state=default_state,
            default_reason=default_reason,
            default_start_time=default_start_time,
            default_time_limit=default_time_limit,
        )
        for job in submitted_jobs
    ]
    payload = {
        "checked_at": checked_at,
        "hpc_root": hpc_root,
        "repo_commit": repo_commit,
        "summary_dir": str(summary_dir),
        "summary_counts": summary_payload.get("counts", {}),
        "job_counts": summary_payload.get("job_counts", _job_counts(submitted_jobs)),
        "jobs": jobs,
        "state_counts": _queue_state_counts(jobs),
        "production_output_counts": {
            "dqmc_output_files": int(dqmc_output_files),
            "ed_results_json": int(ed_results_json),
        },
    }
    path = output_dir / report_name
    path.write_text(_queue_snapshot_html(payload), encoding="utf-8")
    return path


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    subparsers = parser.add_subparsers(dest="command", required=True)

    write_manifest_parser = subparsers.add_parser("write-manifest")
    write_manifest_parser.add_argument("--output", type=Path, required=True)

    stage1_parser = subparsers.add_parser("stage1-report")
    stage1_parser.add_argument("--manifest", type=Path)
    stage1_parser.add_argument("--output-dir", type=Path, required=True)

    init_parser = subparsers.add_parser("init-local")
    init_parser.add_argument("--manifest", type=Path)
    init_parser.add_argument("--output-dir", type=Path, required=True)
    init_parser.add_argument("--seed", type=int)

    summary_parser = subparsers.add_parser("summarize-results")
    summary_parser.add_argument("--manifest", type=Path)
    summary_parser.add_argument("--input-dir", type=Path, required=True)
    summary_parser.add_argument("--output-dir", type=Path, required=True)
    summary_parser.add_argument("--block-size", type=int)
    summary_parser.add_argument("--skip-samples", type=int, default=0)
    summary_parser.add_argument("--min-ipr-number", type=float)
    summary_parser.add_argument("--stderr-tolerance", type=float)
    summary_parser.add_argument("--atol", type=float)
    summary_parser.add_argument("--rtol", type=float)
    summary_parser.add_argument(
        "--report-name",
        default="report_stage_results.html",
    )

    queue_parser = subparsers.add_parser("queue-snapshot-report")
    queue_parser.add_argument("--summary-dir", type=Path, required=True)
    queue_parser.add_argument("--output-dir", type=Path, required=True)
    queue_parser.add_argument("--checked-at", required=True)
    queue_parser.add_argument("--hpc-root", required=True)
    queue_parser.add_argument("--repo-commit", required=True)
    queue_parser.add_argument("--default-state", default="unknown")
    queue_parser.add_argument("--default-reason", default="")
    queue_parser.add_argument("--default-start-time", default="")
    queue_parser.add_argument("--default-time-limit", default="")
    queue_parser.add_argument("--dqmc-output-files", type=int, default=0)
    queue_parser.add_argument("--ed-results-json", type=int, default=0)
    queue_parser.add_argument(
        "--report-name",
        default="report_stage_4_hpc_queue_snapshot.html",
    )

    args = parser.parse_args(argv)
    if args.command == "write-manifest":
        write_manifest(args.output)
        return 0

    manifest = _load_manifest(args.manifest) if getattr(args, "manifest", None) else default_manifest()
    if args.command == "stage1-report":
        write_stage1_report(args.output_dir, manifest)
        return 0
    if args.command == "init-local":
        defaults = manifest["dqmc_defaults"]
        seed = int(args.seed if args.seed is not None else defaults["seeds"][0])
        for case in iter_cases(manifest):
            write_dqmc_case(args.output_dir, case, defaults, seed)
            write_ed_case(args.output_dir, case, manifest["ed_policy"])
        return 0
    if args.command == "summarize-results":
        analysis_defaults = manifest.get("analysis_defaults", {})
        paths = write_results_summary(
            args.input_dir,
            args.output_dir,
            manifest,
            block_size=args.block_size,
            skip_samples=args.skip_samples,
            min_ipr_number=(
                args.min_ipr_number
                if args.min_ipr_number is not None
                else float(analysis_defaults.get("min_ipr_number", 1.0e-8))
            ),
            stderr_tolerance=(
                args.stderr_tolerance
                if args.stderr_tolerance is not None
                else float(analysis_defaults.get("stderr_tolerance", 3.0))
            ),
            atol=(
                args.atol
                if args.atol is not None
                else float(analysis_defaults.get("atol", 0.0))
            ),
            rtol=(
                args.rtol
                if args.rtol is not None
                else float(analysis_defaults.get("rtol", 0.0))
            ),
            report_name=args.report_name,
        )
        for kind, path in paths.items():
            print(f"{kind}: {path}")
        return 0
    if args.command == "queue-snapshot-report":
        path = write_queue_snapshot_report(
            args.summary_dir,
            args.output_dir,
            checked_at=args.checked_at,
            hpc_root=args.hpc_root,
            repo_commit=args.repo_commit,
            default_state=args.default_state,
            default_reason=args.default_reason,
            default_start_time=args.default_start_time,
            default_time_limit=args.default_time_limit,
            dqmc_output_files=args.dqmc_output_files,
            ed_results_json=args.ed_results_json,
            report_name=args.report_name,
        )
        print(f"html: {path}")
        return 0

    raise AssertionError(f"unhandled command {args.command}")


def _load_manifest(path: Path) -> dict[str, Any]:
    return json.loads(path.read_text(encoding="utf-8"))


def _load_optional_json(path: Path) -> dict[str, Any]:
    if not path.exists():
        return {}
    payload = json.loads(path.read_text(encoding="utf-8"))
    if not isinstance(payload, dict):
        raise ValueError(f"expected JSON object: {path}")
    return payload


def _format_number(value: Any) -> str:
    number = float(value)
    if math.isclose(number, round(number), abs_tol=1.0e-12):
        return str(int(round(number)))
    return f"{number:g}"


def _case_label(case: dict[str, Any]) -> str:
    sweep = str(case["sweep"])
    if sweep.startswith("U2"):
        return (
            f"U2_{_format_number(case['U2'])}"
            f"_mu{_format_number(case['mu'])}"
            f"_b{_format_number(case['beta'])}"
        )
    if sweep.startswith("U1"):
        return (
            f"U1_{_format_number(case['U1'])}"
            f"_mu{_format_number(case['mu'])}"
            f"_b{_format_number(case['beta'])}"
        )
    raise ValueError(f"unknown sweep {sweep!r}")


def _ltrot(beta: float, dtau: float) -> int:
    raw = float(beta) / float(dtau)
    ltrot = int(round(raw))
    if not math.isclose(raw, ltrot, rel_tol=0.0, abs_tol=1.0e-9):
        raise ValueError(f"beta={beta} is not divisible by dtau={dtau}")
    return ltrot


def _fortran_bool(value: bool) -> str:
    return ".true." if value else ".false."


def _row(*values: Any) -> str:
    return "".join(f"{value!s:<12}" for value in values).rstrip()


def _param_text(case: dict[str, Any], defaults: dict[str, Any], ltrot: int) -> str:
    rows = [
        _row(case["U1"], case["U2"], case["mu"]),
        _row(case["Lx"], case["Ly"], ltrot, case["beta"]),
        _row(case["Lx"], case["Ly"], ltrot),
        _row(defaults["Nwrap"], defaults["Nbin"], defaults["Nsweep"], defaults["shiftLoc"]),
        _row(_fortran_bool(defaults["is_tau"]), defaults["Nthermal"]),
        _row(
            _fortran_bool(defaults["is_warm"]),
            defaults["Nwarm"],
            defaults["shiftWarm1"],
            defaults["shiftWarm2"],
        ),
        _row(
            defaults["iniType"],
            defaults["iniAmpl"],
            defaults["iniBias1"],
            defaults["iniBias2"],
        ),
        "",
        "RU1         RU2         mu",
        "Nlx         Nly         Ltrot       Beta",
        "NlxTherm    NlyTherm    LtrotTherm",
        "Nwrap       Nbin        Nsweep      shiftLoc",
        "is_tau      Nthermal",
        "is_warm     Nwarm       shiftWarm1  shiftWarm2",
        "iniType     iniAmpl     iniBias1    iniBias2",
    ]
    return "\n".join(rows) + "\n"


def _stage1_html(manifest: dict[str, Any]) -> str:
    defaults = manifest["dqmc_defaults"]
    ed_policy = manifest["ed_policy"]
    hpc = manifest["hpc"]
    rows = "\n".join(_sweep_rows(manifest))
    observables = "\n".join(
        f"<li><code>{html.escape(name)}</code></li>"
        for name in manifest["observables"]
    )
    queues = ", ".join(html.escape(queue) for queue in hpc["queues"])
    cases = list(iter_cases(manifest))
    case_rows = "\n".join(
        "<tr>"
        f"<td>{html.escape(case['label'])}</td>"
        f"<td>{case['U1']}</td><td>{case['U2']}</td>"
        f"<td>{case['mu']}</td><td>{case['beta']}</td>"
        "</tr>"
        for case in cases
    )
    return f"""<!doctype html>
<html lang="zh-CN">
<head>
  <meta charset="utf-8">
  <meta name="viewport" content="width=device-width, initial-scale=1">
  <title>BAFQMC vs ED 3x3 Observable Campaign - Stage 1</title>
  <script>
    window.MathJax = {{tex: {{inlineMath: [['$', '$'], ['\\\\(', '\\\\)']]}}}};
  </script>
  <script defer src="https://cdn.jsdelivr.net/npm/mathjax@3/es5/tex-chtml.js"></script>
  <style>
    body {{ font-family: system-ui, -apple-system, BlinkMacSystemFont, "Segoe UI", sans-serif; margin: 2rem auto; max-width: 1100px; line-height: 1.65; color: #202124; padding: 0 1rem; }}
    h1, h2 {{ line-height: 1.25; }}
    table {{ border-collapse: collapse; width: 100%; margin: 1rem 0; }}
    th, td {{ border: 1px solid #d0d7de; padding: 0.45rem 0.6rem; text-align: left; }}
    th {{ background: #f6f8fa; }}
    code {{ background: #f6f8fa; padding: 0.1rem 0.25rem; }}
    .scroll {{ overflow-x: auto; }}
  </style>
</head>
<body>
  <h1>BAFQMC vs ED：3x3 Observable Campaign Stage 1</h1>
  <p>本报告是三角晶格 3x3 可观测量基准的第一阶段骨架，目标是在正式长时间运行前固定模型、参数、归一化和集群执行约定。当前还没有生产结果；后续 DQMC 和 ED 作业将基于这里的 manifest 生成输入。</p>

  <h2>模型与目标</h2>
  <p>模型为双味 Bose-Hubbard 哈密顿量，$L_x=L_y=3$，$L_q=9$。比较对象为 <strong>BAFQMC vs ED</strong>，重点检查低密度窗口内的密度、能量、占据和结构因子是否一致。</p>
  <div class="scroll">
  $$H = t\\sum_{{\\langle ij\\rangle}}(b_i^\\dagger b_j+c_i^\\dagger c_j+\\mathrm{{h.c.}})
  +U_1\\sum_i(n_{{b,i}}+n_{{c,i}})^2
  +U_2\\sum_i(n_{{b,i}}-n_{{c,i}})^2.$$
  </div>

  <h2>目标可观测量</h2>
  <ul>{observables}</ul>
  <div class="scroll">
  $$\\rho = \\langle N_b+N_c\\rangle/L_q,$$
  $$e = (\\langle E_{{\\rm kin}}\\rangle + \\langle E_{{\\rm int}}\\rangle)/L_q,$$
  $$D = \\frac1{{L_q}}\\sum_i \\langle n_{{b,i}}n_{{c,i}}\\rangle,$$
  $$S_{{\\mathrm{{SF}}}}(K) = \\frac1{{L_q^2}}\\sum_{{ij}} e^{{iK\\cdot(r_i-r_j)}}\\langle b_i^\\dagger b_j+c_i^\\dagger c_j\\rangle,$$
  $$S_{{\\mathrm{{PSF}}}}(\\Gamma) = \\frac1{{L_q^2}}\\sum_{{ij}}\\langle b_i^\\dagger c_i^\\dagger c_j b_j\\rangle,$$
  $$S_{{\\mathrm{{DW}}}}(K) = \\frac1{{L_q^2}}\\sum_{{ij}} e^{{iK\\cdot(r_i-r_j)}}\\langle (n_{{b,i}}+n_{{c,i}})(n_{{b,j}}+n_{{c,j}})\\rangle,$$
  $$\\mathrm{{IPR}}_\\rho = \\frac{{\\sum_i\\rho_i^2}}{{\\left(\\sum_i\\rho_i\\right)^2}},\\qquad
    \\rho_i=\\langle n_{{b,i}}+n_{{c,i}}\\rangle.$$
  </div>
  <p>结构因子使用标准 $S_O(q)=\\langle O_q^\\dagger O_q\\rangle/L_q^2$ 归一化；3x3 的分母为 81。晶格约定为 $a_1=(1,0)$、$a_2=(1/2,\\sqrt{{3}}/2)$，reciprocal basis 采用 $2\\pi$ 约定，$b_1=2\\pi(1,-1/\\sqrt{{3}})$、$b_2=2\\pi(0,2/\\sqrt{{3}})$。3x3 的 K 点为 $K=(2/3)b_1+(1/3)b_2=(4\\pi/3,0)$。IPR 从总密度分布 <code>density_site_total</code> 后处理得到，均匀 3x3 分布给出 $1/9$，单站点局域极限给出 1。</p>

  <h2>稳定性说明</h2>
  <p>3x3 三角晶格当前跃迁约定下单粒子能量下界为 $-3$。因此 $U_1=0$ 的 3x3 grand-canonical 点必须满足 mu < -3；本 manifest 的 U2 sweep 使用 $\\mu=-3.5$。</p>
  <p>负 $U_1$ 的总密度通道在无限粒子壳层下不是严格有界问题。本阶段只定义低密度有限窗口比较，ED 状态可以是 <code>converged</code>、<code>low_density_cutoff_accepted</code> 或 <code>incomplete</code>；负 $U_1$ 接受点需要明确标为 <code>low_density_cutoff_accepted</code>。</p>

  <h2>参数表</h2>
  <table>
    <thead><tr><th>Sweep</th><th>mu</th><th>beta</th><th>U1</th><th>U2</th></tr></thead>
    <tbody>{rows}</tbody>
  </table>
  <table>
    <thead><tr><th>Case</th><th>U1</th><th>U2</th><th>mu</th><th>beta</th></tr></thead>
    <tbody>{case_rows}</tbody>
  </table>

  <h2>DQMC 生产默认值</h2>
  <p><code>dtau={defaults['dtau']}</code>，<code>Ltrot=beta/dtau</code>，<code>Nbin={defaults['Nbin']}</code>，<code>Nsweep={defaults['Nsweep']}</code>，<code>Nwrap={defaults['Nwrap']}</code>，<code>shiftLoc={defaults['shiftLoc']}</code>，<code>block_size={defaults['block_size']}</code>，<code>mpi_np={defaults['mpi_np']}</code>。Warm-up 使用 <code>is_warm={defaults['is_warm']}</code>，<code>Nwarm={defaults['Nwarm']}</code>，<code>shiftWarm1={defaults['shiftWarm1']}</code>，<code>shiftWarm2={defaults['shiftWarm2']}</code>。</p>

  <h2>ED 截断与状态</h2>
  <p>当前 ED 壳层起点为 <code>max_total_particles={ed_policy['max_total_particles']}</code>，尾项阈值为 <code>{ed_policy['tail_tolerance']}</code>，负 $U_1$ 低密度窗口要求 <code>density_total &lt;= {ed_policy['max_density_total']}</code> 且至少完成 <code>{ed_policy['min_completed_shells']}</code> 个非零壳层。后处理的论文可信门槛还会进一步要求负 $U_1$ finite-window 点满足 <code>density_total &lt;= 5e-2</code> 且 <code>completed_shells &gt;= 6</code>；更浅的早停 checkpoint 只能作为诊断，不能进入 trusted ED 曲线。</p>

  <h2>HPC 记录</h2>
  <p>登录命令：<code>{html.escape(hpc['ssh'])}</code>。工作目录：<code>{html.escape(hpc['work_dir'])}</code>。队列：<code>{queues}</code>。Slurm job prefix：<code>{html.escape(hpc['job_prefix'])}</code>，例如 <code>bafqmc/BAFQMC/U2_1.5_mu-3.5_b4</code> 和 <code>bafqmc/ED/U2_1.5_mu-3.5_b4</code>。</p>
  <p>QuSpin 环境：<code>{html.escape(hpc['quspin_env'])}</code>；DQMC/分析环境：<code>{html.escape(hpc['dqmc_env'])}</code>。</p>
</body>
</html>
"""


def _sweep_rows(manifest: dict[str, Any]) -> Iterable[str]:
    for sweep in manifest["sweeps"]:
        yield (
            "<tr>"
            f"<td>{html.escape(sweep['name'])}</td>"
            f"<td>{sweep['mu']}</td>"
            f"<td>{sweep['beta']}</td>"
            f"<td>{html.escape(str(sweep['U1']))}</td>"
            f"<td>{html.escape(str(sweep['U2']))}</td>"
            "</tr>"
        )


def _compact_case_record(
    case: dict[str, Any], comparison: dict[str, Any]
) -> dict[str, Any]:
    ed_payload = comparison.get("ed", {})
    if not isinstance(ed_payload, dict):
        ed_payload = {}
    ed_reliability = comparison.get("ed_reliability", {})
    if (
        not isinstance(ed_reliability, dict)
        or "kind" not in ed_reliability
        or "reason" not in ed_reliability
    ):
        ed_reliability = {
            "reliable": False,
            "kind": "unavailable",
            "reason": comparison.get("reason", "ed_unavailable"),
        }

    record = {
        "case_name": case_name(case),
        "label": case["label"],
        "sweep": case["sweep"],
        "Lx": int(case["Lx"]),
        "Ly": int(case["Ly"]),
        "U1": float(case["U1"]),
        "U2": float(case["U2"]),
        "mu": float(case["mu"]),
        "beta": float(case["beta"]),
        "run_dir": comparison["run_dir"],
        "ed_result_path": comparison["ed_result_path"],
        "raw_paths": _raw_paths(comparison["run_dir"]),
        "trusted": bool(comparison.get("trusted", False)),
        "status": str(comparison.get("status", "unknown")),
        "reason": str(comparison.get("reason", "")),
        "ed_status": str(ed_payload.get("status", "")),
        "ed_completed_shells": ed_payload.get("completed_shells"),
        "ed_reliability": ed_reliability,
        "relative_last_shell": ed_payload.get("relative_last_shell", {}),
        "observables": comparison.get("observables", {}),
    }
    return record


def _raw_paths(run_dir_value: str) -> dict[str, Any]:
    run_dir = Path(run_dir_value)
    return {
        "run_dir": str(run_dir),
        "paramC_sets": str(run_dir / "paramC_sets.txt"),
        "confin": str(run_dir / "confin.txt"),
        "seeds": str(run_dir / "seeds.txt"),
        "ed_results": str(run_dir / "results.json"),
        "dqmc_files": {
            observable: [str(run_dir / name) for name in names]
            for observable, names in OBSERVABLE_DQMC_FILES.items()
        },
    }


def _status_counts(records: list[dict[str, Any]]) -> dict[str, int]:
    counts: dict[str, int] = {}
    for record in records:
        status = str(record["status"])
        counts[status] = counts.get(status, 0) + 1
    return counts


def _load_submitted_jobs(output_dir: Path) -> list[dict[str, Any]]:
    path = output_dir / "submitted_jobs.json"
    if not path.exists():
        return []
    payload = json.loads(path.read_text(encoding="utf-8"))
    if not isinstance(payload, list):
        raise ValueError(f"submitted_jobs.json must contain a list: {path}")
    return [_submitted_job_record(item) for item in payload]


def _submitted_job_record(item: Any) -> dict[str, Any]:
    if not isinstance(item, dict):
        raise ValueError(f"submitted job record must be an object: {item!r}")
    job_id = str(item.get("job_id", ""))
    run_dir = str(item.get("run_dir", ""))
    logs_root = _infer_hpc_logs_root(run_dir)
    record = {
        "job_id": job_id,
        "case_name": str(item.get("case_name", "")),
        "kind": str(item.get("kind", "")),
        "job_name": str(item.get("job_name", "")),
        "partition": str(item.get("partition", "")),
        "time_limit": str(item.get("time_limit", "")),
        "run_dir": run_dir,
        "sbatch_stdout": str(item.get("sbatch_stdout", "")),
    }
    if job_id and logs_root:
        record["stdout_path"] = str(logs_root / f"{job_id}.out")
        record["stderr_path"] = str(logs_root / f"{job_id}.err")
    return record


def _infer_hpc_logs_root(run_dir_value: str) -> Path | None:
    if not run_dir_value:
        return None
    run_dir = Path(run_dir_value)
    if run_dir.parent.name == "production_inputs":
        return run_dir.parent.parent / "logs"
    return None


def _job_counts(records: list[dict[str, Any]]) -> dict[str, int]:
    counts: dict[str, int] = {}
    for record in records:
        kind = str(record.get("kind", "unknown")) or "unknown"
        counts[kind] = counts.get(kind, 0) + 1
    return counts


def _queue_job_record(
    record: dict[str, Any],
    *,
    default_state: str,
    default_reason: str,
    default_start_time: str,
    default_time_limit: str,
) -> dict[str, Any]:
    job = dict(record)
    job["state"] = str(job.get("state", default_state) or default_state)
    job["reason"] = str(job.get("reason", default_reason) or default_reason)
    job["start_time"] = str(job.get("start_time", default_start_time) or default_start_time)
    job["time_limit"] = str(job.get("time_limit", default_time_limit) or default_time_limit)
    return job


def _queue_state_counts(records: list[dict[str, Any]]) -> dict[str, int]:
    counts: dict[str, int] = {}
    for record in records:
        state = str(record.get("state", "unknown")) or "unknown"
        counts[state] = counts.get(state, 0) + 1
    return counts


def _write_results_csv(path: Path, records: list[dict[str, Any]]) -> None:
    fieldnames = [
        "case_name",
        "sweep",
        "U1",
        "U2",
        "mu",
        "beta",
        "run_dir",
        "ed_result_path",
        "paramC_sets_path",
        "dqmc_files",
        "observable",
        "dqmc",
        "stderr",
        "ed",
        "difference",
        "z_score",
        "threshold",
        "imag_actual",
        "imag_stderr",
        "imag_z_score",
        "imag_threshold",
        "pass_imaginary",
        "pass_uncertainty",
        "observable_reliable",
        "observable_reason",
        "case_trusted",
        "case_status",
        "case_reason",
        "ed_status",
        "ed_completed_shells",
        "ed_reliability_kind",
        "ed_reliability_reason",
    ]
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        for record in records:
            base = {
                "case_name": record["case_name"],
                "sweep": record["sweep"],
                "U1": record["U1"],
                "U2": record["U2"],
                "mu": record["mu"],
                "beta": record["beta"],
                "run_dir": record["run_dir"],
                "ed_result_path": record["ed_result_path"],
                "paramC_sets_path": record["raw_paths"]["paramC_sets"],
                "case_trusted": record["trusted"],
                "case_status": record["status"],
                "case_reason": record["reason"],
                "ed_status": record["ed_status"],
                "ed_completed_shells": record["ed_completed_shells"],
                "ed_reliability_kind": record["ed_reliability"]["kind"],
                "ed_reliability_reason": record["ed_reliability"]["reason"],
            }
            observables = record.get("observables", {})
            for observable in OBSERVABLE_ORDER:
                values = observables.get(observable, {})
                writer.writerow(
                    {
                        **base,
                        "observable": observable,
                        "dqmc_files": ";".join(
                            record["raw_paths"]["dqmc_files"].get(observable, [])
                        ),
                        "dqmc": values.get("dqmc"),
                        "stderr": values.get("stderr"),
                        "ed": values.get("ed"),
                        "difference": values.get("difference"),
                        "z_score": values.get("z_score"),
                        "threshold": values.get("threshold"),
                        "imag_actual": values.get("imag_actual"),
                        "imag_stderr": values.get("imag_stderr"),
                        "imag_z_score": values.get("imag_z_score"),
                        "imag_threshold": values.get("imag_threshold"),
                        "pass_imaginary": values.get("pass_imaginary"),
                        "pass_uncertainty": values.get("pass_uncertainty"),
                        "observable_reliable": values.get("reliable"),
                        "observable_reason": values.get("reason"),
                    }
                )


def _write_result_figures(
    output_dir: Path,
    records: list[dict[str, Any]],
    manifest: dict[str, Any],
) -> list[dict[str, Any]]:
    try:
        import matplotlib

        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
    except Exception as exc:
        return [
            {
                "status": "unavailable",
                "reason": f"matplotlib_unavailable: {exc}",
            }
        ]

    figure_dir = output_dir / "figures"
    figure_dir.mkdir(parents=True, exist_ok=True)
    figures: list[dict[str, Any]] = []
    observables = list(manifest.get("observables", OBSERVABLE_ORDER))
    for sweep in manifest["sweeps"]:
        sweep_name = str(sweep["name"])
        vary = str(sweep["vary"])
        sweep_records = [record for record in records if record["sweep"] == sweep_name]
        for observable in observables:
            figure = _plot_observable_sweep(
                plt,
                figure_dir,
                output_dir,
                sweep_records,
                sweep_name=sweep_name,
                vary=vary,
                observable=observable,
            )
            if figure is not None:
                figures.append(figure)
    return figures


def _plot_observable_sweep(
    plt: Any,
    figure_dir: Path,
    output_dir: Path,
    records: list[dict[str, Any]],
    *,
    sweep_name: str,
    vary: str,
    observable: str,
) -> dict[str, Any] | None:
    trusted_points: list[tuple[float, float, float, float]] = []
    diagnostic_points: list[tuple[float, float, float, float]] = []
    for record in records:
        values = record.get("observables", {}).get(observable)
        if not isinstance(values, dict):
            continue
        x_value = _finite_or_none(record.get(vary))
        dqmc = _finite_or_none(values.get("dqmc"))
        ed = _finite_or_none(values.get("ed"))
        stderr = _finite_or_none(values.get("stderr"))
        if x_value is None or dqmc is None or ed is None or stderr is None:
            continue
        point = (x_value, dqmc, ed, stderr)
        if bool(record.get("trusted")) and bool(values.get("reliable")):
            trusted_points.append(point)
        else:
            diagnostic_points.append(point)

    if not trusted_points and not diagnostic_points:
        return None

    trusted_points.sort(key=lambda item: item[0])
    diagnostic_points.sort(key=lambda item: item[0])
    fig, ax = plt.subplots(figsize=(5.6, 3.7), constrained_layout=True)
    if trusted_points:
        xs = [item[0] for item in trusted_points]
        dqmc_values = [item[1] for item in trusted_points]
        ed_values = [item[2] for item in trusted_points]
        errors = [item[3] for item in trusted_points]
        ax.errorbar(
            xs,
            dqmc_values,
            yerr=errors,
            marker="o",
            color="#1f77b4",
            capsize=3,
            linewidth=1.5,
            label="BAFQMC trusted",
        )
        ax.plot(
            xs,
            ed_values,
            marker="s",
            color="#d62728",
            linewidth=1.3,
            label="ED trusted",
        )
    if diagnostic_points:
        xs = [item[0] for item in diagnostic_points]
        dqmc_values = [item[1] for item in diagnostic_points]
        ed_values = [item[2] for item in diagnostic_points]
        ax.scatter(
            xs,
            dqmc_values,
            marker="x",
            color="#7f7f7f",
            label="BAFQMC excluded",
        )
        ax.scatter(
            xs,
            ed_values,
            marker="+",
            color="#7f7f7f",
            label="ED excluded",
        )
    ax.set_xlabel(vary)
    ax.set_ylabel(observable)
    ax.set_title(f"{sweep_name}: {observable}")
    ax.grid(True, alpha=0.25)
    ax.legend(fontsize=8)

    filename = f"{_safe_filename(sweep_name)}_{_safe_filename(observable)}.png"
    path = figure_dir / filename
    fig.savefig(path, dpi=180)
    plt.close(fig)
    return {
        "status": "written",
        "sweep": sweep_name,
        "vary": vary,
        "observable": observable,
        "path": str(path),
        "relative_path": str(path.relative_to(output_dir)),
        "trusted_points": len(trusted_points),
        "diagnostic_points": len(diagnostic_points),
    }


def _results_html(manifest: dict[str, Any], payload: dict[str, Any]) -> str:
    rows = "\n".join(_result_case_rows(payload["cases"]))
    counts = ", ".join(
        f"{html.escape(status)}={count}" for status, count in payload["counts"].items()
    )
    job_counts = ", ".join(
        f"{html.escape(kind)}={count}"
        for kind, count in payload.get("job_counts", {}).items()
    )
    figures = "\n".join(_figure_blocks(payload.get("figures", [])))
    if not figures:
        figures = "<p>当前没有可绘制曲线：所有点仍缺少 DQMC/ED 数值，或尚未通过 trusted gate。</p>"
    observable_table = "\n".join(_observable_status_rows(payload["cases"]))
    raw_index_table = "\n".join(_raw_index_rows(payload["cases"]))
    job_index_table = "\n".join(_job_index_rows(payload.get("submitted_jobs", [])))
    if not job_index_table:
        job_index_table = (
            "<tr><td colspan=\"9\">未在输出目录找到 submitted_jobs.json；"
            "本报告只索引物理输入和结果文件。</td></tr>"
        )
    observable_definitions = r"""
  <div class="scroll">
    $$\rho = \langle N_b+N_c\rangle/L_q,$$
    $$e = (\langle E_{\rm kin}\rangle + \langle E_{\rm int}\rangle)/L_q,$$
    $$D = L_q^{-1}\sum_i \langle n_{b,i}n_{c,i}\rangle,$$
    $$\mathrm{IPR}_{\rho} = \frac{\sum_i \rho_i^2}{\left(\sum_i \rho_i\right)^2},
      \qquad \rho_i=\langle n_{b,i}+n_{c,i}\rangle,$$
    $$S_{\rm SF}(K)=L_q^{-2}\sum_{ij}e^{iK\cdot(r_i-r_j)}
      \langle b_i^\dagger b_j+c_i^\dagger c_j\rangle,$$
    $$S_{\rm PSF}(\Gamma)=L_q^{-2}\sum_{ij}
      \langle b_i^\dagger c_i^\dagger c_j b_j\rangle,$$
    $$S_{\rm DW}(K)=L_q^{-2}\sum_{ij}e^{iK\cdot(r_i-r_j)}
      \langle (n_{b,i}+n_{c,i})(n_{b,j}+n_{c,j})\rangle.$$
  </div>
"""
    return f"""<!doctype html>
<html lang="zh-CN">
<head>
  <meta charset="utf-8">
  <meta name="viewport" content="width=device-width, initial-scale=1">
  <title>BAFQMC vs ED 3x3 Observable Campaign - Results Summary</title>
  <script>
    window.MathJax = {{tex: {{inlineMath: [['$', '$'], ['\\\\(', '\\\\)']]}}}};
  </script>
  <script defer src="https://cdn.jsdelivr.net/npm/mathjax@3/es5/tex-chtml.js"></script>
  <style>
    body {{ font-family: system-ui, -apple-system, BlinkMacSystemFont, "Segoe UI", sans-serif; margin: 2rem auto; max-width: 1200px; line-height: 1.6; color: #202124; padding: 0 1rem; }}
    table {{ border-collapse: collapse; width: 100%; margin: 1rem 0; font-size: 0.92rem; }}
    th, td {{ border: 1px solid #d0d7de; padding: 0.4rem 0.55rem; text-align: left; vertical-align: top; }}
    th {{ background: #f6f8fa; }}
    code {{ background: #f6f8fa; padding: 0.1rem 0.25rem; }}
    figure {{ margin: 1rem 0 1.5rem; }}
    img {{ max-width: 100%; height: auto; border: 1px solid #d0d7de; }}
    .scroll {{ overflow-x: auto; }}
  </style>
</head>
<body>
  <h1>BAFQMC vs ED：3x3 Results Summary</h1>
  <p>本报告由 <code>benchmarks/campaign_3x3.py summarize-results</code> 生成。它使用统一的后处理规则：DQMC 误差来自 block mean 的标准误，ED 结果必须先通过 <code>ed_reference_reliable()</code>；未完成、dry-run、截断不可靠或 DQMC 输出不足的点都不会进入 trusted comparison。</p>
  <p>输入目录：<code>{html.escape(payload['input_dir'])}</code>。block size：<code>{payload['block_size']}</code>；skip samples：<code>{payload['skip_samples']}</code>；stderr tolerance：<code>{payload['stderr_tolerance']}</code>；atol：<code>{payload['atol']}</code>；rtol：<code>{payload['rtol']}</code>。ED 结果还会检查 <code>Lx</code>、<code>Ly</code>、<code>beta</code>、<code>mu</code>、<code>U1</code>、<code>U2</code> 和 <code>max_total_particles</code> 是否与 manifest 对应 case 一致。</p>
  <p>状态计数：{counts}</p>
  <p>提交作业计数：{job_counts if job_counts else "未记录 submitted_jobs.json"}</p>
  <h2>Manifest</h2>
  <p>Campaign：<code>{html.escape(manifest['campaign'])}</code>；lattice：{manifest['lattice']['Lx']}x{manifest['lattice']['Ly']}。</p>
  <h2>Observable Definitions</h2>
  {observable_definitions}
  <p>结构因子的 DQMC 文件为复数；报告用实部与 ED 比较，并要求虚部在同一误差门槛下与零相容。reciprocal basis 采用 $2\\pi$ 约定；当前三角格写作 $a_1=(1,0)$、$a_2=(1/2,\\sqrt{{3}}/2)$，所以 3x3 的 K 点为 $K=(4\\pi/3,0)$。标准结构因子在精确平均中应为实数且非负。</p>
  <h2>Curves</h2>
  <p>曲线只把 trusted 点作为主比较；有数值但被可靠性或误差门槛排除的点会以灰色诊断符号显示，不参与结论。</p>
  {figures}
  <h2>HPC Job Index</h2>
  <p>该表来自 <code>submitted_jobs.json</code>。它记录实际提交的 Slurm job id、job name、队列、时间限制、case 目录，以及按本 campaign 目录约定推断的 stdout/stderr 日志路径。</p>
  <div class="scroll">
    <table>
      <thead>
        <tr><th>Job ID</th><th>Kind</th><th>Case</th><th>Job name</th><th>Partition</th><th>Time limit</th><th>Run dir</th><th>stdout</th><th>stderr</th></tr>
      </thead>
      <tbody>{job_index_table}</tbody>
    </table>
  </div>
  <h2>Case Status</h2>
  <div class="scroll">
    <table>
      <thead>
        <tr><th>Case</th><th>Sweep</th><th>U1</th><th>U2</th><th>mu</th><th>beta</th><th>Status</th><th>Trusted</th><th>ED status</th><th>ED reliability</th><th>Reason</th></tr>
      </thead>
      <tbody>{rows}</tbody>
    </table>
  </div>
  <h2>Observable Table</h2>
  <div class="scroll">
    <table>
      <thead>
        <tr><th>Case</th><th>Observable</th><th>DQMC</th><th>stderr</th><th>ED</th><th>diff</th><th>z</th><th>Imag</th><th>Pass real / imag</th><th>Trusted</th><th>Reason</th></tr>
      </thead>
      <tbody>{observable_table}</tbody>
    </table>
  </div>
  <h2>Raw Index</h2>
  <ul>
    <li><code>summary_campaign_results.json</code>：case 级别结构化汇总。</li>
    <li><code>summary_campaign_results.csv</code>：逐 observable 展开的 DQMC/ED/difference/z-score 表。</li>
    <li><code>submitted_jobs.json</code>：每个 Slurm 作业的 job id、case、运行目录和提交命令回显。</li>
    <li><code>figures/</code>：每个 sweep 和 observable 的 PNG 曲线图；只有存在数值点时才生成。</li>
    <li>每个 case 原始数据保存在输入目录的对应 <code>triangle_3x3_...</code> 子目录中。</li>
  </ul>
  <div class="scroll">
    <table>
      <thead>
        <tr><th>Case</th><th>Run dir</th><th>ED result</th><th>Input files</th><th>DQMC observable files</th></tr>
      </thead>
      <tbody>{raw_index_table}</tbody>
    </table>
  </div>
</body>
</html>
"""


def _queue_snapshot_html(payload: dict[str, Any]) -> str:
    state_counts = ", ".join(
        f"{html.escape(state)}={count}"
        for state, count in payload.get("state_counts", {}).items()
    )
    job_counts = ", ".join(
        f"{html.escape(kind)}={count}"
        for kind, count in payload.get("job_counts", {}).items()
    )
    summary_counts = ", ".join(
        f"{html.escape(status)}={count}"
        for status, count in payload.get("summary_counts", {}).items()
    )
    job_rows = "\n".join(_queue_job_rows(payload.get("jobs", [])))
    if not job_rows:
        job_rows = "<tr><td colspan=\"10\">没有 submitted_jobs.json 记录。</td></tr>"
    output_counts = payload.get("production_output_counts", {})
    return f"""<!doctype html>
<html lang="zh-CN">
<head>
  <meta charset="utf-8">
  <meta name="viewport" content="width=device-width, initial-scale=1">
  <title>BAFQMC vs ED 3x3 Observable Campaign - Stage 4 Queue Snapshot</title>
  <style>
    body {{ font-family: system-ui, -apple-system, BlinkMacSystemFont, "Segoe UI", sans-serif; margin: 2rem auto; max-width: 1200px; line-height: 1.6; color: #202124; padding: 0 1rem; }}
    h1, h2 {{ line-height: 1.25; }}
    table {{ border-collapse: collapse; width: 100%; margin: 1rem 0; font-size: 0.92rem; }}
    th, td {{ border: 1px solid #d0d7de; padding: 0.4rem 0.55rem; text-align: left; vertical-align: top; }}
    th {{ background: #f6f8fa; }}
    code {{ background: #f6f8fa; padding: 0.1rem 0.25rem; }}
    .scroll {{ overflow-x: auto; }}
    .note {{ color: #57606a; }}
  </style>
</head>
<body>
  <h1>BAFQMC vs ED：3x3 HPC Queue Snapshot</h1>
  <p>本报告记录 3x3 production campaign 的 Slurm 等待/运行状态。它不是最终物理分析；目的只是固定当前队列、日志、输出和 summary 索引状态，避免之后回看时无法区分“未完成”和“失败”。</p>

  <h2>Snapshot</h2>
  <table>
    <tbody>
      <tr><th>Checked at</th><td><code>{html.escape(str(payload.get('checked_at', '')))}</code></td></tr>
      <tr><th>HPC root</th><td><code>{html.escape(str(payload.get('hpc_root', '')))}</code></td></tr>
      <tr><th>Repo commit</th><td><code>{html.escape(str(payload.get('repo_commit', '')))}</code></td></tr>
      <tr><th>Summary dir</th><td><code>{html.escape(str(payload.get('summary_dir', '')))}</code></td></tr>
      <tr><th>Case status counts</th><td><code>{html.escape(summary_counts or 'none')}</code></td></tr>
      <tr><th>Job kind counts</th><td><code>{html.escape(job_counts or 'none')}</code></td></tr>
      <tr><th>Job state counts</th><td><code>{html.escape(state_counts or 'none')}</code></td></tr>
      <tr><th>Production outputs</th><td><code>dqmc_output_files={output_counts.get('dqmc_output_files', 0)}</code>, <code>ed_results_json={output_counts.get('ed_results_json', 0)}</code></td></tr>
    </tbody>
  </table>

  <h2>Submitted Jobs</h2>
  <div class="scroll">
    <table>
      <thead>
        <tr><th>Job ID</th><th>Kind</th><th>Case</th><th>State</th><th>Reason</th><th>Start time</th><th>Time limit</th><th>Run dir</th><th>stdout</th><th>stderr</th></tr>
      </thead>
      <tbody>{job_rows}</tbody>
    </table>
  </div>

  <h2>Indexed Artifacts</h2>
  <ul>
    <li><code>summary_campaign_results.json</code>：case 级别结构化汇总。</li>
    <li><code>summary_campaign_results.csv</code>：逐 observable 展开的 DQMC/ED/difference/z-score 表。</li>
    <li><code>report_stage_3_hpc_pending_summary.html</code>：fail-closed 结果页，包含 case status、raw data path 和 Slurm job index。</li>
    <li><code>submitted_jobs.json</code>：原始 Slurm 提交记录。</li>
  </ul>

  <p class="note">下一步：继续轮询这些 job。只有在 DQMC 输出和 ED <code>results.json</code> 出现后，才运行 trusted summary、曲线生成和最终中文物理报告。</p>
</body>
</html>
"""


def _queue_job_rows(records: list[dict[str, Any]]) -> Iterable[str]:
    for record in records:
        yield (
            "<tr>"
            f"<td>{html.escape(str(record.get('job_id', '')))}</td>"
            f"<td>{html.escape(str(record.get('kind', '')))}</td>"
            f"<td>{html.escape(str(record.get('case_name', '')))}</td>"
            f"<td>{html.escape(str(record.get('state', '')))}</td>"
            f"<td>{html.escape(str(record.get('reason', '')))}</td>"
            f"<td>{html.escape(str(record.get('start_time', '')))}</td>"
            f"<td>{html.escape(str(record.get('time_limit', '')))}</td>"
            f"<td>{html.escape(str(record.get('run_dir', '')))}</td>"
            f"<td>{html.escape(str(record.get('stdout_path', '')))}</td>"
            f"<td>{html.escape(str(record.get('stderr_path', '')))}</td>"
            "</tr>"
        )


def _result_case_rows(records: list[dict[str, Any]]) -> Iterable[str]:
    for record in records:
        reliability = record["ed_reliability"]
        yield (
            "<tr>"
            f"<td>{html.escape(record['case_name'])}</td>"
            f"<td>{html.escape(record['sweep'])}</td>"
            f"<td>{record['U1']}</td>"
            f"<td>{record['U2']}</td>"
            f"<td>{record['mu']}</td>"
            f"<td>{record['beta']}</td>"
            f"<td>{html.escape(record['status'])}</td>"
            f"<td>{record['trusted']}</td>"
            f"<td>{html.escape(record['ed_status'])}</td>"
            f"<td>{html.escape(str(reliability.get('kind', '')))}: "
            f"{html.escape(str(reliability.get('reason', '')))}</td>"
            f"<td>{html.escape(record['reason'])}</td>"
            "</tr>"
        )


def _observable_status_rows(records: list[dict[str, Any]]) -> Iterable[str]:
    for record in records:
        observables = record.get("observables", {})
        if not observables:
            yield (
                "<tr>"
                f"<td>{html.escape(record['case_name'])}</td>"
                "<td colspan=\"10\">"
                f"{html.escape(record['status'])}: {html.escape(record['reason'])}"
                "</td></tr>"
            )
            continue
        for observable in OBSERVABLE_ORDER:
            values = observables.get(observable)
            if not isinstance(values, dict):
                continue
            yield (
                "<tr>"
                f"<td>{html.escape(record['case_name'])}</td>"
                f"<td>{html.escape(observable)}</td>"
                f"<td>{_fmt(values.get('dqmc'))}</td>"
                f"<td>{_fmt(values.get('stderr'))}</td>"
                f"<td>{_fmt(values.get('ed'))}</td>"
                f"<td>{_fmt(values.get('difference'))}</td>"
                f"<td>{_fmt(values.get('z_score'))}</td>"
                f"<td>{_fmt(values.get('imag_actual'))} / "
                f"{_fmt(values.get('imag_stderr'))}</td>"
                f"<td>{html.escape(str(values.get('pass_uncertainty')))} / "
                f"{html.escape(str(values.get('pass_imaginary')))}</td>"
                f"<td>{html.escape(str(values.get('reliable')))}</td>"
                f"<td>{html.escape(str(values.get('reason')))}</td>"
                "</tr>"
            )


def _raw_index_rows(records: list[dict[str, Any]]) -> Iterable[str]:
    for record in records:
        raw_paths = record.get("raw_paths", {})
        dqmc_files = raw_paths.get("dqmc_files", {})
        dqmc_text = "<br>".join(
            f"{html.escape(observable)}: "
            + html.escape("; ".join(files))
            for observable, files in dqmc_files.items()
        )
        input_files = "<br>".join(
            html.escape(str(raw_paths.get(name, "")))
            for name in ("paramC_sets", "confin", "seeds")
        )
        yield (
            "<tr>"
            f"<td>{html.escape(record['case_name'])}</td>"
            f"<td>{html.escape(str(raw_paths.get('run_dir', record['run_dir'])))}</td>"
            f"<td>{html.escape(str(raw_paths.get('ed_results', record['ed_result_path'])))}</td>"
            f"<td>{input_files}</td>"
            f"<td>{dqmc_text}</td>"
            "</tr>"
        )


def _job_index_rows(records: list[dict[str, Any]]) -> Iterable[str]:
    for record in records:
        yield (
            "<tr>"
            f"<td>{html.escape(str(record.get('job_id', '')))}</td>"
            f"<td>{html.escape(str(record.get('kind', '')))}</td>"
            f"<td>{html.escape(str(record.get('case_name', '')))}</td>"
            f"<td>{html.escape(str(record.get('job_name', '')))}</td>"
            f"<td>{html.escape(str(record.get('partition', '')))}</td>"
            f"<td>{html.escape(str(record.get('time_limit', '')))}</td>"
            f"<td>{html.escape(str(record.get('run_dir', '')))}</td>"
            f"<td>{html.escape(str(record.get('stdout_path', '')))}</td>"
            f"<td>{html.escape(str(record.get('stderr_path', '')))}</td>"
            "</tr>"
        )


def _figure_blocks(figures: list[dict[str, Any]]) -> Iterable[str]:
    for figure in figures:
        if figure.get("status") != "written":
            yield (
                "<p>"
                f"Figure generation skipped: {html.escape(str(figure.get('reason', '')))}"
                "</p>"
            )
            continue
        caption = (
            f"{figure['sweep']} / {figure['observable']} "
            f"(trusted={figure['trusted_points']}, diagnostics={figure['diagnostic_points']})"
        )
        yield (
            "<figure>"
            f"<img src=\"{html.escape(str(figure['relative_path']))}\" "
            f"alt=\"{html.escape(caption)}\">"
            f"<figcaption>{html.escape(caption)}</figcaption>"
            "</figure>"
        )


def _finite_or_none(value: Any) -> float | None:
    try:
        number = float(value)
    except (TypeError, ValueError):
        return None
    if not math.isfinite(number):
        return None
    return number


def _fmt(value: Any) -> str:
    number = _finite_or_none(value)
    if number is None:
        return ""
    return f"{number:.8g}"


def _safe_filename(value: str) -> str:
    return "".join(char if char.isalnum() or char in ("-", "_") else "_" for char in value)


if __name__ == "__main__":
    raise SystemExit(main())
