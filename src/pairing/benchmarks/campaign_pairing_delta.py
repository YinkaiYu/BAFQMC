#!/usr/bin/env python3
"""Manifest and local-input helpers for the 3x3 pairing Delta campaign."""

from __future__ import annotations

import argparse
import csv
import html
import json
import math
from pathlib import Path
from typing import Any, Iterable, Mapping, Sequence


CAMPAIGN_NAME = "triangle_pairing_delta_3x3"
DEFAULT_DELTAS = [0.0, 0.1, 0.2, 0.3, 0.4]
PARAMC_DATA_ROWS = 7
KNOWN_APPEND_OUTPUTS = (
    "density",
    "density_up",
    "density_do",
    "density_total",
    "density_site_total",
    "num_up",
    "num_do",
    "kinetic",
    "doubleOcc",
    "squareOcc",
    "local_numsquare",
    "numsquare_up",
    "numsquare_do",
    "pair_equal",
    "den_upup_sub11",
    "den_dodo_sub11",
    "den_updo",
    "energy_density",
    "interaction_energy_density",
    "pairing_energy_density",
    "chemical_energy_density",
    "grand_energy_density",
    "onsite_n2_up",
    "onsite_n2_do",
    "sf_K",
    "dw_K",
    "psf_Gamma",
)
KNOWN_OUTPUT_GLOBS = (
    "den_upup_sub*",
    "den_dodo_sub*",
)


def default_manifest(
    *,
    Lx: int = 3,
    Ly: int = 3,
    U1: float = 1.0,
    U2: float = 0.0,
    beta: float = 4.0,
    mu: float = -3.5,
    deltas: Sequence[float] | None = None,
    Nbin: int = 100000,
    block_size: int = 10000,
    nmax: int = 2,
    ncut: int | None = 6,
    basis_cap: int = 200000,
) -> dict[str, Any]:
    """Return the initial 3x3 Delta-sweep manifest."""

    delta_values = [float(value) for value in (DEFAULT_DELTAS if deltas is None else deltas)]
    return {
        "campaign": CAMPAIGN_NAME,
        "description": "Triangular pairing DQMC/ED Delta sweep on a 3x3 lattice",
        "lattice": {"Lx": int(Lx), "Ly": int(Ly)},
        "parameters": {
            "Lx": int(Lx),
            "Ly": int(Ly),
            "U1": float(U1),
            "U2": float(U2),
            "beta": float(beta),
            "mu": float(mu),
            "t": 1.0,
        },
        "deltas": delta_values,
        "observables": [
            "density_total",
            "energy_density",
            "interaction_energy_density",
            "pairing_energy_density",
            "chemical_energy_density",
            "grand_energy_density",
            "onsite_n2_up",
            "onsite_n2_do",
            "IPR",
            "S_SF_K",
            "S_DW_K",
            "S_PSF_Gamma",
        ],
        "dqmc_defaults": {
            "dtau": 0.01,
            "Nwrap": 10,
            "Nbin": int(Nbin),
            "Nsweep": 1,
            "shiftLoc": 1.0,
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
            "block_size": int(block_size),
            "seeds": [13579113],
        },
        "ed_defaults": {
            "nmax": int(nmax),
            "ncut": None if ncut is None else int(ncut),
            "basis_cap": int(basis_cap),
        },
        "analysis_defaults": {
            "block_size": int(block_size),
            "skip_samples": 0,
            "min_ipr_total_density": 1.0e-8,
        },
    }


def iter_cases(manifest: Mapping[str, Any]) -> Iterable[dict[str, Any]]:
    """Yield concrete Delta cases from a manifest."""

    if isinstance(manifest.get("cases"), list):
        for case in manifest["cases"]:
            params = dict(_case_parameters(case))
            params.setdefault("Lx", manifest.get("lattice", {}).get("Lx", 3))
            params.setdefault("Ly", manifest.get("lattice", {}).get("Ly", 3))
            payload = dict(case)
            payload["parameters"] = params
            payload.setdefault("name", case_name(payload))
            yield payload
        return

    base = dict(manifest.get("parameters", {}))
    lattice = manifest.get("lattice", {})
    base.setdefault("Lx", lattice.get("Lx", 3))
    base.setdefault("Ly", lattice.get("Ly", 3))
    for delta in manifest["deltas"]:
        params = dict(base)
        params["Delta"] = float(delta)
        case = {
            "campaign": manifest.get("campaign", CAMPAIGN_NAME),
            "name": "",
            "parameters": params,
            "dqmc": dict(manifest.get("dqmc_defaults", {})),
            "ed": dict(manifest.get("ed_defaults", {})),
        }
        case["name"] = case_name(case)
        yield case


def case_name(case: Mapping[str, Any]) -> str:
    """Return a stable run-directory name for one Delta case."""

    params = _case_parameters(case)
    lx = int(params["Lx"])
    ly = int(params["Ly"])
    delta = _format_label_float(float(_param(params, "Delta", "delta", "RDelta")))
    mu = _format_label_float(float(_param(params, "mu", "Mu")))
    beta = _format_label_float(float(_param(params, "beta", "Beta")))
    return f"triangle_pairing_delta_{lx}x{ly}_D{delta}_mu{mu}_b{beta}"


def write_dqmc_case(
    case: Mapping[str, Any],
    run_dir: str | Path,
    dqmc_defaults: Mapping[str, Any] | None = None,
    *,
    seed: int | None = None,
) -> Path:
    """Write fixed-filename DQMC inputs for one case."""

    run_path = Path(run_dir)
    run_path.mkdir(parents=True, exist_ok=True)
    clean_dqmc_outputs(run_path)

    params = _case_parameters(case)
    defaults = _merged_defaults(case, "dqmc", dqmc_defaults)
    if seed is None:
        seeds = defaults.get("seeds", [13579113])
        seed = int(seeds[0] if isinstance(seeds, Sequence) and not isinstance(seeds, str) else seeds)

    (run_path / "paramC_sets.txt").write_text(
        _paramc_text(params, defaults),
        encoding="utf-8",
    )
    (run_path / "confin.txt").write_text("0\n", encoding="utf-8")
    (run_path / "seeds.txt").write_text(_seeds_text(int(seed)), encoding="utf-8")
    return run_path


def clean_dqmc_outputs(run_dir: str | Path) -> None:
    """Remove known append-only DQMC outputs while preserving unrelated files."""

    run_path = Path(run_dir)
    for name in KNOWN_APPEND_OUTPUTS:
        path = run_path / name
        if path.is_file():
            path.unlink()
    for pattern in KNOWN_OUTPUT_GLOBS:
        for path in run_path.glob(pattern):
            if path.is_file():
                path.unlink()


def write_ed_case(
    case: Mapping[str, Any],
    run_dir: str | Path,
    ed_defaults: Mapping[str, Any] | None = None,
) -> Path:
    """Write a flat ED params JSON compatible with current pairing ED scripts."""

    run_path = Path(run_dir)
    run_path.mkdir(parents=True, exist_ok=True)
    params = _case_parameters(case)
    defaults = _merged_defaults(case, "ed", ed_defaults)
    payload = {
        "case": case_name(case),
        "Lx": int(params["Lx"]),
        "Ly": int(params["Ly"]),
        "t": float(params.get("t", 1.0)),
        "U1": float(_param(params, "U1")),
        "U2": float(_param(params, "U2")),
        "mu": float(_param(params, "mu", "Mu")),
        "Delta": float(_param(params, "Delta", "delta", "RDelta")),
        "beta": float(_param(params, "beta", "Beta")),
        "nmax": int(defaults.get("nmax", 2)),
        "ncut": defaults.get("ncut", 6),
        "basis_cap": int(defaults.get("basis_cap", 200000)),
    }
    if payload["ncut"] is not None:
        payload["ncut"] = int(payload["ncut"])
    (run_path / "params.json").write_text(
        json.dumps(payload, ensure_ascii=False, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )
    return run_path


def write_stage_report(output_dir: str | Path, manifest: Mapping[str, Any]) -> Path:
    """Write a compact Chinese HTML stage report for initialized inputs."""

    output_path = Path(output_dir)
    output_path.mkdir(parents=True, exist_ok=True)
    path = output_path / "report_stage_pairing_delta_inputs.html"
    cases = list(iter_cases(manifest))
    rows = "\n".join(
        "<tr>"
        f"<td>{html.escape(case_name(case))}</td>"
        f"<td>{float(_case_parameters(case)['Delta']):.6g}</td>"
        f"<td>{float(_case_parameters(case)['mu']):.6g}</td>"
        f"<td>{float(_case_parameters(case)['beta']):.6g}</td>"
        "</tr>"
        for case in cases
    )
    path.write_text(
        f"""<!doctype html>
<html lang="zh-CN">
<head>
  <meta charset="utf-8">
  <title>Pairing Delta campaign inputs</title>
  <script>
    window.MathJax = {{tex: {{inlineMath: [['$', '$'], ['\\\\(', '\\\\)']]}}}};
  </script>
  <script defer src="https://cdn.jsdelivr.net/npm/mathjax@3/es5/tex-chtml.js"></script>
  <style>
    body {{ font-family: system-ui, sans-serif; line-height: 1.55; margin: 2rem; }}
    table {{ border-collapse: collapse; width: 100%; }}
    th, td {{ border: 1px solid #bbb; padding: 0.35rem 0.5rem; text-align: left; }}
    code {{ background: #f5f5f5; padding: 0.1rem 0.25rem; }}
  </style>
</head>
<body>
  <h1>三角格 pairing $\\Delta$ 扫描输入</h1>
  <p>本阶段生成固定文件名 DQMC 输入和 ED 参数，用于后续 $3\\times3$ 试运行。</p>
  <p>案例数：{len(cases)}；campaign：<code>{html.escape(str(manifest.get("campaign", CAMPAIGN_NAME)))}</code></p>
  <table>
    <thead><tr><th>case</th><th>Delta</th><th>mu</th><th>beta</th></tr></thead>
    <tbody>{rows}</tbody>
  </table>
</body>
</html>
""",
        encoding="utf-8",
    )
    return path


def write_results_summary(
    input_dir: str | Path,
    output_dir: str | Path,
    manifest: Mapping[str, Any],
    *,
    block_size: int | None = None,
    skip_samples: int | None = None,
    report_name: str = "report_stage_results.html",
) -> dict[str, Path]:
    """Write minimal JSON/CSV/HTML summaries for all manifest cases."""

    analysis = _analysis_module()
    input_path = Path(input_dir)
    output_path = Path(output_dir)
    output_path.mkdir(parents=True, exist_ok=True)
    defaults = manifest.get("analysis_defaults", {})
    effective_block_size = int(block_size if block_size is not None else defaults.get("block_size", 10000))
    effective_skip = int(skip_samples if skip_samples is not None else defaults.get("skip_samples", 0))

    records: list[dict[str, Any]] = []
    for case in iter_cases(manifest):
        params = _case_parameters(case)
        run_dir = input_path / case_name(case)
        record: dict[str, Any] = {
            "case": case_name(case),
            "run_dir": str(run_dir),
            "Delta": float(_param(params, "Delta", "delta", "RDelta")),
            "status": "ok",
        }
        try:
            summary = analysis.summarize_dqmc_run(
                run_dir,
                Lx=int(params["Lx"]),
                Ly=int(params["Ly"]),
                block_size=effective_block_size,
                skip_samples=effective_skip,
            )
            record["reliable"] = bool(summary.get("reliable", True))
            if not record["reliable"]:
                issues = [str(issue) for issue in summary.get("issues", []) if str(issue)]
                record["status"] = "unreliable"
                record["reason"] = "; ".join(issues) if issues else str(
                    summary.get("reason", "dqmc_summary_unreliable")
                )
            for name, observable in summary.get("observables", {}).items():
                if isinstance(observable, Mapping):
                    record[f"{name}_mean"] = _float_or_blank(observable.get("actual", observable.get("value")))
                    record[f"{name}_sem"] = _float_or_blank(observable.get("stderr", observable.get("sem")))
                    if "imag_actual" in observable or "imag_value" in observable:
                        record[f"{name}_imag_mean"] = _float_or_blank(
                            observable.get("imag_actual", observable.get("imag_value"))
                        )
                        record[f"{name}_imag_sem"] = _float_or_blank(
                            observable.get("imag_stderr", observable.get("imag_sem"))
                        )
                        record[f"{name}_imag_threshold"] = _float_or_blank(
                            observable.get("imag_threshold")
                        )
                        record[f"{name}_pass_imaginary"] = bool(
                            observable.get("pass_imaginary", False)
                        )
                        record[f"{name}_imaginary_diagnostic_only"] = bool(
                            observable.get("imaginary_diagnostic_only", False)
                        )
        except (FileNotFoundError, ValueError) as exc:
            record.update({"status": "unavailable", "reliable": False, "reason": str(exc)})
        records.append(record)

    json_path = output_path / "summary_pairing_delta_results.json"
    csv_path = output_path / "summary_pairing_delta_results.csv"
    html_path = output_path / report_name
    payload = {
        "campaign": manifest.get("campaign", CAMPAIGN_NAME),
        "input_dir": str(input_path),
        "block_size": effective_block_size,
        "skip_samples": effective_skip,
        "cases": records,
        "counts": _status_counts(records),
    }
    json_path.write_text(
        json.dumps(payload, ensure_ascii=False, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )
    _write_csv(csv_path, records)
    html_path.write_text(_results_html(payload), encoding="utf-8")
    return {"json": json_path, "csv": csv_path, "html": html_path}


def write_manifest(path: str | Path, manifest: Mapping[str, Any] | None = None) -> Path:
    output = Path(path)
    output.parent.mkdir(parents=True, exist_ok=True)
    data = default_manifest() if manifest is None else dict(manifest)
    output.write_text(
        json.dumps(data, ensure_ascii=False, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )
    return output


def load_manifest(path: str | Path | None) -> dict[str, Any]:
    if path is None:
        return default_manifest()
    payload = json.loads(Path(path).read_text(encoding="utf-8"))
    if not isinstance(payload, dict):
        raise ValueError(f"manifest must be a JSON object: {path}")
    return payload


def init_local(
    manifest: Mapping[str, Any],
    output_dir: str | Path,
    *,
    seed: int | None = None,
) -> list[Path]:
    output_path = Path(output_dir)
    output_path.mkdir(parents=True, exist_ok=True)
    paths: list[Path] = []
    for index, case in enumerate(iter_cases(manifest)):
        run_dir = output_path / case_name(case)
        write_dqmc_case(
            case,
            run_dir,
            manifest.get("dqmc_defaults", {}),
            seed=None if seed is None else int(seed) + index,
        )
        write_ed_case(case, run_dir, manifest.get("ed_defaults", {}))
        paths.append(run_dir)
    return paths


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    subparsers = parser.add_subparsers(dest="command", required=True)

    write_manifest_parser = subparsers.add_parser("write-manifest")
    write_manifest_parser.add_argument("--output", type=Path, required=True)
    _add_manifest_override_args(write_manifest_parser)

    init_parser = subparsers.add_parser("init-local")
    init_parser.add_argument("--manifest", type=Path)
    init_parser.add_argument("--output-dir", type=Path, required=True)
    init_parser.add_argument("--seed", type=int)

    stage_parser = subparsers.add_parser("stage-report")
    stage_parser.add_argument("--manifest", type=Path)
    stage_parser.add_argument("--output-dir", type=Path, required=True)

    summary_parser = subparsers.add_parser("summarize-results")
    summary_parser.add_argument("--manifest", type=Path)
    summary_parser.add_argument("--input-dir", type=Path, required=True)
    summary_parser.add_argument("--output-dir", type=Path, required=True)
    summary_parser.add_argument("--block-size", type=int)
    summary_parser.add_argument("--skip-samples", type=int)
    summary_parser.add_argument("--report-name", default="report_stage_results.html")

    args = parser.parse_args(argv)
    try:
        if args.command == "write-manifest":
            write_manifest(args.output, _manifest_from_args(args))
            return 0
        if args.command == "init-local":
            init_local(load_manifest(args.manifest), args.output_dir, seed=args.seed)
            return 0
        if args.command == "stage-report":
            write_stage_report(args.output_dir, load_manifest(args.manifest))
            return 0
        if args.command == "summarize-results":
            write_results_summary(
                args.input_dir,
                args.output_dir,
                load_manifest(args.manifest),
                block_size=args.block_size,
                skip_samples=args.skip_samples,
                report_name=args.report_name,
            )
            return 0
    except ValueError as exc:
        parser.error(str(exc))
    raise AssertionError(f"unhandled command: {args.command}")


def _add_manifest_override_args(parser: argparse.ArgumentParser) -> None:
    parser.add_argument("--Lx", type=int, default=3)
    parser.add_argument("--Ly", type=int, default=3)
    parser.add_argument("--U1", type=float, default=1.0)
    parser.add_argument("--U2", type=float, default=0.0)
    parser.add_argument("--beta", type=float, default=4.0)
    parser.add_argument("--mu", type=float, default=-3.5)
    parser.add_argument("--deltas", default=",".join(str(value) for value in DEFAULT_DELTAS))
    parser.add_argument("--Nbin", type=int, default=100000)
    parser.add_argument("--block-size", type=int, default=10000)
    parser.add_argument("--nmax", type=int, default=2)
    parser.add_argument("--ncut", default="6")
    parser.add_argument("--basis-cap", type=int, default=200000)


def _manifest_from_args(args: argparse.Namespace) -> dict[str, Any]:
    return default_manifest(
        Lx=args.Lx,
        Ly=args.Ly,
        U1=args.U1,
        U2=args.U2,
        beta=args.beta,
        mu=args.mu,
        deltas=_parse_float_list(args.deltas),
        Nbin=args.Nbin,
        block_size=args.block_size,
        nmax=args.nmax,
        ncut=_parse_optional_int(args.ncut),
        basis_cap=args.basis_cap,
    )


def _parse_float_list(text: str) -> list[float]:
    values = [part.strip() for part in text.split(",") if part.strip()]
    if not values:
        raise ValueError("--deltas must contain at least one value")
    return [float(value) for value in values]


def _parse_optional_int(text: str) -> int | None:
    if text.lower() in {"none", "null"}:
        return None
    return int(text)


def _case_parameters(case: Mapping[str, Any]) -> Mapping[str, Any]:
    for key in ("parameters", "params", "dqmc_parameters"):
        value = case.get(key)
        if isinstance(value, Mapping):
            return value
    return case


def _param(params: Mapping[str, Any], *names: str) -> Any:
    for name in names:
        if name in params:
            return params[name]
    raise KeyError(f"missing any of parameters: {', '.join(names)}")


def _merged_defaults(
    case: Mapping[str, Any],
    key: str,
    explicit: Mapping[str, Any] | None,
) -> dict[str, Any]:
    if key == "dqmc":
        baseline = default_manifest()["dqmc_defaults"]
    else:
        baseline = default_manifest()["ed_defaults"]
    defaults = dict(baseline)
    if explicit is not None:
        defaults.update(explicit)
    nested = case.get(key)
    if isinstance(nested, Mapping):
        defaults.update(nested)
    return defaults


def _paramc_text(params: Mapping[str, Any], defaults: Mapping[str, Any]) -> str:
    lx = int(_param(params, "Lx", "Nlx"))
    ly = int(_param(params, "Ly", "Nly"))
    beta = float(_param(params, "beta", "Beta"))
    ltrot = int(defaults.get("Ltrot", round(beta / float(defaults.get("dtau", 0.01)))))
    if ltrot <= 0:
        raise ValueError("Ltrot must be positive")
    rows = [
        [
            float(_param(params, "U1")), float(_param(params, "U2")),
            float(_param(params, "mu", "Mu")),
            float(_param(params, "Delta", "delta", "RDelta")),
        ],
        [lx, ly, ltrot, beta],
        [
            int(defaults.get("NlxTherm", lx)),
            int(defaults.get("NlyTherm", ly)),
            int(defaults.get("LtrotTherm", ltrot)),
        ],
        [
            int(defaults.get("Nwrap", 10)),
            int(defaults.get("Nbin", 100000)),
            int(defaults.get("Nsweep", 1)),
            float(defaults.get("shiftLoc", 1.0)),
        ],
        [_fortran_bool(bool(defaults.get("is_tau", False))), int(defaults.get("Nthermal", 0))],
        [
            _fortran_bool(bool(defaults.get("is_warm", True))),
            int(defaults.get("Nwarm", 500)),
            float(defaults.get("shiftWarm1", 1.0)),
            float(defaults.get("shiftWarm2", 1.0)),
        ],
        [
            int(defaults.get("iniType", 2)),
            float(defaults.get("iniAmpl", 0.1)),
            float(defaults.get("iniBias1", 0.0)),
            float(defaults.get("iniBias2", 0.0)),
        ],
    ]
    if len(rows) != PARAMC_DATA_ROWS:
        raise AssertionError("paramC numeric block row count changed")
    lines = [" ".join(_format_param_token(value) for value in row) for row in rows]
    lines.extend(
        [
            "",
            "U1          U2          mu          RDelta",
            "Nlx         Nly         Ltrot       Beta",
            "NlxTherm    NlyTherm    LtrotTherm",
            "Nwrap       Nbin        Nsweep      shiftLoc",
            "is_tau      Nthermal",
            "is_warm     Nwarm       shiftWarm1  shiftWarm2",
            "iniType     iniAmpl     iniBias1    iniBias2",
        ]
    )
    return "\n".join(lines) + "\n"


def _format_param_token(value: Any) -> str:
    if isinstance(value, str):
        return value
    if isinstance(value, bool):
        return _fortran_bool(value)
    if isinstance(value, int):
        return str(value)
    number = float(value)
    if not math.isfinite(number):
        raise ValueError(f"non-finite paramC value: {value!r}")
    return f"{number:.12g}"


def _fortran_bool(value: bool) -> str:
    return ".true." if value else ".false."


def _seeds_text(seed: int, count: int = 128) -> str:
    value = int(seed) % 100000000
    if value <= 0:
        value += 1000003
    seeds: list[str] = []
    for _ in range(count):
        value = (1664525 * value + 1013904223) % 100000000
        if value <= 0:
            value += 1
        seeds.append(str(value))
    return "\n".join(seeds) + "\n"


def _format_label_float(value: float) -> str:
    text = f"{value:.12g}"
    return text.replace("+", "").replace(".", "p").replace("-", "m")


def _analysis_module():
    try:
        from . import pairing_delta_analysis as analysis
    except ModuleNotFoundError:
        import pairing_delta_analysis as analysis  # type: ignore
    return analysis


def _float_or_blank(value: Any) -> float | str:
    try:
        number = float(value)
    except (TypeError, ValueError):
        return ""
    return number if math.isfinite(number) else ""


def _status_counts(records: Sequence[Mapping[str, Any]]) -> dict[str, int]:
    counts: dict[str, int] = {}
    for record in records:
        status = str(record.get("status", "unknown"))
        counts[status] = counts.get(status, 0) + 1
    return counts


def _write_csv(path: Path, records: Sequence[Mapping[str, Any]]) -> None:
    keys: list[str] = []
    for record in records:
        for key in record:
            if key not in keys:
                keys.append(key)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=keys)
        writer.writeheader()
        for record in records:
            writer.writerow(record)


def _results_html(payload: Mapping[str, Any]) -> str:
    columns = [
        ("case", "case"),
        ("Delta", "Delta"),
        ("status", "status"),
        ("reason", "reason"),
        ("density_total_mean", "rho"),
        ("energy_density_mean", "e=E/Ns"),
        ("S_SF_K_mean", "Re S_SF_K"),
        ("S_SF_K_imag_mean", "Im S_SF_K"),
        ("S_SF_K_pass_imaginary", "Im pass"),
        ("S_DW_K_mean", "Re S_DW_K"),
        ("S_DW_K_imag_mean", "Im S_DW_K"),
        ("S_DW_K_pass_imaginary", "Im pass"),
        ("S_PSF_Gamma_mean", "Re S_PSF_Gamma"),
        ("S_PSF_Gamma_imag_mean", "Im S_PSF_Gamma"),
        ("S_PSF_Gamma_pass_imaginary", "Im pass"),
    ]
    def row_html(record: Mapping[str, Any]) -> str:
        cells = "".join(
            f"<td>{html.escape(str(record.get(key, 'ok' if key == 'reason' else '')))}</td>"
            for key, _ in columns
        )
        return f"<tr>{cells}</tr>"

    rows = "\n".join(row_html(record) for record in payload.get("cases", []))
    header = "".join(f"<th>{html.escape(label)}</th>" for _, label in columns)
    return f"""<!doctype html>
<html lang="zh-CN">
<head>
  <meta charset="utf-8">
  <title>Pairing Delta results summary</title>
  <style>
    body {{ font-family: system-ui, sans-serif; line-height: 1.55; margin: 2rem; }}
    table {{ border-collapse: collapse; width: 100%; }}
    th, td {{ border: 1px solid #bbb; padding: 0.35rem 0.5rem; text-align: left; }}
    .table-wrap {{ overflow-x: auto; }}
  </style>
</head>
<body>
  <h1>Pairing Delta 结果摘要</h1>
  <p>这是自动生成的阶段性检查表；最终物理报告应根据实际曲线手写。`S_SF_K` 和 `S_DW_K` 显示单 K 点复数估计器的实部与虚部诊断，比较 ED 时使用实部。</p>
  <div class="table-wrap">
    <table>
      <thead><tr>{header}</tr></thead>
      <tbody>{rows}</tbody>
    </table>
  </div>
</body>
</html>
"""


if __name__ == "__main__":
    raise SystemExit(main())
