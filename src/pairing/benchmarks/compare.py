#!/usr/bin/env python3
"""Compare DQMC scalar output files against benchmark references."""

from __future__ import annotations

import argparse
import json
import math
from pathlib import Path
from typing import Any, Iterable


REPO_ROOT = Path(__file__).resolve().parents[1]
SUPPORTED_REFERENCE_KINDS = {
    "analytic",
    "ed_json",
    "nopairing_dqmc_json",
    "live_nopairing_dqmc",
}


def resolve_repo_path(path: str | Path) -> Path:
    candidate = Path(path)
    if candidate.is_absolute():
        return candidate
    return REPO_ROOT / candidate


def read_series(path: Path) -> list[float]:
    if not path.exists():
        raise FileNotFoundError(f"missing DQMC output file: {path}")

    values: list[float] = []
    for line in path.read_text(encoding="utf-8").splitlines():
        stripped = line.strip()
        if not stripped:
            continue
        values.append(float(stripped.split()[0]))

    if not values:
        raise ValueError(f"no numeric values found in {path}")
    return values


def series_values(run_dir: Path, files: Iterable[str]) -> list[list[float]]:
    return [read_series(run_dir / name) for name in files]


def ensure_min_samples(
    series: list[list[float]], files: list[str], min_samples: int
) -> None:
    for name, values in zip(files, series):
        if len(values) < min_samples:
            raise ValueError(
                f"{name} has {len(values)} samples, fewer than required {min_samples}"
            )


def ensure_equal_lengths(series: list[list[float]], files: list[str]) -> None:
    lengths = {len(values) for values in series}
    if len(lengths) != 1:
        rendered = ", ".join(
            f"{name}={len(values)}" for name, values in zip(files, series)
        )
        raise ValueError(f"observable files have different sample counts: {rendered}")


def mean(values: list[float]) -> float:
    return sum(values) / float(len(values))


def parameter_nsite(parameters: dict[str, Any]) -> int:
    if "Nsite" in parameters:
        return int(parameters["Nsite"])
    return int(parameters["Lx"]) * int(parameters["Ly"])


def observable_files(spec: dict[str, Any]) -> list[str]:
    files = spec["files"]
    if not isinstance(files, list) or not files:
        raise ValueError("dqmc.files must be a non-empty list")
    return [str(name) for name in files]


def operation_scale(operation: str, nsite: int) -> float:
    if operation in {"mean_times_nsite", "mean_times_lq"}:
        return float(nsite)
    return 1.0


def compute_observable_samples(
    run_dir: Path,
    dqmc: dict[str, Any],
    nsite: int,
) -> list[float]:
    operation = str(dqmc["operation"])
    files = observable_files(dqmc)
    min_samples = int(dqmc.get("min_samples", 1))
    series = series_values(run_dir, files)
    ensure_min_samples(series, files, min_samples)

    if operation == "mean":
        if len(series) != 1:
            raise ValueError("operation 'mean' requires exactly one file")
        return series[0]

    if operation == "sum_mean":
        ensure_equal_lengths(series, files)
        return [sum(values) for values in zip(*series)]

    if operation in {"mean_times_nsite", "mean_times_lq"}:
        if len(series) != 1:
            raise ValueError(f"operation '{operation}' requires exactly one file")
        scale = operation_scale(operation, nsite)
        return [value * scale for value in series[0]]

    if operation == "last":
        raise ValueError("operation 'last' does not expose per-sample statistics")

    raise ValueError(f"unsupported comparison operation: {operation}")


def compute_actual(run_dir: Path, dqmc: dict[str, Any], nsite: int) -> float:
    operation = str(dqmc["operation"])
    files = observable_files(dqmc)
    min_samples = int(dqmc.get("min_samples", 1))

    if operation == "last":
        series = series_values(run_dir, files)
        ensure_min_samples(series, files, min_samples)
        if len(series) != 1:
            raise ValueError("operation 'last' requires exactly one file")
        return series[0][-1]

    return mean(compute_observable_samples(run_dir, dqmc, nsite))


def block_statistics(
    samples: list[float], block_size: int, skip_samples: int = 0
) -> dict[str, float | int]:
    if block_size <= 0:
        raise ValueError("statistics.block_size must be positive")
    if skip_samples < 0:
        raise ValueError("statistics.skip_samples must be non-negative")

    trimmed = samples[skip_samples:]
    usable_count = (len(trimmed) // block_size) * block_size
    if usable_count < 2 * block_size:
        raise ValueError(
            f"need at least two full blocks; have {len(trimmed)} samples after "
            f"skip_samples={skip_samples}, block_size={block_size}"
        )

    usable = trimmed[:usable_count]
    blocks = [
        mean(usable[index : index + block_size])
        for index in range(0, usable_count, block_size)
    ]
    block_mean = mean(blocks)
    if len(blocks) < 2:
        raise ValueError("need at least two blocks to estimate statistical error")

    variance = sum((value - block_mean) ** 2 for value in blocks) / float(
        len(blocks) - 1
    )
    stderr = math.sqrt(variance / float(len(blocks)))
    return {
        "actual": block_mean,
        "stderr": stderr,
        "blocks": len(blocks),
        "samples_used": usable_count,
        "skip_samples": skip_samples,
    }


def compute_result(
    run_dir: Path, dqmc: dict[str, Any], nsite: int
) -> dict[str, float | int | None]:
    statistics = dqmc.get("statistics")
    if statistics:
        samples = compute_observable_samples(run_dir, dqmc, nsite)
        stats = block_statistics(
            samples,
            int(statistics["block_size"]),
            int(statistics.get("skip_samples", 0)),
        )
        return {
            "value": float(stats["actual"]),
            "stderr": float(stats["stderr"]),
            "blocks": int(stats["blocks"]),
            "samples_used": int(stats["samples_used"]),
            "skip_samples": int(stats["skip_samples"]),
        }

    return {
        "value": compute_actual(run_dir, dqmc, nsite),
        "stderr": None,
        "blocks": None,
        "samples_used": None,
        "skip_samples": None,
    }


def stderr_tolerance_for(
    spec: dict[str, Any],
    dqmc: dict[str, Any],
    fallback_dqmc: dict[str, Any] | None = None,
) -> float:
    if "stderr_tolerance" in spec:
        return float(spec["stderr_tolerance"])

    for candidate in (dqmc, fallback_dqmc):
        if candidate is None:
            continue
        statistics = candidate.get("statistics")
        if not statistics:
            continue
        if "stderr_tolerance" in statistics:
            return float(statistics["stderr_tolerance"])
        if "sigma_tolerance" in statistics:
            return float(statistics["sigma_tolerance"])

    return 3.0


def paired_stderr_tolerance(name: str, spec: dict[str, Any]) -> float:
    if "stderr_tolerance" not in spec:
        raise ValueError(
            f"{name} paired DQMC comparison requires observable-level "
            "stderr_tolerance"
        )

    for side in ("dqmc", "reference_dqmc"):
        statistics = spec[side].get("statistics")
        if not isinstance(statistics, dict) or "block_size" not in statistics:
            raise ValueError(
                f"{name} paired DQMC comparison requires "
                f"{side}.statistics.block_size"
            )

    return float(spec["stderr_tolerance"])


def dotted_value(payload: dict[str, Any], path: str) -> Any:
    current: Any = payload
    for part in path.split("."):
        if not isinstance(current, dict) or part not in current:
            raise KeyError(f"path '{path}' not found at component '{part}'")
        current = current[part]
    return current


def discover_references(reference: Path | None, reference_dir: Path | None) -> list[Path]:
    if reference is not None and reference_dir is not None:
        raise ValueError("use only one of --reference or --reference-dir")

    if reference_dir is not None:
        directory = resolve_repo_path(reference_dir)
        references = sorted(directory.glob("*.json"))
    elif reference is not None:
        reference_path = resolve_repo_path(reference)
        references = (
            sorted(reference_path.glob("*.json"))
            if reference_path.is_dir()
            else [reference_path]
        )
    else:
        references = sorted((REPO_ROOT / "benchmarks" / "references").glob("*.json"))

    if not references:
        raise FileNotFoundError("no benchmark reference JSON files found")
    return references


def resolve_run_dir(
    reference_path: Path, reference: dict[str, Any], explicit_run_dir: Path | None
) -> Path:
    if explicit_run_dir is not None:
        return resolve_repo_path(explicit_run_dir)

    fixture = reference.get("dqmc_fixture")
    if not fixture:
        raise KeyError(f"{reference_path} does not define dqmc_fixture")
    return resolve_repo_path(fixture)


def compare_one_sided(
    name: str,
    run_dir: Path,
    dqmc: dict[str, Any],
    nsite: int,
    expected: float,
    spec: dict[str, Any],
) -> bool:
    result = compute_result(run_dir, dqmc, nsite)
    actual = float(result["value"])
    atol = float(spec.get("atol", 0.0))
    rtol = float(spec.get("rtol", 0.0))
    stderr = result["stderr"]

    if stderr is None:
        tolerance = max(atol, rtol * abs(expected))
        passed = abs(actual - expected) <= tolerance
        status = "PASS" if passed else "FAIL"
        print(
            f"{status} {name}: actual={actual:.16g} expected={expected:.16g} "
            f"abs_diff={abs(actual - expected):.3g} atol={atol:.3g} rtol={rtol:.3g}"
        )
        return passed

    stderr_tolerance = stderr_tolerance_for(spec, dqmc)
    tolerance = max(atol, rtol * abs(expected), stderr_tolerance * float(stderr))
    abs_diff = abs(actual - expected)
    passed = abs_diff <= tolerance
    if float(stderr) == 0.0:
        z_score = 0.0 if abs_diff == 0.0 else math.inf
    else:
        z_score = (actual - expected) / float(stderr)
    status = "PASS" if passed else "FAIL"
    print(
        f"{status} {name}: actual={actual:.16g} expected={expected:.16g} "
        f"abs_diff={abs_diff:.3g} stderr={float(stderr):.3g} z={z_score:.3g} "
        f"blocks={int(result['blocks'])} samples_used={int(result['samples_used'])} "
        f"stderr_tolerance={stderr_tolerance:.3g} threshold={tolerance:.3g}"
    )
    return passed


def compare_paired_dqmc(
    name: str,
    run_dir: Path,
    reference_run_dir: Path,
    spec: dict[str, Any],
    nsite: int,
) -> bool:
    dqmc = spec["dqmc"]
    reference_dqmc = spec["reference_dqmc"]
    stderr_tolerance = paired_stderr_tolerance(name, spec)
    actual_result = compute_result(run_dir, dqmc, nsite)
    expected_result = compute_result(reference_run_dir, reference_dqmc, nsite)

    actual = float(actual_result["value"])
    expected = float(expected_result["value"])
    atol = float(spec.get("atol", 0.0))
    rtol = float(spec.get("rtol", 0.0))
    actual_stderr = actual_result["stderr"]
    expected_stderr = expected_result["stderr"]

    if actual_stderr is None or expected_stderr is None:
        raise ValueError(
            f"{name} paired DQMC comparison requires block statistics on both "
            "dqmc and reference_dqmc"
        )

    stderr_combined = math.sqrt(float(actual_stderr) ** 2 + float(expected_stderr) ** 2)
    tolerance = max(atol, rtol * abs(expected), stderr_tolerance * stderr_combined)
    abs_diff = abs(actual - expected)
    passed = abs_diff <= tolerance
    if stderr_combined == 0.0:
        z_score = 0.0 if abs_diff == 0.0 else math.inf
    else:
        z_score = (actual - expected) / stderr_combined

    status = "PASS" if passed else "FAIL"
    print(
        f"{status} {name}: actual={actual:.16g} expected={expected:.16g} "
        f"abs_diff={abs_diff:.3g} stderr_pairing={float(actual_stderr):.3g} "
        f"stderr_reference={float(expected_stderr):.3g} "
        f"stderr_combined={stderr_combined:.3g} z={z_score:.3g} "
        f"stderr_tolerance={stderr_tolerance:.3g} threshold={tolerance:.3g} "
        f"pairing_blocks={int(actual_result['blocks'])} "
        f"reference_blocks={int(expected_result['blocks'])}"
    )
    return passed


def run_case(reference_path: Path, explicit_run_dir: Path | None = None) -> list[str]:
    reference = json.loads(reference_path.read_text(encoding="utf-8"))
    reference_kind = str(reference.get("reference_kind", "analytic"))
    if reference_kind not in SUPPORTED_REFERENCE_KINDS:
        raise ValueError(f"unsupported reference_kind: {reference_kind}")
    if reference_kind == "live_nopairing_dqmc":
        raise NotImplementedError(
            "reference_kind 'live_nopairing_dqmc' is recognized but not implemented; "
            "use 'nopairing_dqmc_json' committed-output references until a live "
            "dual-repository runner is added"
        )

    params = reference["parameters"]
    nsite = parameter_nsite(params)
    run_dir = resolve_run_dir(reference_path, reference, explicit_run_dir)
    ed_payload: dict[str, Any] | None = None
    reference_run_dir: Path | None = None

    if reference_kind == "ed_json":
        ed_payload = json.loads(resolve_repo_path(reference["ed_json"]).read_text("utf-8"))
    elif reference_kind == "nopairing_dqmc_json":
        reference_run_dir = resolve_repo_path(reference["reference_run_dir"])

    failures: list[str] = []
    case_name = reference.get("case", reference_path.stem)
    print(f"Case: {case_name}")
    print(f"Reference: {reference_path}")
    print(f"Reference kind: {reference_kind}")
    print(f"DQMC run directory: {run_dir}")
    if reference_run_dir is not None:
        print(f"No-pairing reference run directory: {reference_run_dir}")
    print("")

    for name, spec in reference["observables"].items():
        if reference_kind == "nopairing_dqmc_json" and "reference_dqmc" in spec:
            if reference_run_dir is None:
                raise RuntimeError("reference_run_dir was not initialized")
            passed = compare_paired_dqmc(name, run_dir, reference_run_dir, spec, nsite)
        else:
            if reference_kind == "ed_json":
                if ed_payload is None:
                    raise RuntimeError("ed_json payload was not initialized")
                expected = float(dotted_value(ed_payload, spec["expected_from"]))
            else:
                expected = float(spec["value"])
            passed = compare_one_sided(
                name,
                run_dir,
                spec["dqmc"],
                nsite,
                expected,
                spec,
            )

        if not passed:
            failures.append(name)

    print("")
    return failures


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--reference",
        type=Path,
        help="Reference JSON file, or a directory containing reference JSON files",
    )
    parser.add_argument(
        "--reference-dir",
        type=Path,
        help="Directory containing reference JSON files",
    )
    parser.add_argument(
        "--run-dir",
        type=Path,
        help="Directory containing DQMC scalar output files for a single reference",
    )
    args = parser.parse_args()

    try:
        references = discover_references(args.reference, args.reference_dir)
        if args.run_dir is not None and len(references) != 1:
            raise ValueError("--run-dir can only be used with one reference file")

        case_failures: dict[str, list[str]] = {}
        for reference_path in references:
            failures = run_case(reference_path, args.run_dir)
            if failures:
                case_failures[reference_path.stem] = failures
    except Exception as exc:
        print(f"ERROR: {exc}")
        return 2

    if case_failures:
        print("Failed benchmark cases:")
        for case_name, failures in case_failures.items():
            print(f"- {case_name}: {', '.join(failures)}")
        return 1

    print("All benchmark cases passed.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
