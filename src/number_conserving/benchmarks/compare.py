#!/usr/bin/env python3
"""Compare DQMC scalar output files against documented ED reference values."""

from __future__ import annotations

import argparse
import json
import math
from pathlib import Path
from typing import Any, Iterable


REPO_ROOT = Path(__file__).resolve().parents[1]


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


def ensure_min_samples(series: list[list[float]], files: list[str], min_samples: int) -> None:
    for name, values in zip(files, series):
        if len(values) < min_samples:
            raise ValueError(
                f"{name} has {len(values)} samples, fewer than required {min_samples}"
            )


def mean(values: list[float]) -> float:
    return sum(values) / float(len(values))


def ensure_equal_lengths(series: list[list[float]], files: list[str]) -> None:
    lengths = {len(values) for values in series}
    if len(lengths) != 1:
        rendered = ", ".join(
            f"{name}={len(values)}" for name, values in zip(files, series)
        )
        raise ValueError(f"observable files have different sample counts: {rendered}")


def compute_observable_samples(
    run_dir: Path,
    operation: str,
    files: list[str],
    lq: int,
    min_samples: int = 1,
) -> list[float]:
    series = series_values(run_dir, files)
    ensure_min_samples(series, files, min_samples)

    if operation == "mean":
        if len(series) != 1:
            raise ValueError("operation 'mean' requires exactly one file")
        return series[0]
    if operation == "sum_mean":
        ensure_equal_lengths(series, files)
        return [sum(values) for values in zip(*series)]
    if operation == "mean_times_lq":
        if len(series) != 1:
            raise ValueError("operation 'mean_times_lq' requires exactly one file")
        return [value * float(lq) for value in series[0]]
    if operation == "mean_times_lq_over_2":
        if len(series) != 1:
            raise ValueError("operation 'mean_times_lq_over_2' requires exactly one file")
        return [value * float(lq) / 2.0 for value in series[0]]

    raise ValueError(
        f"operation '{operation}' does not expose per-sample values for statistics"
    )


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


def compute_actual(
    run_dir: Path,
    operation: str,
    files: list[str],
    lq: int,
    min_samples: int = 1,
) -> float:
    series = series_values(run_dir, files)
    ensure_min_samples(series, files, min_samples)
    last_values = [values[-1] for values in series]
    mean_values = [mean(values) for values in series]
    if operation == "last":
        if len(last_values) != 1:
            raise ValueError("operation 'last' requires exactly one file")
        return last_values[0]
    if operation == "sum_last":
        return sum(last_values)
    if operation == "last_times_lq":
        if len(last_values) != 1:
            raise ValueError("operation 'last_times_lq' requires exactly one file")
        return last_values[0] * float(lq)
    if operation == "last_times_lq_over_2":
        if len(last_values) != 1:
            raise ValueError("operation 'last_times_lq_over_2' requires exactly one file")
        return last_values[0] * float(lq) / 2.0
    if operation == "mean":
        return mean(compute_observable_samples(run_dir, operation, files, lq, min_samples))
    if operation == "sum_mean":
        return mean(compute_observable_samples(run_dir, operation, files, lq, min_samples))
    if operation == "mean_times_lq":
        return mean(compute_observable_samples(run_dir, operation, files, lq, min_samples))
    if operation == "mean_times_lq_over_2":
        return mean(compute_observable_samples(run_dir, operation, files, lq, min_samples))
    raise ValueError(f"unsupported comparison operation: {operation}")


def discover_references(reference: Path | None, reference_dir: Path | None) -> list[Path]:
    if reference is not None and reference_dir is not None:
        raise ValueError("use only one of --reference or --reference-dir")
    if reference_dir is not None:
        references = sorted(reference_dir.glob("*.json"))
    elif reference is not None:
        references = sorted(reference.glob("*.json")) if reference.is_dir() else [reference]
    else:
        references = sorted((REPO_ROOT / "benchmarks" / "references").glob("*.json"))

    if not references:
        raise FileNotFoundError("no benchmark reference JSON files found")
    return references


def resolve_run_dir(
    reference_path: Path, reference: dict[str, Any], explicit_run_dir: Path | None
) -> Path:
    if explicit_run_dir is not None:
        return explicit_run_dir

    fixture = reference.get("dqmc_fixture")
    if not fixture:
        raise KeyError(f"{reference_path} does not define dqmc_fixture")

    fixture_path = Path(fixture)
    if fixture_path.is_absolute():
        return fixture_path
    return REPO_ROOT / fixture_path


def run_case(reference_path: Path, explicit_run_dir: Path | None = None) -> list[str]:
    reference = json.loads(reference_path.read_text(encoding="utf-8"))
    params = reference["parameters"]
    lq = int(params["Lx"]) * int(params["Ly"])
    run_dir = resolve_run_dir(reference_path, reference, explicit_run_dir)

    failures: list[str] = []
    case_name = reference.get("case", reference_path.stem)
    print(f"Case: {case_name}")
    print(f"Reference: {reference_path}")
    print(f"DQMC run directory: {run_dir}")
    print("")

    for name, spec in reference["observables"].items():
        dqmc = spec["dqmc"]
        operation = dqmc["operation"]
        files = dqmc["files"]
        min_samples = int(dqmc.get("min_samples", 1))
        statistics = dqmc.get("statistics")
        expected = float(spec["value"])
        atol = float(spec.get("atol", 0.0))
        rtol = float(spec.get("rtol", 0.0))

        if statistics:
            samples = compute_observable_samples(
                run_dir, operation, files, lq, min_samples
            )
            stats = block_statistics(
                samples,
                int(statistics["block_size"]),
                int(statistics.get("skip_samples", 0)),
            )
            actual = float(stats["actual"])
            stderr = float(stats["stderr"])
            stderr_tolerance = float(
                statistics.get(
                    "stderr_tolerance", statistics.get("sigma_tolerance", 3.0)
                )
            )
            tolerance = max(atol, abs(expected) * rtol, stderr_tolerance * stderr)
            abs_diff = abs(actual - expected)
            passed = abs_diff <= tolerance
            if stderr == 0.0:
                z_score = 0.0 if abs_diff == 0.0 else math.inf
            else:
                z_score = (actual - expected) / stderr
            status = "PASS" if passed else "FAIL"
            print(
                f"{status} {name}: actual={actual:.16g} expected={expected:.16g} "
                f"abs_diff={abs_diff:.3g} stderr={stderr:.3g} z={z_score:.3g} "
                f"blocks={stats['blocks']} samples_used={stats['samples_used']} "
                f"stderr_tolerance={stderr_tolerance:.3g} threshold={tolerance:.3g}"
            )
        else:
            actual = compute_actual(run_dir, operation, files, lq, min_samples)
            passed = math.isclose(actual, expected, abs_tol=atol, rel_tol=rtol)
            status = "PASS" if passed else "FAIL"
            print(
                f"{status} {name}: actual={actual:.16g} expected={expected:.16g} "
                f"abs_diff={abs(actual - expected):.3g} atol={atol:.3g} rtol={rtol:.3g}"
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
