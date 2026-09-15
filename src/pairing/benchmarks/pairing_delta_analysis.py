#!/usr/bin/env python3
"""Fail-closed helpers for triangular pairing Delta-sweep DQMC output."""

from __future__ import annotations

import json
import math
from pathlib import Path
from typing import Any, Iterable, Sequence


DEFAULT_REQUIRED_SCALAR_FILES = (
    "density_total",
    "energy_density",
    "interaction_energy_density",
    "pairing_energy_density",
    "chemical_energy_density",
    "grand_energy_density",
    "onsite_n2_up",
    "onsite_n2_do",
)
OPTIONAL_SCALAR_FILES = (
    "density",
    "density_up",
    "density_do",
    "num_up",
    "num_do",
    "kinetic",
    "doubleOcc",
    "squareOcc",
    "local_numsquare",
    "numsquare_up",
    "numsquare_do",
    "pair_equal",
)
COMPLEX_FILES = ("sf_K", "dw_K", "psf_Gamma")
STRUCTURE_FACTOR_NAMES = {"sf_K", "dw_K"}
OBSERVABLE_ALIASES = {
    "sf_K": "S_SF_K",
    "dw_K": "S_DW_K",
    "psf_Gamma": "S_PSF_Gamma",
    "IPR": "IPR",
}
ED_TARGET_OBSERVABLES = (
    "density_total",
    "energy_density",
    "doubleOcc",
    "interaction_energy_density",
    "pairing_energy_density",
    "chemical_energy_density",
    "grand_energy_density",
    "onsite_n2_up",
    "onsite_n2_do",
    "S_SF_K",
    "S_DW_K",
    "S_PSF_Gamma",
)
IMAGINARY_ROUNDOFF_ATOL = 1.0e-14


def read_scalar_series(path: str | Path) -> list[float]:
    """Read a scalar output file with exactly one finite numeric column."""

    values: list[float] = []
    path = Path(path)
    for line_number, line in _numeric_lines(path):
        parts = line.split()
        if len(parts) != 1:
            raise ValueError(
                f"malformed scalar row in {path} at line {line_number}: "
                f"expected 1 column, found {len(parts)}"
            )
        values.append(_parse_finite_float(parts[0], path, line_number, "scalar"))
    if not values:
        raise ValueError(f"no scalar samples found in {path}")
    return values


def read_complex_series(path: str | Path) -> list[complex]:
    """Read a complex output file with exactly two finite numeric columns."""

    values: list[complex] = []
    path = Path(path)
    for line_number, line in _numeric_lines(path):
        parts = line.split()
        if len(parts) != 2:
            raise ValueError(
                f"malformed complex row in {path} at line {line_number}: "
                f"expected 2 columns, found {len(parts)}"
            )
        real = _parse_finite_float(parts[0], path, line_number, "complex real")
        imag = _parse_finite_float(parts[1], path, line_number, "complex imaginary")
        values.append(complex(real, imag))
    if not values:
        raise ValueError(f"no complex samples found in {path}")
    return values


def read_vector_series(
    path: str | Path,
    Lq: int | None = None,
    *,
    expected_length: int | None = None,
) -> list[list[float]]:
    """Read fixed-width site-vector rows with exactly ``Lq`` finite values."""

    width = expected_length if expected_length is not None else Lq
    if width is None:
        raise ValueError("read_vector_series requires Lq or expected_length")
    width = int(width)
    if width <= 0:
        raise ValueError("vector row length must be positive")

    rows: list[list[float]] = []
    path = Path(path)
    for line_number, line in _numeric_lines(path):
        parts = line.split()
        if len(parts) != width:
            raise ValueError(
                f"malformed density_site_total row in {path} at line {line_number}: "
                f"expected {width} columns, found {len(parts)}"
            )
        rows.append(
            [_parse_finite_float(token, path, line_number, "vector") for token in parts]
        )
    if not rows:
        raise ValueError(f"no vector samples found in {path}")
    return rows


def block_statistics(
    samples: Sequence[float],
    block_size: int,
    skip_samples: int = 0,
) -> dict[str, float | int | bool | str]:
    """Return mean and block SEM from complete block means."""

    values = [_finite_number(value, "sample") for value in samples]
    block_means, samples_used = _block_means(values, block_size, skip_samples)
    return _statistics_from_block_values(
        block_means,
        samples_used=samples_used,
        skip_samples=skip_samples,
        block_size=block_size,
    )


def complex_block_statistics(
    samples: Sequence[complex | Sequence[float]],
    block_size: int,
    skip_samples: int = 0,
) -> dict[str, Any]:
    """Return independent block statistics for real and imaginary columns."""

    real_values, imag_values = _complex_components(samples)
    real_stats = block_statistics(real_values, block_size, skip_samples)
    imag_stats = block_statistics(imag_values, block_size, skip_samples)
    return {
        "real": real_stats,
        "real_stats": real_stats,
        "imag": imag_stats,
        "imaginary": imag_stats,
        "imag_stats": imag_stats,
        "reliable": True,
        "is_reliable": True,
        "available": True,
        "valid": True,
        "status": "ok",
        "reason": "ok",
    }


def compute_ipr(
    site_density: Iterable[float],
    min_total_density: float = 1.0e-8,
) -> dict[str, float | bool | str]:
    """Compute ``sum_i rho_i^2 / (sum_i rho_i)^2`` from a site-density row."""

    try:
        values = [_finite_number(value, "site density") for value in site_density]
    except ValueError as exc:
        return _ipr_result(math.nan, False, str(exc))
    if not values:
        return _ipr_result(math.nan, False, "empty_density_profile")

    total = sum(values)
    if total <= float(min_total_density):
        return _ipr_result(math.nan, False, "density_total_below_threshold")

    ipr = sum(value * value for value in values) / (total * total)
    return _ipr_result(ipr, True, "ok")


def summarize_dqmc_run(
    run_dir: str | Path,
    *,
    Lx: int | None = None,
    Ly: int | None = None,
    Lq: int | None = None,
    block_size: int = 1000,
    skip_samples: int = 0,
    required_files: Sequence[str] | None = None,
    scalar_files: Sequence[str] | None = None,
    min_ipr_total_density: float = 1.0e-8,
) -> dict[str, Any]:
    """Summarize one DQMC run directory using strict row validation."""

    run_path = Path(run_dir)
    if Lq is None and Lx is not None and Ly is not None:
        Lq = int(Lx) * int(Ly)
    if Lq is not None and int(Lq) <= 0:
        raise ValueError("Lq must be positive")

    if required_files is not None:
        summary = _summarize_explicit_required(
            run_path,
            required_files,
            Lq=Lq,
            block_size=block_size,
            skip_samples=skip_samples,
        )
        summary.update({"run_dir": str(run_path), "Lx": Lx, "Ly": Ly, "Lq": Lq})
        return summary

    if Lq is None:
        Lq = _infer_vector_width(run_path / "density_site_total")

    requested_scalars = list(
        scalar_files if scalar_files is not None else DEFAULT_REQUIRED_SCALAR_FILES
    )
    if scalar_files is None:
        requested_scalars.extend(
            name for name in OPTIONAL_SCALAR_FILES if (run_path / name).exists()
        )

    scalar_series: dict[str, list[float]] = {
        name: read_scalar_series(run_path / name) for name in requested_scalars
    }
    density_site_total = read_vector_series(run_path / "density_site_total", int(Lq))

    compatible_k = _k_geometry_compatible(Lx, Ly)
    complex_series: dict[str, list[complex]] = {}
    for name in COMPLEX_FILES:
        if name in STRUCTURE_FACTOR_NAMES and compatible_k is False:
            continue
        complex_series[name] = read_complex_series(run_path / name)

    _ensure_equal_lengths(
        {**scalar_series, **complex_series, "density_site_total": density_site_total}
    )

    summary: dict[str, Any] = {
        "run_dir": str(run_path),
        "Lx": Lx,
        "Ly": Ly,
        "Lq": int(Lq),
        "block_size": int(block_size),
        "skip_samples": int(skip_samples),
        "reliable": True,
        "is_reliable": True,
        "ok": True,
        "success": True,
        "status": "ok",
        "issues": [],
        "observables": {},
    }

    for name, samples in scalar_series.items():
        stats = block_statistics(samples, block_size, skip_samples)
        if name == "density_total" and float(stats["actual"]) <= 0.0:
            _mark_observable_unreliable(stats, "non_positive_density_total")
        summary[name] = stats
        summary["observables"][name] = dict(stats)
        if not bool(stats.get("reliable", True)):
            _mark_summary_unreliable(summary, f"{name}:{stats.get('reason', 'unreliable')}")

    for name in COMPLEX_FILES:
        if name in STRUCTURE_FACTOR_NAMES and compatible_k is False:
            unavailable = _unavailable_observable(
                "K-point observable requires Lx and Ly to be multiples of 3"
            )
            summary[name] = unavailable
            summary["observables"][OBSERVABLE_ALIASES[name]] = dict(unavailable)
            _mark_summary_unreliable(summary, f"{name}:{unavailable['reason']}")
            continue
        stats = complex_block_statistics(complex_series[name], block_size, skip_samples)
        observable = _complex_observable(
            stats,
            positive_semidefinite=name in STRUCTURE_FACTOR_NAMES,
        )
        if not bool(observable.get("reliable", True)):
            _mark_observable_unreliable(stats, str(observable.get("reason", "complex_unreliable")))
        summary[name] = stats
        summary["observables"][OBSERVABLE_ALIASES[name]] = observable
        if not bool(observable.get("reliable", True)):
            _mark_summary_unreliable(summary, f"{name}:{observable.get('reason', 'unreliable')}")

    ipr_stats = _ipr_summary(
        density_site_total,
        block_size=block_size,
        skip_samples=skip_samples,
        min_total_density=min_ipr_total_density,
    )
    summary["IPR"] = ipr_stats
    summary["ipr"] = ipr_stats
    summary["observables"]["IPR"] = dict(ipr_stats)
    if not bool(ipr_stats.get("reliable", True)):
        _mark_summary_unreliable(summary, f"IPR:{ipr_stats.get('reason', 'unreliable')}")
    return summary


def load_ed_result(path: str | Path) -> dict[str, Any]:
    payload = json.loads(Path(path).read_text(encoding="utf-8"))
    if not isinstance(payload, dict):
        raise ValueError(f"ED result must be a JSON object: {path}")
    return payload


def ed_target_observables(result: dict[str, Any]) -> dict[str, dict[str, Any]]:
    """Return ED observables in the DQMC comparison shape when fields exist."""

    observables = result.get("observables", result)
    if not isinstance(observables, dict):
        raise ValueError("ED result is missing an observables object")

    targets: dict[str, dict[str, Any]] = {}
    for name in ED_TARGET_OBSERVABLES:
        if name in observables:
            value = _observable_value(observables[name], f"ED observable {name}")
            targets[name] = _ed_observable(value)
    return targets


def compare_dqmc_ed_case(
    run_dir: str | Path,
    *,
    block_size: int,
    Lx: int | None = None,
    Ly: int | None = None,
    Lq: int | None = None,
    ed_result_path: str | Path | None = None,
    stderr_tolerance: float = 3.0,
    atol: float = 0.0,
    rtol: float = 0.0,
    skip_samples: int = 0,
) -> dict[str, Any]:
    """Compare one DQMC case against an ED JSON result if available."""

    run_path = Path(run_dir)
    result_path = Path(ed_result_path) if ed_result_path is not None else run_path / "results.json"
    payload: dict[str, Any] = {
        "run_dir": str(run_path),
        "ed_result_path": str(result_path),
        "observables": {},
    }

    try:
        dqmc = summarize_dqmc_run(
            run_path,
            Lx=Lx,
            Ly=Ly,
            Lq=Lq,
            block_size=block_size,
            skip_samples=skip_samples,
        )
    except (FileNotFoundError, ValueError) as exc:
        payload.update({"trusted": False, "status": "dqmc_unavailable", "reason": str(exc)})
        return payload
    payload["dqmc"] = dqmc

    try:
        ed_result = load_ed_result(result_path)
        ed_reliable, ed_reason = _ed_payload_reliability(ed_result)
        ed = ed_target_observables(ed_result)
    except (FileNotFoundError, ValueError, json.JSONDecodeError) as exc:
        payload.update({"trusted": False, "status": "ed_unavailable", "reason": str(exc)})
        return payload

    trusted = True
    reasons: list[str] = []
    for name, dqmc_observable in dqmc["observables"].items():
        if name not in ed:
            continue
        ed_observable = dict(ed[name])
        if not ed_reliable:
            ed_observable["reliable"] = False
            ed_observable["is_reliable"] = False
            ed_observable["reason"] = ed_reason
        comparison = _compare_observable(
            dqmc_observable,
            ed_observable,
            stderr_tolerance=stderr_tolerance,
            atol=atol,
            rtol=rtol,
        )
        payload["observables"][name] = comparison
        if not comparison["reliable"]:
            trusted = False
            reasons.append(f"{name}:{comparison['reason']}")

    payload.update(
        {
            "trusted": trusted,
            "status": "ok" if trusted else "comparison_failed",
            "reason": "ok" if trusted else "; ".join(reasons),
        }
    )
    return payload


def _numeric_lines(path: Path) -> list[tuple[int, str]]:
    if not path.exists():
        raise FileNotFoundError(f"missing DQMC output file: {path}")
    return [
        (line_number, stripped)
        for line_number, line in enumerate(path.read_text(encoding="utf-8").splitlines(), 1)
        if (stripped := line.strip())
    ]


def _parse_finite_float(token: str, path: Path, line_number: int, value_name: str) -> float:
    try:
        value = float(token.replace("D", "E").replace("d", "e"))
    except ValueError as exc:
        raise ValueError(
            f"malformed {value_name} value in {path} at line {line_number}: {token!r}"
        ) from exc
    if not math.isfinite(value):
        raise ValueError(
            f"non-finite {value_name} value in {path} at line {line_number}: {token!r}"
        )
    return value


def _finite_number(value: Any, name: str) -> float:
    try:
        number = float(value)
    except (TypeError, ValueError) as exc:
        raise ValueError(f"{name} is not numeric: {value!r}") from exc
    if not math.isfinite(number):
        raise ValueError(f"{name} is not finite: {value!r}")
    return number


def _mark_observable_unreliable(observable: dict[str, Any], reason: str) -> None:
    observable["reliable"] = False
    observable["is_reliable"] = False
    observable["available"] = False
    observable["valid"] = False
    observable["status"] = "unreliable"
    observable["reason"] = reason


def _mark_summary_unreliable(summary: dict[str, Any], reason: str) -> None:
    summary["reliable"] = False
    summary["is_reliable"] = False
    summary["ok"] = False
    summary["success"] = False
    summary["status"] = "unreliable"
    issues = summary.setdefault("issues", [])
    if isinstance(issues, list):
        issues.append(reason)


def _observable_value(value: Any, name: str) -> float:
    if isinstance(value, dict):
        for key in ("value", "actual", "mean"):
            if key in value:
                return _finite_number(value[key], name)
        raise ValueError(f"{name} is missing value/actual/mean field: {value!r}")
    return _finite_number(value, name)


def _mean(values: Sequence[float]) -> float:
    if not values:
        raise ValueError("cannot compute mean of empty sequence")
    return sum(values) / float(len(values))


def _block_means(
    samples: Sequence[float],
    block_size: int,
    skip_samples: int,
) -> tuple[list[float], int]:
    block_size = int(block_size)
    skip_samples = int(skip_samples)
    if block_size <= 0:
        raise ValueError("block_size must be positive")
    if skip_samples < 0:
        raise ValueError("skip_samples must be non-negative")

    trimmed = list(samples)[skip_samples:]
    samples_used = (len(trimmed) // block_size) * block_size
    if samples_used < 2 * block_size:
        raise ValueError(
            f"need at least two complete blocks; have {len(trimmed)} samples "
            f"after skip_samples={skip_samples}, block_size={block_size}"
        )
    usable = trimmed[:samples_used]
    return (
        [
            _mean(usable[index : index + block_size])
            for index in range(0, samples_used, block_size)
        ],
        samples_used,
    )


def _statistics_from_block_values(
    block_values: Sequence[float],
    *,
    samples_used: int,
    skip_samples: int,
    block_size: int,
) -> dict[str, float | int | bool | str]:
    values = [_finite_number(value, "block mean") for value in block_values]
    if len(values) < 2:
        raise ValueError("need at least two block means")
    average = _mean(values)
    variance = sum((value - average) ** 2 for value in values) / float(len(values) - 1)
    sem = math.sqrt(variance / float(len(values)))
    return {
        "mean": average,
        "actual": average,
        "value": average,
        "sem": sem,
        "stderr": sem,
        "blocks": len(values),
        "n_blocks": len(values),
        "samples_used": int(samples_used),
        "skip_samples": int(skip_samples),
        "block_size": int(block_size),
        "reliable": True,
        "is_reliable": True,
        "available": True,
        "valid": True,
        "status": "ok",
        "reason": "ok",
    }


def _complex_components(
    samples: Sequence[complex | Sequence[float]],
) -> tuple[list[float], list[float]]:
    real_values: list[float] = []
    imag_values: list[float] = []
    for index, value in enumerate(samples):
        if isinstance(value, complex):
            real = value.real
            imag = value.imag
        else:
            try:
                if len(value) != 2:  # type: ignore[arg-type]
                    raise ValueError
                real = value[0]  # type: ignore[index]
                imag = value[1]  # type: ignore[index]
            except (TypeError, ValueError, IndexError) as exc:
                raise ValueError(f"complex sample {index} must have two columns") from exc
        real_values.append(_finite_number(real, "complex real sample"))
        imag_values.append(_finite_number(imag, "complex imaginary sample"))
    if not real_values:
        raise ValueError("cannot summarize empty complex series")
    return real_values, imag_values


def _infer_vector_width(path: Path) -> int:
    lines = _numeric_lines(path)
    first_line_number, first_line = lines[0]
    parts = first_line.split()
    if not parts:
        raise ValueError(f"empty vector row in {path} at line {first_line_number}")
    return len(parts)


def _ensure_equal_lengths(series: dict[str, Sequence[Any]]) -> None:
    lengths = {name: len(values) for name, values in series.items()}
    if not lengths:
        return
    expected = next(iter(lengths.values()))
    mismatched = {name: length for name, length in lengths.items() if length != expected}
    if mismatched:
        raise ValueError(f"DQMC output sample counts differ: {mismatched}; expected {expected}")


def _k_geometry_compatible(Lx: int | None, Ly: int | None) -> bool | None:
    if Lx is None or Ly is None:
        return None
    return int(Lx) % 3 == 0 and int(Ly) % 3 == 0


def _unavailable_observable(reason: str) -> dict[str, Any]:
    return {
        "mean": math.nan,
        "actual": math.nan,
        "value": math.nan,
        "sem": math.nan,
        "stderr": math.nan,
        "reliable": False,
        "is_reliable": False,
        "available": False,
        "valid": False,
        "status": "unavailable",
        "reason": reason,
    }


def _complex_observable(
    stats: dict[str, Any],
    *,
    imaginary_diagnostic_only: bool = False,
    positive_semidefinite: bool = False,
) -> dict[str, Any]:
    real = stats["real"]
    imag = stats["imag"]
    real_value = float(real["actual"])
    real_stderr = float(real["stderr"])
    imag_value = float(imag["actual"])
    imag_stderr = float(imag["stderr"])
    imag_threshold = max(IMAGINARY_ROUNDOFF_ATOL, 3.0 * imag_stderr)
    pass_imaginary = math.isfinite(imag_value) and abs(imag_value) <= imag_threshold
    nonnegative_threshold = max(IMAGINARY_ROUNDOFF_ATOL, 3.0 * real_stderr)
    pass_nonnegative = (
        not positive_semidefinite
        or (math.isfinite(real_value) and real_value >= -nonnegative_threshold)
    )
    reliable = pass_imaginary and pass_nonnegative
    reason = "ok"
    if not pass_imaginary:
        reason = "imaginary_part_not_zero"
    elif not pass_nonnegative:
        reason = "negative_structure_factor"
    return {
        "mean": real["mean"],
        "actual": real["actual"],
        "value": real["value"],
        "sem": real["sem"],
        "stderr": real["stderr"],
        "imag_mean": imag["mean"],
        "imag_actual": imag["actual"],
        "imag_value": imag["value"],
        "imag_sem": imag["sem"],
        "imag_stderr": imag["stderr"],
        "imag_threshold": imag_threshold,
        "pass_imaginary": pass_imaginary,
        "nonnegative_threshold": nonnegative_threshold,
        "positive_semidefinite": bool(positive_semidefinite),
        "pass_nonnegative": pass_nonnegative,
        "imaginary_diagnostic_only": bool(imaginary_diagnostic_only),
        "blocks": real["blocks"],
        "n_blocks": real["n_blocks"],
        "samples_used": real["samples_used"],
        "reliable": reliable,
        "is_reliable": reliable,
        "available": reliable,
        "valid": reliable,
        "status": "ok" if reliable else "unreliable",
        "reason": reason,
    }


def _ipr_summary(
    rows: Sequence[Sequence[float]],
    *,
    block_size: int,
    skip_samples: int,
    min_total_density: float,
) -> dict[str, Any]:
    block_size = int(block_size)
    skip_samples = int(skip_samples)
    if block_size <= 0:
        raise ValueError("block_size must be positive")
    trimmed = list(rows)[skip_samples:]
    samples_used = (len(trimmed) // block_size) * block_size
    if samples_used < 2 * block_size:
        raise ValueError(
            f"need at least two complete density_site_total blocks; have "
            f"{len(trimmed)} samples after skip_samples={skip_samples}, "
            f"block_size={block_size}"
        )
    usable = trimmed[:samples_used]
    width = len(usable[0])
    block_ipr_values: list[float] = []
    reliability_reasons: list[str] = []
    for index in range(0, samples_used, block_size):
        block = usable[index : index + block_size]
        block_mean = [_mean([row[column] for row in block]) for column in range(width)]
        ipr = compute_ipr(block_mean, min_total_density)
        if not bool(ipr["reliable"]):
            reliability_reasons.append(str(ipr["reason"]))
        else:
            block_ipr_values.append(float(ipr["value"]))

    if reliability_reasons:
        result = _unavailable_observable("; ".join(sorted(set(reliability_reasons))))
        result.update({"ipr": math.nan, "is_reliable": False})
        return result

    stats = _statistics_from_block_values(
        block_ipr_values,
        samples_used=samples_used,
        skip_samples=skip_samples,
        block_size=block_size,
    )
    stats.update({"ipr": stats["value"], "is_reliable": True})
    return stats


def _ipr_result(value: float, reliable: bool, reason: str) -> dict[str, float | bool | str]:
    return {
        "ipr": value,
        "value": value,
        "reliable": reliable,
        "is_reliable": reliable,
        "available": reliable,
        "valid": reliable,
        "status": "ok" if reliable else "unreliable",
        "reason": reason,
    }


def _summarize_explicit_required(
    run_path: Path,
    required_files: Sequence[str],
    *,
    Lq: int | None,
    block_size: int,
    skip_samples: int,
) -> dict[str, Any]:
    for name in required_files:
        path = run_path / name
        if not path.exists():
            raise FileNotFoundError(f"missing required DQMC output file: {path}")

    summary: dict[str, Any] = {
        "block_size": int(block_size),
        "skip_samples": int(skip_samples),
        "reliable": True,
        "is_reliable": True,
        "ok": True,
        "success": True,
        "status": "ok",
        "issues": [],
        "observables": {},
    }
    for name in required_files:
        path = run_path / name
        if name == "density_site_total":
            width = int(Lq) if Lq is not None else _infer_vector_width(path)
            rows = read_vector_series(path, width)
            summary[name] = {"rows": len(rows), "columns": width, "reliable": True}
            continue
        if name in COMPLEX_FILES:
            stats = complex_block_statistics(
                read_complex_series(path),
                block_size=block_size,
                skip_samples=skip_samples,
            )
            observable = _complex_observable(
                stats,
                positive_semidefinite=name in STRUCTURE_FACTOR_NAMES,
            )
            if not bool(observable.get("reliable", True)):
                _mark_observable_unreliable(stats, str(observable.get("reason", "complex_unreliable")))
                _mark_summary_unreliable(summary, f"{name}:{observable.get('reason', 'unreliable')}")
            summary[name] = stats
            summary["observables"][OBSERVABLE_ALIASES[name]] = observable
            continue
        stats = block_statistics(
            read_scalar_series(path),
            block_size=block_size,
            skip_samples=skip_samples,
        )
        summary[name] = stats
        summary["observables"][name] = dict(stats)
    return summary


def _ed_observable(value: float) -> dict[str, Any]:
    return {
        "actual": value,
        "value": value,
        "stderr": 0.0,
        "sem": 0.0,
        "reliable": True,
        "is_reliable": True,
        "reason": "ed_reference",
    }


def _ed_payload_reliability(result: dict[str, Any]) -> tuple[bool, str]:
    if result.get("reliable") is False or result.get("is_reliable") is False:
        return False, "ed_reference_unreliable"
    status = str(result.get("status", "ok")).lower()
    if any(word in status for word in ("unreliable", "failed", "blocked", "dry_run")):
        return False, f"ed_status_{status}"
    return True, "ok"


def _compare_observable(
    dqmc: dict[str, Any],
    ed: dict[str, Any],
    *,
    stderr_tolerance: float,
    atol: float,
    rtol: float,
) -> dict[str, Any]:
    dqmc_value = _finite_or_nan(dqmc.get("actual", dqmc.get("value")))
    ed_value = _finite_or_nan(ed.get("actual", ed.get("value")))
    stderr = _finite_or_nan(dqmc.get("stderr", dqmc.get("sem")))
    reliable = bool(dqmc.get("reliable", True)) and bool(ed.get("reliable", True))
    diff = dqmc_value - ed_value
    z_score = diff / stderr if math.isfinite(diff) and math.isfinite(stderr) and stderr > 0.0 else math.nan
    threshold = max(float(atol), abs(ed_value) * float(rtol), float(stderr_tolerance) * stderr)
    pass_uncertainty = math.isfinite(diff) and math.isfinite(threshold) and abs(diff) <= threshold
    imag = _imaginary_check(dqmc, stderr_tolerance=stderr_tolerance, atol=atol)
    imag_gate = bool(imag["pass_imaginary"])
    final_reliable = reliable and pass_uncertainty and imag_gate
    reason = "ok"
    if not reliable:
        reason = str(dqmc.get("reason", ed.get("reason", "unreliable")))
    elif not pass_uncertainty:
        reason = "outside_dqmc_uncertainty"
    elif not imag_gate:
        reason = "imaginary_part_not_zero"
    return {
        "dqmc": dqmc_value,
        "ed": ed_value,
        "difference": diff,
        "stderr": stderr,
        "z_score": z_score,
        "threshold": threshold,
        "pass_uncertainty": pass_uncertainty,
        **imag,
        "imaginary_diagnostic_only": bool(dqmc.get("imaginary_diagnostic_only", False)),
        "reliable": final_reliable,
        "reason": reason,
    }


def _imaginary_check(
    observable: dict[str, Any],
    *,
    stderr_tolerance: float,
    atol: float,
) -> dict[str, Any]:
    if "imag_actual" not in observable and "imag_value" not in observable:
        return {
            "imag_actual": math.nan,
            "imag_stderr": math.nan,
            "imag_threshold": math.nan,
            "pass_imaginary": True,
        }
    imag = _finite_or_nan(observable.get("imag_actual", observable.get("imag_value")))
    imag_sem = _finite_or_nan(observable.get("imag_stderr", observable.get("imag_sem")))
    threshold = max(float(atol), IMAGINARY_ROUNDOFF_ATOL, float(stderr_tolerance) * imag_sem)
    return {
        "imag_actual": imag,
        "imag_stderr": imag_sem,
        "imag_threshold": threshold,
        "pass_imaginary": math.isfinite(imag) and abs(imag) <= threshold,
    }


def _finite_or_nan(value: Any) -> float:
    try:
        number = float(value)
    except (TypeError, ValueError):
        return math.nan
    return number if math.isfinite(number) else math.nan
