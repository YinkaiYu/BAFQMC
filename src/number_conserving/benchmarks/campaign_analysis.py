#!/usr/bin/env python3
"""Helpers for 3x3 BAFQMC-vs-ED campaign analysis."""

from __future__ import annotations

from dataclasses import dataclass
import json
import math
from pathlib import Path
from typing import Any, Iterable


TARGET_OBSERVABLES = {
    "density_total",
    "energy_density",
    "doubleOcc",
    "IPR",
    "S_SF_K",
    "S_PSF_Gamma",
    "S_DW_K",
}
STRUCTURE_FACTOR_OBSERVABLES = {
    "S_SF_K",
    "S_PSF_Gamma",
    "S_DW_K",
}
EXPECTED_ED_PARAMETER_FIELDS = (
    "Lx",
    "Ly",
    "beta",
    "mu",
    "U1",
    "U2",
)

ED_REFERENCE_OBSERVABLES = {
    "density_total",
    "total_kinetic",
    "doubleOcc",
    "onsite_n2_up",
    "onsite_n2_do",
    "interaction_energy_density",
    "energy_density",
    "S_SF_K",
    "S_PSF_Gamma",
    "S_DW_K",
    "num_up",
    "num_do",
}

NEGATIVE_U2_TRUSTED_LOW_DENSITY_MAX_TOTAL = 5.0e-2
NEGATIVE_U2_TRUSTED_LOW_DENSITY_MIN_SHELLS = 6
IMAGINARY_ROUNDOFF_ATOL = 1.0e-14


@dataclass(frozen=True)
class RatioResult:
    value: float
    reliable: bool
    reason: str


@dataclass(frozen=True)
class ReliabilityDecision:
    reliable: bool
    kind: str
    reason: str


def _numeric_lines(path: Path) -> list[tuple[int, str]]:
    if not path.exists():
        raise FileNotFoundError(f"missing DQMC output file: {path}")

    text = path.read_text(encoding="utf-8")
    return [
        (line_number, stripped)
        for line_number, line in enumerate(text.splitlines(), start=1)
        if (stripped := line.strip())
    ]


def _parse_finite_float(
    token: str, path: Path, line_number: int, value_name: str
) -> float:
    try:
        value = float(token)
    except ValueError as exc:
        raise ValueError(
            f"malformed {value_name} value in {path} at line {line_number}: {token!r}"
        ) from exc
    if not math.isfinite(value):
        raise ValueError(
            f"non-finite {value_name} value in {path} at line {line_number}: {token!r}"
        )
    return value


def read_scalar_series(path: Path) -> list[float]:
    """Read the first numeric column from a scalar DQMC output file."""
    values: list[float] = []
    for line_number, line in _numeric_lines(path):
        token = line.split()[0]
        values.append(_parse_finite_float(token, path, line_number, "scalar"))

    if not values:
        raise ValueError(f"no numeric values found in {path}")
    return values


def read_complex_series(path: Path) -> list[complex]:
    """Read Re/Im numeric columns from a complex DQMC output file."""
    values: list[complex] = []
    for line_number, line in _numeric_lines(path):
        parts = line.split()
        if len(parts) < 2:
            raise ValueError(
                f"malformed complex value in {path} at line {line_number}: {line!r}"
            )
        real = _parse_finite_float(parts[0], path, line_number, "complex real")
        imag = _parse_finite_float(parts[1], path, line_number, "complex imaginary")
        values.append(complex(real, imag))

    if not values:
        raise ValueError(f"no numeric values found in {path}")
    return values


def read_vector_series(path: Path, expected_length: int) -> list[list[float]]:
    """Read fixed-length numeric rows from a vector DQMC output file."""

    if expected_length <= 0:
        raise ValueError("expected_length must be positive")

    rows: list[list[float]] = []
    for line_number, line in _numeric_lines(path):
        parts = line.split()
        if len(parts) != expected_length:
            raise ValueError(
                f"expected {expected_length} columns in {path} at line "
                f"{line_number}; found {len(parts)}"
            )
        rows.append(
            [
                _parse_finite_float(token, path, line_number, "vector")
                for token in parts
            ]
        )

    if not rows:
        raise ValueError(f"no numeric values found in {path}")
    return rows


def mean(values: Iterable[float]) -> float:
    total = 0.0
    count = 0
    for value in values:
        total += value
        count += 1
    if count == 0:
        raise ValueError("cannot compute mean of empty values")
    return total / float(count)


def _block_means(
    samples: list[float], block_size: int, skip_samples: int = 0
) -> tuple[list[float], int]:
    if block_size <= 0:
        raise ValueError("block_size must be positive")
    if skip_samples < 0:
        raise ValueError("skip_samples must be non-negative")

    trimmed = samples[skip_samples:]
    usable_count = (len(trimmed) // block_size) * block_size
    if usable_count < 2 * block_size:
        raise ValueError(
            f"need at least two full blocks; have {len(trimmed)} samples after "
            f"skip_samples={skip_samples}, block_size={block_size}"
        )

    usable = trimmed[:usable_count]
    return [
        mean(usable[index : index + block_size])
        for index in range(0, usable_count, block_size)
    ], usable_count


def _vector_block_means(
    samples: list[list[float]], block_size: int, skip_samples: int = 0
) -> tuple[list[list[float]], int]:
    if block_size <= 0:
        raise ValueError("block_size must be positive")
    if skip_samples < 0:
        raise ValueError("skip_samples must be non-negative")

    trimmed = samples[skip_samples:]
    usable_count = (len(trimmed) // block_size) * block_size
    if usable_count < 2 * block_size:
        raise ValueError(
            f"need at least two full blocks; have {len(trimmed)} samples after "
            f"skip_samples={skip_samples}, block_size={block_size}"
        )

    usable = trimmed[:usable_count]
    width = len(usable[0])
    blocks: list[list[float]] = []
    for index in range(0, usable_count, block_size):
        block = usable[index : index + block_size]
        blocks.append(
            [mean(row[column] for row in block) for column in range(width)]
        )
    return blocks, usable_count


def block_statistics(
    samples: list[float], block_size: int, skip_samples: int = 0
) -> dict[str, float | int]:
    block_means, usable_count = _block_means(samples, block_size, skip_samples)
    actual = mean(block_means)
    variance = sum((value - actual) ** 2 for value in block_means) / float(
        len(block_means) - 1
    )

    return {
        "actual": actual,
        "stderr": math.sqrt(variance / float(len(block_means))),
        "blocks": len(block_means),
        "samples_used": usable_count,
        "skip_samples": skip_samples,
    }


def complex_block_statistics(
    samples: list[complex], block_size: int, skip_samples: int = 0
) -> dict[str, dict[str, float | int]]:
    return {
        "real": block_statistics(
            [value.real for value in samples],
            block_size=block_size,
            skip_samples=skip_samples,
        ),
        "imag": block_statistics(
            [value.imag for value in samples],
            block_size=block_size,
            skip_samples=skip_samples,
        ),
    }


def k_point_cell_index(lx: int, ly: int) -> tuple[int, int]:
    if lx % 3 != 0 or ly % 3 != 0:
        raise ValueError("K-point index requires lattice dimensions to be multiples of 3")
    return (2 * lx // 3 + 1, ly // 3 + 1)


def compute_ipr(site_density: Iterable[float], min_number: float = 1.0e-8) -> RatioResult:
    """Compute the normalized density-distribution IPR.

    The input is a site-resolved total-density profile rho_i. The returned
    value is sum_i rho_i^2 / (sum_i rho_i)^2, so a uniform 3x3 profile gives
    1/9 and a profile localized on one site gives 1.
    """

    values = [float(value) for value in site_density]
    if not values:
        return RatioResult(value=math.nan, reliable=False, reason="empty_density_profile")

    total_density = sum(values)
    if abs(total_density) < min_number:
        return RatioResult(
            value=math.nan,
            reliable=False,
            reason="density_below_minimum",
        )

    value = sum(value * value for value in values) / (total_density * total_density)
    return RatioResult(value=value, reliable=True, reason="ok")


def load_ed_result(path: Path) -> dict[str, Any]:
    """Load an ED result JSON payload from a production case directory."""

    if not path.exists():
        raise FileNotFoundError(f"missing ED result file: {path}")
    payload = json.loads(path.read_text(encoding="utf-8"))
    if not isinstance(payload, dict):
        raise ValueError(f"ED result payload must be a JSON object: {path}")
    return payload


def ed_target_observables(
    result: dict[str, Any],
    min_ipr_number: float = 1.0e-8,
) -> dict[str, dict[str, float | bool | str]]:
    """Return ED values on the same target-observable names as DQMC summaries."""

    observables = result.get("observables", {})
    if not isinstance(observables, dict):
        raise ValueError("ED result is missing an observables object")

    targets: dict[str, dict[str, float | bool | str]] = {}
    for name in (
        "density_total",
        "energy_density",
        "doubleOcc",
        "S_SF_K",
        "S_PSF_Gamma",
        "S_DW_K",
    ):
        targets[name] = {
            "actual": _finite_float(observables[name], f"ED observable {name}"),
            "reliable": True,
            "reason": "direct_ed_observable",
        }

    parameters = result.get("parameters", {})
    lq = int(parameters["Lx"]) * int(parameters["Ly"])
    density_total = _finite_float(
        observables["density_total"], "ED observable density_total"
    )
    ipr = compute_ipr([density_total] * lq, min_number=min_ipr_number)
    targets["IPR"] = {
        "actual": ipr.value,
        "reliable": ipr.reliable,
        "reason": ipr.reason,
    }
    return targets


def free_boson_reference_observables(parameters: dict[str, Any]) -> dict[str, float]:
    """Return exact grand-canonical observables for the stable U1=U2=0 case."""

    lx = int(parameters["Lx"])
    ly = int(parameters["Ly"])
    lq = lx * ly
    beta = float(parameters["beta"])
    mu = float(parameters["mu"])
    t = float(parameters.get("t", 1.0))
    u1 = float(parameters.get("U1", 0.0))
    u2 = float(parameters.get("U2", 0.0))
    if not math.isclose(u1, 0.0, rel_tol=0.0, abs_tol=1.0e-14):
        raise ValueError("free-boson reference requires U1=U2=0")
    if not math.isclose(u2, 0.0, rel_tol=0.0, abs_tol=1.0e-14):
        raise ValueError("free-boson reference requires U2=0")
    if beta <= 0.0:
        raise ValueError("free-boson reference requires beta > 0")

    occupations: dict[tuple[int, int], float] = {}
    energies: dict[tuple[int, int], float] = {}
    for kx in range(lx):
        for ky in range(ly):
            energy = _triangular_single_particle_energy(lx, ly, kx, ky, t)
            gap = energy - mu
            if gap <= 0.0:
                raise ValueError(
                    "free-boson reference is unstable because min(epsilon)-mu <= 0"
                )
            energies[(kx, ky)] = energy
            occupations[(kx, ky)] = 1.0 / math.expm1(beta * gap)

    total_per_flavor = sum(occupations.values())
    total_variance_per_flavor = sum(
        occupation * (1.0 + occupation) for occupation in occupations.values()
    )
    total_square_per_flavor = (
        total_per_flavor * total_per_flavor + total_variance_per_flavor
    )
    density_total = 2.0 * total_per_flavor / float(lq)
    total_kinetic = 2.0 * sum(
        energies[momentum] * occupation
        for momentum, occupation in occupations.items()
    )
    site_density_per_flavor = total_per_flavor / float(lq)
    k_point = (2 * lx // 3, ly // 3) if lx % 3 == 0 and ly % 3 == 0 else (0, 0)
    sf_k = 2.0 * occupations[k_point] / float(lq)
    psf_gamma = sum(
        occupation * occupations[((-kx) % lx, (-ky) % ly)]
        for (kx, ky), occupation in occupations.items()
    ) / float(lq * lq)
    dw_k = 2.0 * sum(
        occupations[((kx + k_point[0]) % lx, (ky + k_point[1]) % ly)]
        * (1.0 + occupation)
        for (kx, ky), occupation in occupations.items()
    ) / float(lq * lq)
    double_occ = site_density_per_flavor * site_density_per_flavor
    same_flavor_n2_sum = total_per_flavor + 2.0 * double_occ * float(lq)

    return {
        "density_total": density_total,
        "total_kinetic": total_kinetic,
        "doubleOcc": double_occ,
        "onsite_n2_up": same_flavor_n2_sum,
        "onsite_n2_do": same_flavor_n2_sum,
        "interaction_energy_density": 0.0,
        "energy_density": total_kinetic / float(lq),
        "S_SF_K": sf_k,
        "S_PSF_Gamma": psf_gamma,
        "S_DW_K": dw_k,
        "num_up": total_per_flavor,
        "num_do": total_per_flavor,
        "numsquare_up": total_square_per_flavor,
        "numsquare_do": total_square_per_flavor,
    }


def free_boson_reference_result(
    parameters: dict[str, Any],
    expected_policy: dict[str, Any] | None = None,
) -> dict[str, Any]:
    """Build an ED-shaped exact reference payload for stable U1=U2=0 cases."""

    reference_parameters = {
        "Lx": int(parameters["Lx"]),
        "Ly": int(parameters["Ly"]),
        "beta": float(parameters["beta"]),
        "mu": float(parameters["mu"]),
        "U1": float(parameters.get("U1", 0.0)),
        "U2": float(parameters.get("U2", 0.0)),
        "t": float(parameters.get("t", 1.0)),
    }
    if expected_policy is not None and "max_total_particles" in expected_policy:
        reference_parameters["max_total_particles"] = int(
            expected_policy["max_total_particles"]
        )
    observables = free_boson_reference_observables(reference_parameters)
    return {
        "schema_version": 1,
        "parameters": reference_parameters,
        "observables": observables,
        "last_shell_contribution": {name: 0.0 for name in ED_REFERENCE_OBSERVABLES},
        "relative_last_shell": {name: 0.0 for name in ED_REFERENCE_OBSERVABLES},
        "convergence_policy": {
            "reference": "exact_free_boson",
        },
        "status": "exact_free_boson",
        "completed_shells": 0,
    }


def ed_result_matches_case(
    result: dict[str, Any],
    expected_parameters: dict[str, Any],
    expected_policy: dict[str, Any] | None = None,
    *,
    abs_tol: float = 1.0e-10,
) -> ReliabilityDecision:
    """Check that an ED result JSON belongs to the requested manifest case."""

    parameters = result.get("parameters", {})
    if not isinstance(parameters, dict):
        return ReliabilityDecision(
            reliable=False,
            kind="unreliable",
            reason="missing_or_invalid_parameters",
        )

    for field in EXPECTED_ED_PARAMETER_FIELDS:
        if field not in expected_parameters:
            continue
        if field not in parameters:
            return ReliabilityDecision(
                reliable=False,
                kind="unreliable",
                reason=f"ed_parameter_missing: {field}",
            )
        if not _numbers_match(parameters[field], expected_parameters[field], abs_tol):
            return ReliabilityDecision(
                reliable=False,
                kind="unreliable",
                reason=(
                    f"ed_parameter_mismatch: {field}="
                    f"{parameters[field]!r} expected {expected_parameters[field]!r}"
                ),
            )

    if expected_policy is not None and "max_total_particles" in expected_policy:
        field = "max_total_particles"
        if field not in parameters:
            return ReliabilityDecision(
                reliable=False,
                kind="unreliable",
                reason=f"ed_parameter_missing: {field}",
            )
        if int(parameters[field]) != int(expected_policy[field]):
            return ReliabilityDecision(
                reliable=False,
                kind="unreliable",
                reason=(
                    f"ed_parameter_mismatch: {field}="
                    f"{parameters[field]!r} expected {expected_policy[field]!r}"
                ),
            )

    return ReliabilityDecision(
        reliable=True,
        kind="case_match",
        reason="ed_parameters_match_manifest_case",
    )


def compare_observable(
    dqmc_observable: dict[str, Any],
    ed_observable: dict[str, Any],
    *,
    stderr_tolerance: float = 3.0,
    atol: float = 0.0,
    rtol: float = 0.0,
) -> dict[str, float | bool | str]:
    """Compare one DQMC observable against an ED value using block SEM."""

    dqmc_reliable = bool(dqmc_observable.get("reliable", True))
    ed_reliable = bool(ed_observable.get("reliable", True))
    dqmc_value = _finite_or_none(dqmc_observable.get("actual"))
    ed_value = _finite_or_none(ed_observable.get("actual"))
    stderr = _finite_or_none(dqmc_observable.get("stderr"))
    if dqmc_value is None or ed_value is None or stderr is None or stderr < 0.0:
        return {
            **_comparison_base(
                dqmc_observable,
                ed_observable,
                stderr_tolerance=stderr_tolerance,
                atol=atol,
                rtol=rtol,
            ),
            "reason": _first_unreliable_reason(
                dqmc_observable,
                ed_observable,
                dqmc_reliable=dqmc_reliable,
                ed_reliable=ed_reliable,
                fallback="missing_or_invalid_comparison_value",
            ),
        }

    diff = dqmc_value - ed_value
    threshold = max(
        float(atol),
        abs(ed_value) * float(rtol),
        float(stderr_tolerance) * stderr,
    )
    if stderr > 0.0:
        z_score = diff / stderr
    elif diff == 0.0:
        z_score = 0.0
    else:
        z_score = math.copysign(math.inf, diff)

    imag_diagnostics = _imaginary_diagnostics(
        dqmc_observable,
        stderr_tolerance=stderr_tolerance,
        atol=atol,
    )
    pass_uncertainty = abs(diff) <= threshold
    pass_imaginary = bool(imag_diagnostics["pass_imaginary"])
    reliable = dqmc_reliable and ed_reliable and pass_uncertainty and pass_imaginary
    reason = "ok"
    if not dqmc_reliable:
        reason = str(dqmc_observable.get("reason", "dqmc_observable_unreliable"))
    elif not ed_reliable:
        reason = str(ed_observable.get("reason", "ed_observable_unreliable"))
    elif not pass_uncertainty:
        reason = "outside_dqmc_uncertainty"
    elif not pass_imaginary:
        reason = "imaginary_part_not_zero"

    return {
        "dqmc": dqmc_value,
        "ed": ed_value,
        "difference": diff,
        "stderr": stderr,
        "z_score": z_score,
        "threshold": threshold,
        "stderr_tolerance": float(stderr_tolerance),
        "atol": float(atol),
        "rtol": float(rtol),
        "pass_uncertainty": pass_uncertainty,
        **imag_diagnostics,
        "reliable": reliable,
        "reason": reason,
    }


def compare_dqmc_ed_case(
    run_dir: Path,
    *,
    block_size: int,
    lq: int,
    expected_parameters: dict[str, Any] | None = None,
    expected_policy: dict[str, Any] | None = None,
    ed_result_path: Path | None = None,
    skip_samples: int = 0,
    min_ipr_number: float = 1.0e-8,
    stderr_tolerance: float = 3.0,
    atol: float = 0.0,
    rtol: float = 0.0,
) -> dict[str, Any]:
    """Summarize one case and fail closed unless DQMC and ED are both accepted."""

    ed_path = ed_result_path if ed_result_path is not None else run_dir / "results.json"
    result: dict[str, Any] = {
        "run_dir": str(run_dir),
        "ed_result_path": str(ed_path),
        "observables": {},
    }

    dqmc_summary: dict[str, Any] | None = None
    dqmc_error: str | None = None
    try:
        dqmc_summary = summarize_dqmc_run(
            run_dir,
            block_size=block_size,
            lq=lq,
            skip_samples=skip_samples,
            min_ipr_number=min_ipr_number,
        )
    except (FileNotFoundError, ValueError) as exc:
        dqmc_error = str(exc)
    if dqmc_summary is not None:
        result["dqmc"] = dqmc_summary

    ed_result: dict[str, Any] | None = None
    ed_error: str | None = None
    try:
        ed_result = load_ed_result(ed_path)
    except (FileNotFoundError, ValueError, json.JSONDecodeError) as exc:
        ed_error = str(exc)

    if expected_parameters is not None and _is_free_boson_case(expected_parameters):
        use_exact_reference = ed_result is None
        if ed_result is not None:
            candidate_match = ed_result_matches_case(
                ed_result,
                expected_parameters,
                expected_policy,
            )
            candidate_decision = (
                ed_reference_reliable(ed_result)
                if candidate_match.reliable
                else candidate_match
            )
            use_exact_reference = not (
                candidate_match.reliable and candidate_decision.reliable
            )
        if use_exact_reference:
            try:
                ed_result = free_boson_reference_result(
                    expected_parameters,
                    expected_policy,
                )
                ed_error = None
                result["ed_result_path"] = "exact_free_boson"
            except (KeyError, ValueError) as exc:
                ed_error = str(exc)

    ed_targets: dict[str, dict[str, Any]] | None = None
    ed_target_error: str | None = None
    ed_match = ReliabilityDecision(True, "case_match", "not_checked")
    ed_decision = ReliabilityDecision(False, "unavailable", ed_error or "ed_unavailable")
    if ed_result is not None:
        result["ed"] = ed_result
        if expected_parameters is not None:
            ed_match = ed_result_matches_case(
                ed_result,
                expected_parameters,
                expected_policy,
            )
        ed_decision = ed_reference_reliable(ed_result) if ed_match.reliable else ed_match
        try:
            ed_targets = ed_target_observables(
                ed_result,
                min_ipr_number=min_ipr_number,
            )
        except (KeyError, ValueError) as exc:
            ed_target_error = str(exc)

    result["ed_reliability"] = {
        "reliable": ed_decision.reliable,
        "kind": ed_decision.kind,
        "reason": ed_decision.reason,
    }
    result["ed_case_match"] = {
        "reliable": ed_match.reliable,
        "kind": ed_match.kind,
        "reason": ed_match.reason,
    }

    if dqmc_summary is None or ed_targets is None:
        reason = (
            dqmc_error
            if dqmc_summary is None
            else ed_target_error or ed_error or "ed_unavailable"
        )
        for name in sorted(TARGET_OBSERVABLES):
            result["observables"][name] = compare_observable(
                _unavailable_observable(
                    None if dqmc_summary is None else dqmc_summary["observables"].get(name),
                    reliable=False,
                    reason=dqmc_error or "dqmc_unavailable",
                ),
                _unavailable_observable(
                    None if ed_targets is None else ed_targets.get(name),
                    reliable=False,
                    reason=ed_target_error or ed_error or "ed_unavailable",
                ),
                stderr_tolerance=stderr_tolerance,
                atol=atol,
                rtol=rtol,
            )
        result.update(
            {
                "trusted": False,
                "status": "dqmc_unavailable" if dqmc_summary is None else "ed_unavailable",
                "reason": reason or "unavailable",
            }
        )
        return result

    trusted = True
    reasons: list[str] = []
    for name in sorted(TARGET_OBSERVABLES):
        ed_observable = dict(ed_targets[name])
        if not ed_decision.reliable:
            ed_observable["reliable"] = False
            ed_observable["reason"] = ed_decision.reason
        dqmc_observable = _dqmc_observable_for_comparison(
            name, dqmc_summary["observables"][name]
        )
        comparison = compare_observable(
            dqmc_observable,
            ed_observable,
            stderr_tolerance=stderr_tolerance,
            atol=atol,
            rtol=rtol,
        )
        result["observables"][name] = comparison
        if not bool(comparison["reliable"]):
            trusted = False
            reasons.append(f"{name}:{comparison['reason']}")

    result.update(
        {
            "trusted": trusted,
            "status": _case_status(
                trusted=trusted,
                ed_match=ed_match,
                ed_decision=ed_decision,
            ),
            "reason": "ok" if trusted else "; ".join(reasons),
        }
    )
    return result


def _dqmc_observable_for_comparison(
    name: str, observable: dict[str, Any]
) -> dict[str, Any]:
    values = dict(observable)
    if name == "density_total":
        actual = _finite_or_none(values.get("actual"))
        if actual is not None and actual <= 0.0:
            values["reliable"] = False
            values["reason"] = "non_positive_density"
    return values


def ed_reference_reliable(result: dict[str, Any]) -> ReliabilityDecision:
    status = str(result.get("status", ""))
    parameters = result.get("parameters", {})
    if not isinstance(parameters, dict):
        return ReliabilityDecision(
            reliable=False,
            kind="unreliable",
            reason="missing_or_invalid_parameters",
        )
    if status == "dry_run" or bool(parameters.get("dry_run", False)):
        return ReliabilityDecision(
            reliable=False,
            kind="unreliable",
            reason="dry_run_payload",
        )
    if result.get("schema_version") != 1:
        return ReliabilityDecision(
            reliable=False,
            kind="unreliable",
            reason="unsupported_or_missing_schema_version",
        )
    if status == "exact_free_boson":
        if not _is_free_boson_case(parameters):
            return ReliabilityDecision(
                reliable=False,
                kind="unreliable",
                reason="exact_free_boson_parameter_mismatch",
            )
        observables = result.get("observables", {})
        if not isinstance(observables, dict):
            return ReliabilityDecision(
                reliable=False,
                kind="unreliable",
                reason="missing_or_invalid_observables",
            )
        observable_check = _finite_required_values(
            observables,
            ED_REFERENCE_OBSERVABLES,
        )
        if observable_check is not None:
            return ReliabilityDecision(False, "unreliable", observable_check)
        return ReliabilityDecision(
            reliable=True,
            kind="exact_free_boson",
            reason="stable_noninteracting_grand_canonical_reference",
        )
    try:
        completed_shells = int(result.get("completed_shells"))
    except (TypeError, ValueError):
        return ReliabilityDecision(
            reliable=False,
            kind="unreliable",
            reason="missing_or_invalid_completed_shells",
        )

    observables = result.get("observables", {})
    relative_last_shell = result.get("relative_last_shell", {})
    policy = result.get("convergence_policy", {})
    if not isinstance(observables, dict):
        return ReliabilityDecision(
            reliable=False,
            kind="unreliable",
            reason="missing_or_invalid_observables",
        )
    if not isinstance(relative_last_shell, dict):
        return ReliabilityDecision(
            reliable=False,
            kind="unreliable",
            reason="missing_or_invalid_relative_last_shell",
        )
    if not isinstance(policy, dict):
        return ReliabilityDecision(
            reliable=False,
            kind="unreliable",
            reason="missing_or_invalid_convergence_policy",
        )
    observable_check = _finite_required_values(observables, ED_REFERENCE_OBSERVABLES)
    if observable_check is not None:
        return ReliabilityDecision(False, "unreliable", observable_check)
    tail_check = _finite_required_values(relative_last_shell, ED_REFERENCE_OBSERVABLES)
    if tail_check is not None:
        return ReliabilityDecision(False, "unreliable", tail_check)
    try:
        u2 = float(parameters.get("U2"))
    except (TypeError, ValueError):
        return ReliabilityDecision(
            reliable=False,
            kind="unreliable",
            reason="missing_or_invalid_u2",
        )

    if status == "converged" and u2 >= 0.0:
        try:
            tolerance = float(policy["tail_tolerance"])
            min_completed_shells = int(policy.get("min_completed_shells", 0))
        except (KeyError, TypeError, ValueError):
            return ReliabilityDecision(
                reliable=False,
                kind="unreliable",
                reason="missing_or_invalid_tail_tolerance",
            )
        if completed_shells < min_completed_shells:
            return ReliabilityDecision(
                reliable=False,
                kind="unreliable",
                reason="insufficient_completed_shells",
            )
        max_tail = max(float(relative_last_shell[name]) for name in ED_REFERENCE_OBSERVABLES)
        if max_tail > tolerance:
            return ReliabilityDecision(
                reliable=False,
                kind="unreliable",
                reason="tail_not_converged",
            )
        return ReliabilityDecision(
            reliable=True,
            kind="tail_converged",
            reason="converged_nonnegative_u2",
        )
    if status == "low_density_cutoff_accepted" and u2 < 0.0:
        try:
            max_density_total = float(policy["max_density_total"])
            min_completed_shells = int(policy["min_completed_shells"])
            density_total = float(observables["density_total"])
        except (KeyError, TypeError, ValueError):
            return ReliabilityDecision(
                reliable=False,
                kind="unreliable",
                reason="missing_or_invalid_low_density_cutoff_policy",
            )
        if density_total > max_density_total:
            return ReliabilityDecision(
                reliable=False,
                kind="unreliable",
                reason="density_exceeds_cutoff",
            )
        if completed_shells < min_completed_shells:
            return ReliabilityDecision(
                reliable=False,
                kind="unreliable",
                reason="insufficient_completed_shells",
            )
        if completed_shells < NEGATIVE_U2_TRUSTED_LOW_DENSITY_MIN_SHELLS:
            return ReliabilityDecision(
                reliable=False,
                kind="unreliable",
                reason="insufficient_low_density_window_shells",
            )
        if density_total > NEGATIVE_U2_TRUSTED_LOW_DENSITY_MAX_TOTAL:
            return ReliabilityDecision(
                reliable=False,
                kind="unreliable",
                reason="density_exceeds_trusted_low_density_window",
            )
        return ReliabilityDecision(
            reliable=True,
            kind="finite_window",
            reason="negative_u2_low_density_cutoff_accepted",
        )
    return ReliabilityDecision(
        reliable=False,
        kind="unreliable",
        reason=f"unsupported_ed_status: status={status}, U2={u2}",
    )


def _finite_required_values(values: dict[str, Any], required: set[str]) -> str | None:
    missing = sorted(required - set(values))
    if missing:
        return f"missing_required_values: {','.join(missing)}"
    for name in sorted(required):
        try:
            value = float(values[name])
        except (TypeError, ValueError):
            return f"invalid_numeric_value: {name}"
        if not math.isfinite(value):
            return f"non_finite_value: {name}"
    return None


def _is_free_boson_case(parameters: dict[str, Any]) -> bool:
    try:
        return math.isclose(
            float(parameters.get("U1")),
            0.0,
            rel_tol=0.0,
            abs_tol=1.0e-14,
        ) and math.isclose(
            float(parameters.get("U2")),
            0.0,
            rel_tol=0.0,
            abs_tol=1.0e-14,
        )
    except (TypeError, ValueError):
        return False


def _triangular_single_particle_energy(
    lx: int,
    ly: int,
    kx: int,
    ky: int,
    t: float,
) -> float:
    phase_x = float(kx) / float(lx)
    phase_y = float(ky) / float(ly)
    return 2.0 * t * (
        math.cos(2.0 * math.pi * phase_x)
        + math.cos(2.0 * math.pi * phase_y)
        + math.cos(2.0 * math.pi * (phase_y - phase_x))
    )


def _finite_float(value: Any, name: str) -> float:
    try:
        number = float(value)
    except (TypeError, ValueError) as exc:
        raise ValueError(f"{name} is not numeric: {value!r}") from exc
    if not math.isfinite(number):
        raise ValueError(f"{name} is not finite: {value!r}")
    return number


def _float_or_nan(value: Any) -> float:
    try:
        number = float(value)
    except (TypeError, ValueError):
        return math.nan
    return number if math.isfinite(number) else math.nan


def _finite_or_none(value: Any) -> float | None:
    try:
        number = float(value)
    except (TypeError, ValueError):
        return None
    return number if math.isfinite(number) else None


def _numbers_match(actual: Any, expected: Any, abs_tol: float) -> bool:
    try:
        actual_number = float(actual)
        expected_number = float(expected)
    except (TypeError, ValueError):
        return actual == expected
    return math.isclose(
        actual_number,
        expected_number,
        rel_tol=0.0,
        abs_tol=abs_tol,
    )


def _unavailable_observable(
    observable: dict[str, Any] | None,
    *,
    reliable: bool,
    reason: str,
) -> dict[str, Any]:
    payload = dict(observable) if isinstance(observable, dict) else {}
    payload.setdefault("actual", math.nan)
    payload.setdefault("stderr", math.nan)
    payload["reliable"] = reliable
    payload["reason"] = reason
    return payload


def _comparison_base(
    dqmc_observable: dict[str, Any],
    ed_observable: dict[str, Any],
    *,
    stderr_tolerance: float,
    atol: float,
    rtol: float,
) -> dict[str, float | bool]:
    return {
        "dqmc": _float_or_nan(dqmc_observable.get("actual")),
        "ed": _float_or_nan(ed_observable.get("actual")),
        "difference": math.nan,
        "stderr": _float_or_nan(dqmc_observable.get("stderr")),
        "z_score": math.nan,
        "threshold": math.nan,
        "stderr_tolerance": float(stderr_tolerance),
        "atol": float(atol),
        "rtol": float(rtol),
        "pass_uncertainty": False,
        **_imaginary_diagnostics(
            dqmc_observable,
            stderr_tolerance=stderr_tolerance,
            atol=atol,
        ),
        "reliable": False,
    }


def _imaginary_diagnostics(
    dqmc_observable: dict[str, Any],
    *,
    stderr_tolerance: float,
    atol: float,
) -> dict[str, float | bool]:
    if "imag_actual" not in dqmc_observable:
        return {
            "imag_actual": math.nan,
            "imag_stderr": math.nan,
            "imag_z_score": math.nan,
            "imag_threshold": math.nan,
            "pass_imaginary": True,
        }

    imag_actual = _finite_or_none(dqmc_observable.get("imag_actual"))
    imag_stderr = _finite_or_none(dqmc_observable.get("imag_stderr"))
    if imag_actual is None or imag_stderr is None or imag_stderr < 0.0:
        return {
            "imag_actual": _float_or_nan(dqmc_observable.get("imag_actual")),
            "imag_stderr": _float_or_nan(dqmc_observable.get("imag_stderr")),
            "imag_z_score": math.nan,
            "imag_threshold": math.nan,
            "pass_imaginary": False,
        }
    threshold = max(
        float(atol),
        IMAGINARY_ROUNDOFF_ATOL,
        float(stderr_tolerance) * imag_stderr,
    )
    if imag_stderr > 0.0:
        z_score = imag_actual / imag_stderr
    elif imag_actual == 0.0:
        z_score = 0.0
    else:
        z_score = math.copysign(math.inf, imag_actual)
    return {
        "imag_actual": imag_actual,
        "imag_stderr": imag_stderr,
        "imag_z_score": z_score,
        "imag_threshold": threshold,
        "pass_imaginary": abs(imag_actual) <= threshold,
    }


def _first_unreliable_reason(
    dqmc_observable: dict[str, Any],
    ed_observable: dict[str, Any],
    *,
    dqmc_reliable: bool,
    ed_reliable: bool,
    fallback: str,
) -> str:
    if not dqmc_reliable:
        return str(dqmc_observable.get("reason", "dqmc_observable_unreliable"))
    if not ed_reliable:
        return str(ed_observable.get("reason", "ed_observable_unreliable"))
    return fallback


def _case_status(
    *,
    trusted: bool,
    ed_match: ReliabilityDecision,
    ed_decision: ReliabilityDecision,
) -> str:
    if trusted:
        return "trusted"
    if not ed_match.reliable:
        return "ed_parameter_mismatch"
    if not ed_decision.reliable:
        return "ed_unreliable"
    return "comparison_failed"


def _structure_factor_summary(
    samples: list[complex], block_size: int, skip_samples: int
) -> dict[str, float | int]:
    stats = complex_block_statistics(samples, block_size, skip_samples)
    real = stats["real"]
    imag = stats["imag"]
    return {
        **real,
        "imag_actual": imag["actual"],
        "imag_stderr": imag["stderr"],
    }


def _ensure_equal_lengths(series_by_name: dict[str, list[Any]]) -> None:
    lengths = {len(values) for values in series_by_name.values()}
    if len(lengths) == 1:
        return

    rendered = ", ".join(
        f"{name}={len(values)}" for name, values in series_by_name.items()
    )
    raise ValueError(f"DQMC output files have different sample counts: {rendered}")


def _ipr_summary(
    density_site_series: list[list[float]],
    block_size: int,
    skip_samples: int,
    min_ipr_number: float,
) -> dict[str, float | int | bool | str]:
    block_profiles, samples_used = _vector_block_means(
        density_site_series, block_size, skip_samples
    )
    block_iprs: list[float] = []
    for index, profile in enumerate(block_profiles, start=1):
        block_ipr = compute_ipr(
            profile,
            min_number=min_ipr_number,
        )
        if not block_ipr.reliable:
            return {
                "actual": math.nan,
                "stderr": math.nan,
                "reliable": False,
                "reason": f"block_{index}_{block_ipr.reason}",
                "blocks": len(block_profiles),
                "samples_used": samples_used,
                "skip_samples": skip_samples,
            }
        block_iprs.append(block_ipr.value)

    width = len(block_profiles[0])
    mean_profile = [
        mean(profile[column] for profile in block_profiles) for column in range(width)
    ]
    actual_ipr = compute_ipr(
        mean_profile,
        min_number=min_ipr_number,
    )
    if not actual_ipr.reliable:
        return {
            "actual": math.nan,
            "stderr": math.nan,
            "reliable": False,
            "reason": f"blocked_mean_{actual_ipr.reason}",
            "blocks": len(block_iprs),
            "samples_used": samples_used,
            "skip_samples": skip_samples,
        }

    block_ipr_mean = mean(block_iprs)
    variance = sum((value - block_ipr_mean) ** 2 for value in block_iprs) / float(
        len(block_iprs) - 1
    )
    return {
        "actual": actual_ipr.value,
        "stderr": math.sqrt(variance / float(len(block_iprs))),
        "reliable": True,
        "reason": actual_ipr.reason,
        "blocks": len(block_iprs),
        "samples_used": samples_used,
        "skip_samples": skip_samples,
    }


def summarize_dqmc_run(
    run_dir: Path,
    block_size: int,
    lq: int,
    skip_samples: int = 0,
    min_ipr_number: float = 1.0e-8,
) -> dict[str, Any]:
    if lq <= 0:
        raise ValueError("lq must be positive")

    scalar_series = {
        name: read_scalar_series(run_dir / name)
        for name in (
            "density_total",
            "energy_density",
            "doubleOcc",
            "num_up",
            "num_do",
            "onsite_n2_up",
            "onsite_n2_do",
        )
    }
    density_site_series = read_vector_series(run_dir / "density_site_total", lq)
    complex_series = {
        name: read_complex_series(run_dir / name)
        for name in ("sf_K", "psf_Gamma", "dw_K")
    }
    _ensure_equal_lengths(
        {**scalar_series, **complex_series, "density_site_total": density_site_series}
    )

    observables: dict[str, dict[str, float | int | bool | str]] = {
        "density_total": block_statistics(
            scalar_series["density_total"], block_size, skip_samples
        ),
        "energy_density": block_statistics(
            scalar_series["energy_density"], block_size, skip_samples
        ),
        "doubleOcc": block_statistics(
            scalar_series["doubleOcc"], block_size, skip_samples
        ),
        "S_SF_K": _structure_factor_summary(
            complex_series["sf_K"], block_size, skip_samples
        ),
        "S_PSF_Gamma": _structure_factor_summary(
            complex_series["psf_Gamma"], block_size, skip_samples
        ),
        "S_DW_K": _structure_factor_summary(
            complex_series["dw_K"], block_size, skip_samples
        ),
    }

    observables["IPR"] = _ipr_summary(
        density_site_series, block_size, skip_samples, min_ipr_number
    )

    return {
        "run_dir": str(run_dir),
        "lq": lq,
        "block_size": block_size,
        "skip_samples": skip_samples,
        "observables": observables,
    }
