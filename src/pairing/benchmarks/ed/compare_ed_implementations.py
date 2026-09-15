#!/usr/bin/env python3
from __future__ import annotations

import argparse
import math

if __package__:
    from . import ed_pairing_triangle_general as general
    from . import ed_pairing_triangle_tensor as tensor
else:
    import ed_pairing_triangle_general as general
    import ed_pairing_triangle_tensor as tensor

ABS_TOL = 1.0e-9
REL_TOL = 1.0e-9

REQUIRED_OBSERVABLE_KEYS = {
    "density_total",
    "kinetic_total",
    "interaction_energy_total",
    "interaction_energy_density",
    "pairing_total",
    "pairing_energy_density",
    "chemical_energy_total",
    "chemical_energy_density",
    "energy_total",
    "energy_density",
    "grand_energy_total",
    "grand_energy_density",
    "onsite_n2_up",
    "onsite_n2_do",
    "S_PSF_Gamma",
    "same_flavor_pair",
}

TINY_PARAMS = {
    "case": "triangle_pairing_2x1_tensor_equivalence",
    "Lx": 2,
    "Ly": 1,
    "t": 1.0,
    "U1": -0.1,
    "U2": 1.0,
    "mu": -5.0,
    "Delta": 0.2,
    "beta": 1.0,
    "nmax": 2,
    "ncut": None,
    "basis_cap": 500,
}


def assert_observables_close(combined: dict[str, object], tensor_result: dict[str, object]) -> None:
    combined_observables = combined["observables"]
    tensor_observables = tensor_result["observables"]
    combined_keys = set(combined_observables)
    tensor_keys = set(tensor_observables)
    missing_required = REQUIRED_OBSERVABLE_KEYS - (combined_keys & tensor_keys)
    if missing_required:
        raise AssertionError(
            f"missing required ED observable keys: {sorted(missing_required)}"
        )
    if combined_keys != tensor_keys:
        raise AssertionError(
            f"observable key mismatch: combined_only={sorted(combined_keys - tensor_keys)}, "
            f"tensor_only={sorted(tensor_keys - combined_keys)}"
        )

    failures: list[str] = []
    for key in sorted(combined_keys):
        left = combined_observables[key]
        right = tensor_observables[key]
        if left["normalization"] != right["normalization"]:
            failures.append(f"{key}: normalization {left['normalization']} != {right['normalization']}")
            continue
        left_value = float(left["value"])
        right_value = float(right["value"])
        if not math.isclose(left_value, right_value, rel_tol=REL_TOL, abs_tol=ABS_TOL):
            failures.append(f"{key}: combined={left_value:.16g}, tensor={right_value:.16g}")

    if not math.isclose(
        float(combined["free_energy"]),
        float(tensor_result["free_energy"]),
        rel_tol=REL_TOL,
        abs_tol=ABS_TOL,
    ):
        failures.append(
            f"free_energy: combined={float(combined['free_energy']):.16g}, "
            f"tensor={float(tensor_result['free_energy']):.16g}"
        )

    if failures:
        joined = "\n".join(failures)
        raise AssertionError(f"combined and tensor ED results differ:\n{joined}")


def run_tiny() -> None:
    combined = general.run_ed(TINY_PARAMS, basis_cap=int(TINY_PARAMS["basis_cap"]), allow_large_basis=False)
    tensor_result = tensor.run_tensor_ed(TINY_PARAMS, basis_cap=int(TINY_PARAMS["basis_cap"]), allow_large_basis=False)
    assert_observables_close(combined, tensor_result)
    print(f"combined and tensor ED implementations agree for {len(combined['observables'])} observables")


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--case", choices=["tiny"], default="tiny")
    args = parser.parse_args()
    if args.case == "tiny":
        run_tiny()
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
