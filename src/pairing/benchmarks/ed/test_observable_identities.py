#!/usr/bin/env python3
from __future__ import annotations

import argparse

import numpy as np

import ed_pairing_triangle_general as general

FD_STEP = 1.0e-4
TOL = 5.0e-5

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


def observable_value(observables: dict[str, object], key: str) -> float:
    if key not in observables:
        raise AssertionError(
            f"missing ED observable {key!r}; available keys: {sorted(observables)}"
        )
    return float(observables[key]["value"])


def free_energy_at(params: dict[str, object], key: str, value: float) -> float:
    shifted = dict(params)
    shifted[key] = value
    result = general.run_ed(shifted, basis_cap=int(shifted["basis_cap"]), allow_large_basis=False)
    return float(result["free_energy"])


def assert_close(name: str, actual: float, expected: float) -> None:
    if abs(actual - expected) > TOL:
        raise AssertionError(f"{name}: actual={actual:.16g}, expected={expected:.16g}, tol={TOL:.1e}")


def run_tiny() -> None:
    result = general.run_ed(TINY_PARAMS, basis_cap=int(TINY_PARAMS["basis_cap"]), allow_large_basis=False)
    observables = result["observables"]
    nsite = int(TINY_PARAMS["Lx"]) * int(TINY_PARAMS["Ly"])
    for unavailable_k in ("S_SF_K", "S_DW_K"):
        if unavailable_k in observables:
            raise AssertionError(f"{unavailable_k} should be unavailable for non-3-compatible ED")

    num_up = observable_value(observables, "num_up")
    num_do = observable_value(observables, "num_do")
    total_ne = observable_value(observables, "total_NE")
    total_density = observable_value(observables, "total_density")
    density_total = observable_value(observables, "density_total")
    pair_equal = observable_value(observables, "pair_equal")
    kinetic_total = observable_value(observables, "kinetic_total")
    interaction_energy_total = observable_value(observables, "interaction_energy_total")
    pairing_total = observable_value(observables, "pairing_total")
    energy_total = observable_value(observables, "energy_total")
    energy_density = observable_value(observables, "energy_density")
    chemical_energy_total = observable_value(observables, "chemical_energy_total")
    chemical_energy_density = observable_value(observables, "chemical_energy_density")
    grand_energy_total = observable_value(observables, "grand_energy_total")
    grand_energy_density = observable_value(observables, "grand_energy_density")
    onsite_n2_up = observable_value(observables, "onsite_n2_up")
    onsite_n2_do = observable_value(observables, "onsite_n2_do")
    same_flavor_pair = observable_value(observables, "same_flavor_pair")

    mu = float(TINY_PARAMS["mu"])

    assert_close("total_NE", total_ne, num_up + num_do)
    assert_close("total_density", total_density, total_ne / nsite)
    assert_close("density_total", density_total, total_ne / nsite)
    assert_close(
        "energy_total",
        energy_total,
        kinetic_total + interaction_energy_total + pairing_total,
    )
    assert_close(
        "energy_density",
        energy_density,
        (kinetic_total + interaction_energy_total + pairing_total) / nsite,
    )
    assert_close("energy_density from total", energy_density, energy_total / nsite)
    assert_close("chemical_energy_total", chemical_energy_total, -mu * total_ne)
    assert_close(
        "chemical_energy_density from total",
        chemical_energy_density,
        chemical_energy_total / nsite,
    )
    assert_close("chemical_energy_density", chemical_energy_density, -mu * density_total)
    assert_close(
        "grand_energy_total",
        grand_energy_total,
        energy_total + chemical_energy_total,
    )
    assert_close(
        "grand_energy_density from total",
        grand_energy_density,
        grand_energy_total / nsite,
    )
    assert_close(
        "grand_energy_density",
        grand_energy_density,
        energy_density + chemical_energy_density,
    )
    assert_close(
        "same_flavor_pair",
        same_flavor_pair,
        0.5 * ((onsite_n2_up - num_up) + (onsite_n2_do - num_do)) / nsite,
    )

    delta = float(TINY_PARAMS["Delta"])
    d_delta = (
        free_energy_at(TINY_PARAMS, "Delta", delta + FD_STEP)
        - free_energy_at(TINY_PARAMS, "Delta", delta - FD_STEP)
    ) / (2.0 * FD_STEP)
    assert_close("Nsite * pair_equal", nsite * pair_equal, d_delta)

    t = float(TINY_PARAMS["t"])
    d_t = (
        free_energy_at(TINY_PARAMS, "t", t + FD_STEP)
        - free_energy_at(TINY_PARAMS, "t", t - FD_STEP)
    ) / (2.0 * FD_STEP)
    assert_close("kinetic_total / t", kinetic_total / t, d_t)

    basis = general.build_basis(1, 1, 2, None)
    operator = general.number_operator(1, basis, "b")
    rng = np.random.default_rng(12345)
    trial = rng.normal(size=(basis.Ns, 2))
    evecs, _ = np.linalg.qr(trial)
    dense_diag = np.einsum(
        "ia,ij,ja->a",
        evecs.conj(),
        operator.toarray(),
        evecs,
        optimize=True,
    )
    sparse_diag = general.operator_diag_in_eigenbasis_sparse_safe(operator, evecs)
    if not np.allclose(sparse_diag, dense_diag, rtol=1.0e-12, atol=1.0e-12):
        raise AssertionError(
            "sparse-safe operator diagonal disagrees with dense toarray reference"
        )


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--case", choices=["tiny"], default="tiny")
    args = parser.parse_args()
    if args.case == "tiny":
        run_tiny()
    print("observable identity tests passed")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
