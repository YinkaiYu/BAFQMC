#!/usr/bin/env python3
from __future__ import annotations

import argparse
import json
from pathlib import Path
from typing import Any

import numpy as np
from quspin.basis import boson_basis_general, tensor_basis

if __package__:
    from . import ed_pairing_triangle_general as general
    from .geometry_triangle import triangular_bonds
else:
    import ed_pairing_triangle_general as general
    from geometry_triangle import triangular_bonds


def build_tensor_basis(Lx: int, Ly: int, nmax: int):
    nsite = Lx * Ly
    basis_b = boson_basis_general(nsite, sps=nmax + 1)
    basis_c = boson_basis_general(nsite, sps=nmax + 1)
    return tensor_basis(basis_b, basis_c)


def tensor_basis_size_estimate(Lx: int, Ly: int, nmax: int) -> int:
    """Return the exact tensor-basis size before constructing QuSpin bases."""

    nsite = int(Lx) * int(Ly)
    sps = int(nmax) + 1
    if nsite <= 0:
        raise ValueError("Lx*Ly must be positive")
    if sps <= 0:
        raise ValueError("nmax must be non-negative")
    return int(sps ** (2 * nsite))


def assert_tensor_prebasis_budget(
    params: dict[str, Any],
    *,
    basis_cap: int,
    allow_large_basis: bool,
    dense_memory_cap_gib: float | None,
) -> int:
    """Fail closed before QuSpin constructs large tensor bases."""

    estimated_basis_size = tensor_basis_size_estimate(
        int(params["Lx"]),
        int(params["Ly"]),
        int(params["nmax"]),
    )
    if estimated_basis_size > int(basis_cap) and not allow_large_basis:
        raise SystemExit(
            f"estimated tensor basis.Ns={estimated_basis_size} exceeds basis_cap={int(basis_cap)}; "
            "use --allow-large-basis for manual references"
        )
    general.assert_dense_memory_budget(
        estimated_basis_size,
        dense_memory_cap_gib=dense_memory_cap_gib,
    )
    return estimated_basis_size


def build_tensor_hamiltonian(params: dict[str, Any], basis):
    nsite = int(params["Lx"]) * int(params["Ly"])
    bonds = triangular_bonds(int(params["Lx"]), int(params["Ly"]))
    t = float(params["t"])
    u1 = float(params["U1"])
    u2 = float(params["U2"])
    mu = float(params["mu"])
    delta = float(params["Delta"])
    hop = [[t, i, j] for i, j in bonds] + [[t, j, i] for i, j in bonds]
    static = [
        ["+-|", hop],
        ["|+-", hop],
        ["n|", [[-mu, i] for i in range(nsite)]],
        ["|n", [[-mu, i] for i in range(nsite)]],
        ["nn|", [[u1 + u2, i, i] for i in range(nsite)]],
        ["|nn", [[u1 + u2, i, i] for i in range(nsite)]],
        # Paper convention: U1 multiplies relative density and U2 total density.
        ["n|n", [[2.0 * (u2 - u1), i, i] for i in range(nsite)]],
        ["+|+", [[delta, i, i] for i in range(nsite)]],
        ["-|-", [[delta, i, i] for i in range(nsite)]],
    ]
    return general.make_operator(static, basis)


def make_tensor_operator(
    static: list[list[Any]],
    basis,
    dtype=np.float64,
    check_herm: bool = True,
):
    return general.make_operator(static, basis, dtype=dtype, check_herm=check_herm)


def tensor_number_operator(nsite: int, basis, layer: str):
    op = "n|" if layer == "b" else "|n"
    return make_tensor_operator([[op, [[1.0, i] for i in range(nsite)]]], basis)


def tensor_kinetic_operator(params: dict[str, Any], basis):
    nsite = int(params["Lx"]) * int(params["Ly"])
    bonds = triangular_bonds(int(params["Lx"]), int(params["Ly"]))
    t = float(params["t"])
    hop = [[t, i, j] for i, j in bonds] + [[t, j, i] for i, j in bonds]
    return make_tensor_operator([["+-|", hop], ["|+-", hop]], basis)


def tensor_double_occ_operator(nsite: int, basis):
    return make_tensor_operator([["n|n", [[1.0 / nsite, i, i] for i in range(nsite)]]], basis)


def tensor_local_numsquare_operator(nsite: int, basis):
    up_terms = [[1.0 / nsite, i, i] for i in range(nsite)]
    do_terms = [[1.0 / nsite, i, i] for i in range(nsite)]
    return make_tensor_operator([["nn|", up_terms], ["|nn", do_terms]], basis)


def tensor_numsquare_operator(nsite: int, basis, layer: str):
    op = "nn|" if layer == "b" else "|nn"
    terms = [[1.0, i, j] for i in range(nsite) for j in range(nsite)]
    return make_tensor_operator([[op, terms]], basis)


def tensor_onsite_n2_operator(nsite: int, basis, layer: str):
    op = "nn|" if layer == "b" else "|nn"
    terms = [[1.0, i, i] for i in range(nsite)]
    return make_tensor_operator([[op, terms]], basis)


def tensor_pair_equal_operator(nsite: int, basis):
    terms = [[1.0 / nsite, i, i] for i in range(nsite)]
    return make_tensor_operator([["+|+", terms], ["-|-", terms]], basis)


def tensor_sf_structure_operator(params: dict[str, Any], basis):
    lx = int(params["Lx"])
    ly = int(params["Ly"])
    nsite = lx * ly
    norm = 1.0 / float(nsite * nsite)
    phases = general.triangular_phase_matrix(lx, ly, general.structure_k_index(lx, ly))
    terms = [[norm * phases[i, j], i, j] for i in range(nsite) for j in range(nsite)]
    return make_tensor_operator(
        [["+-|", terms], ["|+-", terms]],
        basis,
        dtype=np.complex128,
        check_herm=False,
    )


def tensor_dw_structure_operator(params: dict[str, Any], basis):
    lx = int(params["Lx"])
    ly = int(params["Ly"])
    nsite = lx * ly
    norm = 1.0 / float(nsite * nsite)
    phases = general.triangular_phase_matrix(lx, ly, general.structure_k_index(lx, ly))
    same_terms = []
    cross_terms = []
    for i in range(nsite):
        for j in range(nsite):
            coef = norm * phases[i, j]
            same_terms.append([coef, i, j])
            cross_terms.append([coef, i, j])
            cross_terms.append([coef, j, i])
    return make_tensor_operator(
        [["nn|", same_terms], ["|nn", same_terms], ["n|n", cross_terms]],
        basis,
        dtype=np.complex128,
        check_herm=False,
    )


def tensor_psf_gamma_operator(nsite: int, basis):
    norm = 1.0 / float(nsite * nsite)
    terms = [[norm, i, j, i, j] for i in range(nsite) for j in range(nsite)]
    return make_tensor_operator(
        [["+-|+-", terms]],
        basis,
        dtype=np.complex128,
        check_herm=False,
    )


def tensor_observable_values(
    params: dict[str, Any],
    basis,
    evals: np.ndarray,
    evecs: np.ndarray,
) -> dict[str, float]:
    beta = float(params["beta"])
    nsite = int(params["Lx"]) * int(params["Ly"])

    def avg(operator, name: str) -> float:
        op_diag = general.operator_diag_in_eigenbasis(operator, evecs)
        value, _, _ = general.thermal_average(evals, op_diag, beta)
        return general.real_physical_average(value, name)

    num_up = avg(tensor_number_operator(nsite, basis, "b"), "num_up")
    num_do = avg(tensor_number_operator(nsite, basis, "c"), "num_do")
    total_ne = num_up + num_do
    density_total = total_ne / nsite
    kinetic_total = avg(tensor_kinetic_operator(params, basis), "kinetic_total")
    double_occ = avg(tensor_double_occ_operator(nsite, basis), "doubleOcc")
    onsite_n2_up = avg(tensor_onsite_n2_operator(nsite, basis, "b"), "onsite_n2_up")
    onsite_n2_do = avg(tensor_onsite_n2_operator(nsite, basis, "c"), "onsite_n2_do")
    pair_equal = avg(tensor_pair_equal_operator(nsite, basis), "pair_equal")
    u1 = float(params["U1"])
    u2 = float(params["U2"])
    mu = float(params["mu"])
    delta = float(params["Delta"])
    interaction_energy_total = (
        (u1 + u2) * (onsite_n2_up + onsite_n2_do)
        + 2.0 * (u2 - u1) * float(nsite) * double_occ
    )
    pairing_total = delta * float(nsite) * pair_equal
    chemical_energy_total = -mu * total_ne
    energy_total = kinetic_total + interaction_energy_total + pairing_total
    grand_energy_total = energy_total + chemical_energy_total
    values = {
        "num_up": num_up,
        "num_do": num_do,
        "total_NE": total_ne,
        "total_density": density_total,
        "density_total": density_total,
        "kinetic_total": kinetic_total,
        "doubleOcc": double_occ,
        "local_numsquare": (onsite_n2_up + onsite_n2_do) / nsite,
        "numsquare_up": avg(tensor_numsquare_operator(nsite, basis, "b"), "numsquare_up"),
        "numsquare_do": avg(tensor_numsquare_operator(nsite, basis, "c"), "numsquare_do"),
        "pair_equal": pair_equal,
        "onsite_n2_up": onsite_n2_up,
        "onsite_n2_do": onsite_n2_do,
        "interaction_energy_total": interaction_energy_total,
        "interaction_energy_density": interaction_energy_total / nsite,
        "pairing_total": pairing_total,
        "pairing_energy_density": pairing_total / nsite,
        "chemical_energy_total": chemical_energy_total,
        "chemical_energy_density": chemical_energy_total / nsite,
        "energy_total": energy_total,
        "energy_density": energy_total / nsite,
        "grand_energy_total": grand_energy_total,
        "grand_energy_density": grand_energy_total / nsite,
        "S_PSF_Gamma": avg(tensor_psf_gamma_operator(nsite, basis), "S_PSF_Gamma"),
        "same_flavor_pair": 0.5 * ((onsite_n2_up - num_up) + (onsite_n2_do - num_do)) / nsite,
    }
    if general.is_k_geometry_compatible(int(params["Lx"]), int(params["Ly"])):
        values.update(
            {
                "S_SF_K": avg(tensor_sf_structure_operator(params, basis), "S_SF_K"),
                "S_DW_K": avg(tensor_dw_structure_operator(params, basis), "S_DW_K"),
            }
        )
    return values


def tensor_cutoff_diagnostics(basis, evecs: np.ndarray, evals: np.ndarray, params: dict[str, Any]) -> dict[str, float | int]:
    """Return ordering-invariant cutoff diagnostics for tensor-product states.

    If a future diagnostic becomes site-resolved or flavor-resolved, do not reuse
    integer decoding blindly; first add a self-test against basis.int_to_state()
    or use QuSpin operator expectation values.
    """
    nsite = int(params["Lx"]) * int(params["Ly"])
    nmax = int(params["nmax"])
    sps = nmax + 1
    left = np.array(
        [general.decode_state(int(state), nsite, sps) for state in basis.basis_left.states],
        dtype=np.int16,
    )
    right = np.array(
        [general.decode_state(int(state), nsite, sps) for state in basis.basis_right.states],
        dtype=np.int16,
    )

    total_particles = []
    local_boundary = []
    for occ_b in left:
        for occ_c in right:
            occ = np.concatenate((occ_b, occ_c))
            total_particles.append(int(occ.sum()))
            local_boundary.append(float(np.any(occ == nmax)))

    total_particles_arr = np.array(total_particles, dtype=np.int16)
    local_boundary_arr = np.array(local_boundary, dtype=np.float64)
    local_boundary_diag = general.fock_diag_in_eigenbasis(local_boundary_arr, evecs)
    local_weight, _, _ = general.thermal_average(evals, local_boundary_diag, float(params["beta"]))
    return {
        "max_total_particles_in_basis": int(total_particles_arr.max()),
        "total_cutoff_boundary_weight": 0.0,
        "local_nmax_boundary_weight": float(local_weight),
    }


def run_tensor_ed(
    params: dict[str, Any],
    basis_cap: int | None = None,
    allow_large_basis: bool = False,
    dense_memory_cap_gib: float | None = general.DEFAULT_DENSE_MEMORY_CAP_GIB,
) -> dict[str, Any]:
    normalized = general.normalize_params(params)
    if normalized["ncut"] is not None:
        raise SystemExit("tensor_basis does not impose a combined N_b+N_c cutoff; use ncut=null")
    if basis_cap is not None:
        normalized["basis_cap"] = int(basis_cap)

    estimated_basis_size = assert_tensor_prebasis_budget(
        normalized,
        basis_cap=int(normalized["basis_cap"]),
        allow_large_basis=allow_large_basis,
        dense_memory_cap_gib=dense_memory_cap_gib,
    )
    basis = build_tensor_basis(int(normalized["Lx"]), int(normalized["Ly"]), int(normalized["nmax"]))
    if basis.Ns != estimated_basis_size:
        raise AssertionError(
            f"tensor basis estimate {estimated_basis_size} disagrees with QuSpin basis.Ns={basis.Ns}"
        )
    if basis.Ns > int(normalized["basis_cap"]) and not allow_large_basis:
        raise SystemExit(
            f"basis.Ns={basis.Ns} exceeds basis_cap={normalized['basis_cap']}; "
            "use --allow-large-basis for manual references"
        )
    general.assert_dense_memory_budget(basis.Ns, dense_memory_cap_gib=dense_memory_cap_gib)

    hamiltonian_mu = build_tensor_hamiltonian(normalized, basis)
    evals, evecs = hamiltonian_mu.eigh()
    observables = tensor_observable_values(normalized, basis, evals, evecs)
    diagnostics = tensor_cutoff_diagnostics(basis, evecs, evals, normalized)
    diagnostics.update(
        general.dense_memory_diagnostics(basis.Ns, dense_memory_cap_gib=dense_memory_cap_gib)
    )
    return general.result_payload(normalized, basis.Ns, evals, observables, diagnostics, "dense_H_eigh_tensor_basis")


def load_params(path: str | Path) -> dict[str, Any]:
    return general.load_params(path)


def write_result(result: dict[str, Any], output: str | None) -> None:
    if output is None:
        print(json.dumps(result, indent=2))
        return
    path = Path(output)
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8") as handle:
        json.dump(result, handle, indent=2)
        handle.write("\n")


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--params", required=True, help="Input ED parameter JSON")
    parser.add_argument("--output", help="Output JSON path; stdout is used when omitted")
    parser.add_argument("--basis-cap", type=int, help="Override params basis_cap")
    parser.add_argument("--allow-large-basis", action="store_true")
    parser.add_argument(
        "--dense-memory-cap-gib",
        type=float,
        default=general.DEFAULT_DENSE_MEMORY_CAP_GIB,
        help=(
            "Refuse dense ED when estimated peak RSS exceeds this cap; "
            "set 0 to disable after reserving sufficient memory"
        ),
    )
    args = parser.parse_args()

    params = load_params(args.params)
    result = run_tensor_ed(
        params,
        basis_cap=args.basis_cap,
        allow_large_basis=args.allow_large_basis,
        dense_memory_cap_gib=args.dense_memory_cap_gib,
    )
    write_result(result, args.output)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
