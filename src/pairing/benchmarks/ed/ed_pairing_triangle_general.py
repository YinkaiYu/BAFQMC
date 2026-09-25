#!/usr/bin/env python3
from __future__ import annotations

import argparse
import contextlib
import io
import json
import sys
from pathlib import Path
from typing import Any

import numpy as np
from quspin.basis import boson_basis_general
from quspin.operators import hamiltonian

if __package__:
    from .geometry_triangle import A1, A2, B1, B2, triangular_bonds
else:
    from geometry_triangle import A1, A2, B1, B2, triangular_bonds

OBSERVABLE_KEYS = (
    "num_up",
    "num_do",
    "total_NE",
    "total_density",
    "density_total",
    "kinetic_total",
    "doubleOcc",
    "local_numsquare",
    "numsquare_up",
    "numsquare_do",
    "pair_equal",
    "onsite_n2_up",
    "onsite_n2_do",
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
    "S_SF_K",
    "S_DW_K",
    "S_PSF_Gamma",
    "same_flavor_pair",
)

OBSERVABLE_NORMALIZATION = {
    "num_up": "total",
    "num_do": "total",
    "total_NE": "total",
    "total_density": "per_site",
    "density_total": "per_site",
    "kinetic_total": "total",
    "doubleOcc": "per_site",
    "local_numsquare": "per_site",
    "numsquare_up": "total",
    "numsquare_do": "total",
    "pair_equal": "per_site",
    "onsite_n2_up": "total",
    "onsite_n2_do": "total",
    "interaction_energy_total": "total",
    "interaction_energy_density": "per_site",
    "pairing_total": "total",
    "pairing_energy_density": "per_site",
    "chemical_energy_total": "total",
    "chemical_energy_density": "per_site",
    "energy_total": "total",
    "energy_density": "per_site",
    "grand_energy_total": "total",
    "grand_energy_density": "per_site",
    "S_SF_K": "per_site_squared",
    "S_DW_K": "per_site_squared",
    "S_PSF_Gamma": "per_site_squared",
    "same_flavor_pair": "per_site",
}

PARAMETER_KEYS = ("Lx", "Ly", "t", "U1", "U2", "mu", "Delta", "beta")
DEFAULT_BASIS_CAP = 2000
DENSE_EIGH_RSS_MATRIX_COPIES = 4.0
DEFAULT_DENSE_MEMORY_CAP_GIB = 12.0
ED_IMAGINARY_ATOL = 1.0e-10
ED_POSITIVE_OBSERVABLE_ATOL = 1.0e-10
K_OBSERVABLES = {"S_SF_K", "S_DW_K"}
POSITIVE_OBSERVABLES = K_OBSERVABLES | {"S_PSF_Gamma"}


def dense_eigh_memory_estimate_gib(basis_size: int) -> float:
    """Estimate peak dense ``eigh`` memory from the complex matrix footprint.

    The factor is calibrated to the 3x3 ``nmax=3,ncut=4`` full-trace run:
    a 7297-dimensional basis used 3.18 GiB peak RSS while one complex128 dense
    matrix is 0.79 GiB.  It is a guardrail, not a scheduler reservation.
    """

    basis_size = int(basis_size)
    if basis_size <= 0:
        raise ValueError("basis_size must be positive")
    one_dense_matrix_gib = 16.0 * basis_size * basis_size / (1024.0**3)
    return DENSE_EIGH_RSS_MATRIX_COPIES * one_dense_matrix_gib


def basis_size_estimate(num_modes: int, nmax: int, ncut: int | None) -> int:
    """Return the exact constrained boson-basis size before QuSpin allocation."""

    num_modes = int(num_modes)
    nmax = int(nmax)
    if num_modes <= 0:
        raise ValueError("num_modes must be positive")
    if nmax < 0:
        raise ValueError("nmax must be non-negative")
    if ncut is None:
        return (nmax + 1) ** num_modes

    ncut = int(ncut)
    if ncut < 0:
        raise ValueError("ncut must be non-negative or null")
    max_total = min(ncut, nmax * num_modes)
    counts = [0] * (max_total + 1)
    counts[0] = 1
    for _ in range(num_modes):
        next_counts = [0] * (max_total + 1)
        for total, value in enumerate(counts):
            if value == 0:
                continue
            for local_occ in range(min(nmax, max_total - total) + 1):
                next_counts[total + local_occ] += value
        counts = next_counts
    return sum(counts)


def assert_dense_memory_budget(
    basis_size: int,
    *,
    dense_memory_cap_gib: float | None = DEFAULT_DENSE_MEMORY_CAP_GIB,
) -> None:
    """Fail before constructing dense ED matrices when the estimate is unsafe."""

    if dense_memory_cap_gib is None or float(dense_memory_cap_gib) <= 0.0:
        return
    estimate = dense_eigh_memory_estimate_gib(int(basis_size))
    if estimate > float(dense_memory_cap_gib):
        raise SystemExit(
            f"estimated dense ED RSS {estimate:.1f} GiB for basis.Ns={int(basis_size)} "
            f"exceeds cap {float(dense_memory_cap_gib):.1f} GiB; "
            "increase --dense-memory-cap-gib only on a machine with enough RAM, "
            "or use the sparse checkpoint/HPC workflow"
        )


def dense_memory_diagnostics(
    basis_size: int,
    *,
    dense_memory_cap_gib: float | None,
) -> dict[str, float | str]:
    cap = None if dense_memory_cap_gib is None else float(dense_memory_cap_gib)
    return {
        "dense_eigh_memory_estimate_gib": float(dense_eigh_memory_estimate_gib(basis_size)),
        "dense_eigh_memory_cap_gib": float(cap) if cap is not None else 0.0,
        "dense_eigh_memory_model": (
            f"{DENSE_EIGH_RSS_MATRIX_COPIES:g} complex128 dense matrices"
        ),
    }


def normalize_params(params: dict[str, Any]) -> dict[str, Any]:
    required = ("case", *PARAMETER_KEYS, "nmax", "ncut")
    missing = [key for key in required if key not in params]
    if missing:
        raise ValueError(f"missing required parameter keys: {', '.join(missing)}")

    normalized = dict(params)
    normalized["case"] = str(params["case"])
    normalized["Lx"] = int(params["Lx"])
    normalized["Ly"] = int(params["Ly"])
    normalized["t"] = float(params["t"])
    normalized["U1"] = float(params["U1"])
    normalized["U2"] = float(params["U2"])
    normalized["mu"] = float(params["mu"])
    normalized["Delta"] = float(params["Delta"])
    normalized["beta"] = float(params["beta"])
    normalized["nmax"] = int(params["nmax"])
    normalized["ncut"] = None if params["ncut"] is None else int(params["ncut"])
    normalized["basis_cap"] = int(params.get("basis_cap", DEFAULT_BASIS_CAP))

    if normalized["Lx"] <= 0 or normalized["Ly"] <= 0:
        raise ValueError("Lx and Ly must be positive")
    if normalized["nmax"] < 0:
        raise ValueError("nmax must be non-negative")
    if normalized["ncut"] is not None and normalized["ncut"] < 0:
        raise ValueError("ncut must be non-negative or null")
    if normalized["beta"] <= 0.0:
        raise ValueError("beta must be positive")
    if normalized["basis_cap"] <= 0:
        raise ValueError("basis_cap must be positive")
    return normalized


def load_params(path: str | Path) -> dict[str, Any]:
    with Path(path).open("r", encoding="utf-8") as handle:
        return normalize_params(json.load(handle))


def build_basis(Lx: int, Ly: int, nmax: int, ncut: int | None):
    nsite = Lx * Ly
    nb = list(range(ncut + 1)) if ncut is not None else None
    return boson_basis_general(2 * nsite, sps=nmax + 1, Nb=nb)


def _checked_hamiltonian(
    static: list[list[Any]],
    basis,
    dtype=np.float64,
    check_herm: bool = True,
):
    with contextlib.redirect_stdout(io.StringIO()):
        return hamiltonian(
            static,
            [],
            basis=basis,
            dtype=dtype,
            check_symm=False,
            check_herm=check_herm,
            check_pcon=False,
        )


def build_hamiltonian(params: dict[str, Any], basis):
    nsite = int(params["Lx"]) * int(params["Ly"])
    t = float(params["t"])
    u1 = float(params["U1"])
    u2 = float(params["U2"])
    mu = float(params["mu"])
    delta = float(params["Delta"])
    bonds = triangular_bonds(int(params["Lx"]), int(params["Ly"]))

    hop_b = [[t, i, j] for i, j in bonds] + [[t, j, i] for i, j in bonds]
    hop_c = [[coef, i + nsite, j + nsite] for coef, i, j in hop_b]
    chem = [[-mu, i] for i in range(2 * nsite)]
    self_int = [[u1 + u2, i, i] for i in range(2 * nsite)]
    # Paper convention: U1 multiplies relative density and U2 total density.
    cross_int = [[2.0 * (u2 - u1), i, i + nsite] for i in range(nsite)]
    pair_create = [[delta, i, i + nsite] for i in range(nsite)]
    pair_annihilate = [[delta, i, i + nsite] for i in range(nsite)]
    static = [
        ["+-", hop_b],
        ["+-", hop_c],
        ["n", chem],
        ["nn", self_int],
        ["nn", cross_int],
        ["++", pair_create],
        ["--", pair_annihilate],
    ]
    return _checked_hamiltonian(static, basis)


def operator_diag_in_eigenbasis_sparse_safe(operator, evecs: np.ndarray) -> np.ndarray:
    """Return diagonal entries of an operator in an eigenvector basis.

    Avoid `operator.toarray()` so sparse low-energy checkpoints do not allocate
    dense observable matrices at large basis size.
    """

    applied = operator.dot(evecs)
    return np.einsum("ia,ia->a", evecs.conj(), applied, optimize=True)


def operator_diag_in_eigenbasis(operator, evecs: np.ndarray) -> np.ndarray:
    return operator_diag_in_eigenbasis_sparse_safe(operator, evecs)


def thermal_average(evals: np.ndarray, op_diag: np.ndarray, beta: float) -> tuple[complex, float, float]:
    energy_shift = float(evals.min())
    weights = np.exp(-beta * (evals - energy_shift))
    z_shifted = float(weights.sum())
    value = np.sum(weights * op_diag) / z_shifted
    return value, z_shifted, energy_shift


def free_energy(evals: np.ndarray, beta: float) -> float:
    energy_shift = float(evals.min())
    weights = np.exp(-beta * (evals - energy_shift))
    return energy_shift - float(np.log(weights.sum())) / beta


def real_physical_average(value: complex, name: str) -> float:
    real = float(np.real(value))
    imag = float(np.imag(value))
    tolerance = ED_IMAGINARY_ATOL * max(1.0, abs(real))
    if abs(imag) > tolerance:
        raise ValueError(f"ED observable {name} has nonzero imaginary part: {imag:.16g}")
    if name in POSITIVE_OBSERVABLES and real < -ED_POSITIVE_OBSERVABLE_ATOL:
        raise ValueError(f"ED observable {name} is negative: {real:.16g}")
    return 0.0 if name in POSITIVE_OBSERVABLES and real < 0.0 else real


def make_operator(
    static: list[list[Any]],
    basis,
    dtype=np.float64,
    check_herm: bool = True,
):
    return _checked_hamiltonian(static, basis, dtype=dtype, check_herm=check_herm)


def number_operator(nsite: int, basis, layer: str):
    offset = 0 if layer == "b" else nsite
    return make_operator([["n", [[1.0, offset + i] for i in range(nsite)]]], basis)


def kinetic_operator(params: dict[str, Any], basis):
    nsite = int(params["Lx"]) * int(params["Ly"])
    t = float(params["t"])
    bonds = triangular_bonds(int(params["Lx"]), int(params["Ly"]))
    hop_b = [[t, i, j] for i, j in bonds] + [[t, j, i] for i, j in bonds]
    hop_c = [[coef, i + nsite, j + nsite] for coef, i, j in hop_b]
    return make_operator([["+-", hop_b], ["+-", hop_c]], basis)


def double_occ_operator(nsite: int, basis):
    return make_operator([["nn", [[1.0 / nsite, i, i + nsite] for i in range(nsite)]]], basis)


def local_numsquare_operator(nsite: int, basis):
    terms = [[1.0 / nsite, i, i] for i in range(nsite)]
    terms += [[1.0 / nsite, i + nsite, i + nsite] for i in range(nsite)]
    return make_operator([["nn", terms]], basis)


def numsquare_operator(nsite: int, basis, layer: str):
    offset = 0 if layer == "b" else nsite
    terms = [[1.0, offset + i, offset + j] for i in range(nsite) for j in range(nsite)]
    return make_operator([["nn", terms]], basis)


def onsite_n2_operator(nsite: int, basis, layer: str):
    offset = 0 if layer == "b" else nsite
    terms = [[1.0, offset + i, offset + i] for i in range(nsite)]
    return make_operator([["nn", terms]], basis)


def pair_equal_operator(nsite: int, basis):
    pair = [[1.0 / nsite, i, i + nsite] for i in range(nsite)]
    return make_operator([["++", pair], ["--", pair]], basis)


def _npbc(nr: int, length: int) -> int:
    return nr % length + 1


def _cell_coord_1based(index: int, lx: int) -> tuple[int, int]:
    return index % lx + 1, index // lx + 1


def _cell_index_1based(x: int, y: int, lx: int) -> int:
    return (y - 1) * lx + (x - 1)


def is_k_geometry_compatible(lx: int, ly: int) -> bool:
    return int(lx) % 3 == 0 and int(ly) % 3 == 0


def structure_k_index(lx: int, ly: int) -> int:
    if not is_k_geometry_compatible(lx, ly):
        raise ValueError("K-point structure factors require Lx and Ly to be multiples of 3")
    return _cell_index_1based(2 * lx // 3 + 1, ly // 3 + 1, lx)


def _aimj_vector(index: int, lx: int) -> tuple[float, float]:
    x, y = _cell_coord_1based(index, lx)
    return (
        float(x - 1) * A1[0] + float(y - 1) * A2[0],
        float(x - 1) * A1[1] + float(y - 1) * A2[1],
    )


def _k_vector(index: int, lx: int, ly: int) -> tuple[float, float]:
    x, y = _cell_coord_1based(index, lx)
    return (
        float(x - 1) * B1[0] / float(lx) + float(y - 1) * B2[0] / float(ly),
        float(x - 1) * B1[1] / float(lx) + float(y - 1) * B2[1] / float(ly),
    )


def triangular_phase_matrix(lx: int, ly: int, k_index: int) -> np.ndarray:
    nsite = lx * ly
    kx, ky = _k_vector(k_index, lx, ly)
    phases = np.empty((nsite, nsite), dtype=np.complex128)
    for i in range(nsite):
        ix, iy = _cell_coord_1based(i, lx)
        for j in range(nsite):
            jx, jy = _cell_coord_1based(j, lx)
            imjx = _npbc(ix - jx, lx)
            imjy = _npbc(iy - jy, ly)
            imj = _cell_index_1based(imjx, imjy, lx)
            rx, ry = _aimj_vector(imj, lx)
            phases[i, j] = np.exp(1j * (kx * rx + ky * ry))
    return phases


def sf_structure_operator(params: dict[str, Any], basis):
    lx = int(params["Lx"])
    ly = int(params["Ly"])
    nsite = lx * ly
    norm = 1.0 / float(nsite * nsite)
    phases = triangular_phase_matrix(lx, ly, structure_k_index(lx, ly))
    terms_b = [[norm * phases[i, j], i, j] for i in range(nsite) for j in range(nsite)]
    terms_c = [[coef, i + nsite, j + nsite] for coef, i, j in terms_b]
    return make_operator(
        [["+-", terms_b], ["+-", terms_c]],
        basis,
        dtype=np.complex128,
        check_herm=False,
    )


def dw_structure_operator(params: dict[str, Any], basis):
    lx = int(params["Lx"])
    ly = int(params["Ly"])
    nsite = lx * ly
    norm = 1.0 / float(nsite * nsite)
    phases = triangular_phase_matrix(lx, ly, structure_k_index(lx, ly))
    terms = []
    for i in range(nsite):
        for j in range(nsite):
            coef = norm * phases[i, j]
            terms.append([coef, i, j])
            terms.append([coef, i + nsite, j + nsite])
            terms.append([coef, i, j + nsite])
            terms.append([coef, i + nsite, j])
    return make_operator(
        [["nn", terms]],
        basis,
        dtype=np.complex128,
        check_herm=False,
    )


def psf_gamma_operator(nsite: int, basis):
    norm = 1.0 / float(nsite * nsite)
    terms = [
        [norm, i, i + nsite, j + nsite, j]
        for i in range(nsite)
        for j in range(nsite)
    ]
    return make_operator(
        [["++--", terms]],
        basis,
        dtype=np.complex128,
        check_herm=False,
    )


def observable_values(params: dict[str, Any], basis, evals: np.ndarray, evecs: np.ndarray) -> dict[str, float]:
    beta = float(params["beta"])
    nsite = int(params["Lx"]) * int(params["Ly"])

    def avg(operator, name: str) -> float:
        op_diag = operator_diag_in_eigenbasis(operator, evecs)
        value, _, _ = thermal_average(evals, op_diag, beta)
        return real_physical_average(value, name)

    num_up = avg(number_operator(nsite, basis, "b"), "num_up")
    num_do = avg(number_operator(nsite, basis, "c"), "num_do")
    total_ne = num_up + num_do
    density_total = total_ne / nsite
    kinetic_total = avg(kinetic_operator(params, basis), "kinetic_total")
    double_occ = avg(double_occ_operator(nsite, basis), "doubleOcc")
    onsite_n2_up = avg(onsite_n2_operator(nsite, basis, "b"), "onsite_n2_up")
    onsite_n2_do = avg(onsite_n2_operator(nsite, basis, "c"), "onsite_n2_do")
    pair_equal = avg(pair_equal_operator(nsite, basis), "pair_equal")
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
        "numsquare_up": avg(numsquare_operator(nsite, basis, "b"), "numsquare_up"),
        "numsquare_do": avg(numsquare_operator(nsite, basis, "c"), "numsquare_do"),
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
        "S_PSF_Gamma": avg(psf_gamma_operator(nsite, basis), "S_PSF_Gamma"),
        "same_flavor_pair": 0.5 * ((onsite_n2_up - num_up) + (onsite_n2_do - num_do)) / nsite,
    }
    if is_k_geometry_compatible(int(params["Lx"]), int(params["Ly"])):
        values.update(
            {
                "S_SF_K": avg(sf_structure_operator(params, basis), "S_SF_K"),
                "S_DW_K": avg(dw_structure_operator(params, basis), "S_DW_K"),
            }
        )
    return values


def decode_state(state: int, length: int, sps: int) -> list[int]:
    occ: list[int] = []
    value = int(state)
    for _ in range(length):
        occ.append(value % sps)
        value //= sps
    return occ


def fock_diag_in_eigenbasis(fock_diag: np.ndarray, evecs: np.ndarray) -> np.ndarray:
    weights = np.abs(evecs) ** 2
    return np.einsum("ia,i->a", weights, fock_diag, optimize=True)


def cutoff_diagnostics(basis, evecs: np.ndarray, evals: np.ndarray, params: dict[str, Any]) -> dict[str, float | int]:
    """Return ordering-invariant cutoff diagnostics from integer Fock encodings.

    If a future diagnostic becomes site-resolved or flavor-resolved, do not reuse
    decode_state blindly; first add a self-test against basis.int_to_state() or
    use QuSpin operator expectation values.
    """
    length = 2 * int(params["Lx"]) * int(params["Ly"])
    sps = int(params["nmax"]) + 1
    occupations = np.array([decode_state(int(state), length, sps) for state in basis.states], dtype=np.int16)
    total_particles = occupations.sum(axis=1)
    ncut = params.get("ncut")

    total_boundary = np.zeros_like(total_particles, dtype=np.float64)
    if ncut is not None:
        total_boundary = (total_particles == int(ncut)).astype(np.float64)
    local_boundary = (occupations == int(params["nmax"])).any(axis=1).astype(np.float64)

    total_boundary_diag = fock_diag_in_eigenbasis(total_boundary, evecs)
    local_boundary_diag = fock_diag_in_eigenbasis(local_boundary, evecs)
    total_weight, _, _ = thermal_average(evals, total_boundary_diag, float(params["beta"]))
    local_weight, _, _ = thermal_average(evals, local_boundary_diag, float(params["beta"]))
    return {
        "max_total_particles_in_basis": int(total_particles.max()),
        "total_cutoff_boundary_weight": float(total_weight),
        "local_nmax_boundary_weight": float(local_weight),
    }


def wrap_observables(values: dict[str, float]) -> dict[str, dict[str, float | str]]:
    return {
        key: {
            "value": float(values[key]),
            "normalization": OBSERVABLE_NORMALIZATION[key],
        }
        for key in OBSERVABLE_KEYS
        if key in values
    }


def result_payload(
    params: dict[str, Any],
    basis_size: int,
    evals: np.ndarray,
    observables: dict[str, float],
    diagnostics: dict[str, float | int],
    diagonalization_method: str,
) -> dict[str, Any]:
    bonds = [[int(i), int(j)] for i, j in triangular_bonds(int(params["Lx"]), int(params["Ly"]))]
    return {
        "case": params["case"],
        "parameters": {key: params[key] for key in PARAMETER_KEYS},
        "cutoffs": {
            "nmax": int(params["nmax"]),
            "ncut": params["ncut"],
            "basis_cap": int(params["basis_cap"]),
        },
        "basis_size": int(basis_size),
        "diagonalization_method": diagonalization_method,
        "energy_shift": float(evals.min()),
        "free_energy": free_energy(evals, float(params["beta"])),
        "observables": wrap_observables(observables),
        "cutoff_diagnostics": diagnostics,
        "bond_list": bonds,
    }


def run_ed(
    params: dict[str, Any],
    basis_cap: int | None = None,
    allow_large_basis: bool = False,
    dense_memory_cap_gib: float | None = DEFAULT_DENSE_MEMORY_CAP_GIB,
) -> dict[str, Any]:
    normalized = normalize_params(params)
    if basis_cap is not None:
        normalized["basis_cap"] = int(basis_cap)

    estimated_basis_size = basis_size_estimate(
        2 * int(normalized["Lx"]) * int(normalized["Ly"]),
        int(normalized["nmax"]),
        normalized["ncut"],
    )
    if estimated_basis_size > int(normalized["basis_cap"]) and not allow_large_basis:
        raise SystemExit(
            f"estimated basis.Ns={estimated_basis_size} exceeds basis_cap={normalized['basis_cap']}; "
            "use --allow-large-basis for manual references"
        )
    assert_dense_memory_budget(estimated_basis_size, dense_memory_cap_gib=dense_memory_cap_gib)

    basis = build_basis(int(normalized["Lx"]), int(normalized["Ly"]), int(normalized["nmax"]), normalized["ncut"])
    if basis.Ns != estimated_basis_size:
        raise SystemExit(
            f"estimated basis.Ns={estimated_basis_size} disagrees with QuSpin basis.Ns={basis.Ns}"
        )
    if basis.Ns > int(normalized["basis_cap"]) and not allow_large_basis:
        raise SystemExit(
            f"basis.Ns={basis.Ns} exceeds basis_cap={normalized['basis_cap']}; "
            "use --allow-large-basis for manual references"
        )
    assert_dense_memory_budget(basis.Ns, dense_memory_cap_gib=dense_memory_cap_gib)

    hamiltonian_mu = build_hamiltonian(normalized, basis)
    evals, evecs = hamiltonian_mu.eigh()
    observables = observable_values(normalized, basis, evals, evecs)
    diagnostics = cutoff_diagnostics(basis, evecs, evals, normalized)
    diagnostics.update(
        dense_memory_diagnostics(basis.Ns, dense_memory_cap_gib=dense_memory_cap_gib)
    )
    return result_payload(normalized, basis.Ns, evals, observables, diagnostics, "dense_H_eigh")


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
        default=DEFAULT_DENSE_MEMORY_CAP_GIB,
        help=(
            "Refuse dense ED when estimated peak RSS exceeds this cap; "
            "set 0 to disable after reserving sufficient memory"
        ),
    )
    parser.add_argument("--print-bonds", action="store_true")
    args = parser.parse_args()

    params = load_params(args.params)
    result = run_ed(
        params,
        basis_cap=args.basis_cap,
        allow_large_basis=args.allow_large_basis,
        dense_memory_cap_gib=args.dense_memory_cap_gib,
    )
    if args.print_bonds:
        print(json.dumps(result["bond_list"]), file=sys.stderr)
    write_result(result, args.output)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
