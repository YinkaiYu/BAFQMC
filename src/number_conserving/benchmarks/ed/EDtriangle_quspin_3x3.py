#!/usr/bin/env python3
"""3x3 QuSpin exact-diagonalization reference schema helper.

This module is intentionally import-light: importing it must not require
QuSpin. The heavy ED implementation is isolated in ``run_ed_reference`` so
campaign tooling can import the schema helpers with the system Python.
"""

from __future__ import annotations

import argparse
import json
import math
import os
import tempfile
from datetime import datetime, timezone
from pathlib import Path
from typing import Any


SCHEMA_VERSION = 1
DEFAULT_DENSE_MEMORY_CAP_GIB = 12.0

ED_OBSERVABLES = (
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
    "numsquare_up",
    "numsquare_do",
)


def empty_weighted_observables() -> dict[str, float]:
    """Return zero-filled weighted observable accumulators."""

    return {name: 0.0 for name in ED_OBSERVABLES}


def normalize_weighted_observables(
    weighted: dict[str, float], partition: float
) -> dict[str, float]:
    """Normalize weighted observable sums by the grand-canonical partition."""

    if partition == 0.0:
        raise ValueError("partition must be nonzero")
    return {name: float(weighted.get(name, 0.0)) / partition for name in ED_OBSERVABLES}


def relative_tail(
    last_shell: dict[str, float],
    observables: dict[str, float],
    floor: float = 1e-12,
) -> dict[str, float]:
    """Return relative last-shell contributions for every ED observable."""

    if floor <= 0.0:
        raise ValueError("floor must be positive")
    return {
        name: abs(float(last_shell.get(name, 0.0)))
        / max(abs(float(observables.get(name, 0.0))), floor)
        for name in ED_OBSERVABLES
    }


def classify_status(
    parameters: dict[str, Any],
    relative_last_shell: dict[str, float],
    policy: dict[str, Any],
    completed_shells: int = 0,
    observables: dict[str, float] | None = None,
) -> str:
    """Classify ED convergence/cutoff status with explicit U1<0 semantics."""

    if bool(parameters.get("dry_run", False)):
        return "dry_run"

    u1 = float(parameters.get("U1", 0.0))
    tolerance = float(policy.get("tail_tolerance", 0.0))
    max_relative_tail = max(
        (float(relative_last_shell.get(name, math.inf)) for name in ED_OBSERVABLES),
        default=math.inf,
    )
    min_completed_shells = int(policy.get("min_completed_shells", 0))

    if (
        u1 >= 0.0
        and completed_shells >= min_completed_shells
        and max_relative_tail <= tolerance
    ):
        return "converged"

    if u1 < 0.0:
        density_total = math.inf
        if observables is not None:
            density_total = float(observables.get("density_total", math.inf))
        max_density_total = float(policy.get("max_density_total", -math.inf))
        if (
            density_total <= max_density_total
            and completed_shells >= min_completed_shells
        ):
            return "low_density_cutoff_accepted"

    return "incomplete"


def result_payload(
    *,
    parameters: dict[str, Any],
    observables: dict[str, float],
    last_shell_contribution: dict[str, float],
    relative_last_shell: dict[str, float],
    convergence_policy: dict[str, Any],
    status: str,
    completed_shells: int,
    timestamp_utc: str | None = None,
) -> dict[str, Any]:
    """Build a JSON-serializable ED result payload."""

    return {
        "schema_version": SCHEMA_VERSION,
        "parameters": _json_clean(parameters),
        "observables": _ordered_observable_dict(observables),
        "last_shell_contribution": _ordered_observable_dict(last_shell_contribution),
        "relative_last_shell": _ordered_observable_dict(relative_last_shell),
        "convergence_policy": _json_clean(convergence_policy),
        "status": str(status),
        "completed_shells": int(completed_shells),
        "timestamp_utc": timestamp_utc or _utc_timestamp(),
    }


def write_results(path: Path, payload: dict[str, Any]) -> None:
    """Write JSON results atomically through a temporary file."""

    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    fd, tmp_name = tempfile.mkstemp(
        prefix=f".{path.name}.", suffix=".tmp", dir=str(path.parent)
    )
    try:
        with os.fdopen(fd, "w", encoding="utf-8") as handle:
            json.dump(payload, handle, indent=2, sort_keys=True)
            handle.write("\n")
        os.replace(tmp_name, path)
    except Exception:
        try:
            os.unlink(tmp_name)
        except FileNotFoundError:
            pass
        raise


def run_ed_reference(
    parameters: dict[str, Any],
    policy: dict[str, Any],
    output_path: Path | None = None,
    dense_memory_cap_gib: float | None = DEFAULT_DENSE_MEMORY_CAP_GIB,
) -> dict[str, Any]:
    """Run the 3x3 QuSpin ED reference calculation.

    QuSpin and Numba are imported here, not at module import time, so campaign
    tooling can import this module with the system Python for schema handling.
    """

    _validate_parameters(parameters)
    runtime = _load_quspin_runtime()
    runtime["dense_memory_cap_gib"] = dense_memory_cap_gib
    return _run_ed_reference_runtime(
        runtime,
        parameters=parameters,
        policy=policy,
        output_path=output_path,
    )


def _load_quspin_runtime() -> dict[str, Any]:
    """Load QuSpin/Numba and compile the 64-bit user_basis callbacks."""

    try:
        import numpy as np
        from numba import carray, cfunc, float64, uint64
        from quspin.basis.user import (
            map_sig_64,
            next_state_sig_64,
            op_sig_64,
            user_basis,
        )
        from quspin.operators import hamiltonian
    except ImportError as exc:
        raise RuntimeError(
            "QuSpin is required for 3x3 ED execution. Use "
            "the Python interpreter in your QuSpin environment, or install QuSpin in the "
            "active environment."
        ) from exc

    @cfunc(op_sig_64, locals=dict(n=uint64, new_state=uint64, fac=float64))
    def op_two_species(op_struct_ptr, op_str, ind, N, args):
        arr = carray(args, 2 * N)
        sps = arr[:N]
        full_factor = arr[N : 2 * N]
        op_struct = carray(op_struct_ptr, 1)[0]
        local_sps = sps[ind]
        base = full_factor[ind]
        n = (op_struct.state // base) % local_sps
        new_state = op_struct.state
        fac = 1.0
        if op_str == ord("+"):
            if n < local_sps - 1:
                fac = math.sqrt(n + 1)
                new_state = op_struct.state + base
            else:
                fac = 0.0
        elif op_str == ord("-"):
            if n > 0:
                fac = math.sqrt(n)
                new_state = op_struct.state - base
            else:
                fac = 0.0
        elif op_str == ord("n"):
            fac = float(n)
        else:
            fac = 0.0
        op_struct.state = new_state
        op_struct.matrix_ele *= fac
        return 0

    @cfunc(
        map_sig_64,
        locals=dict(
            L_total=uint64,
            offset1=uint64,
            offset2=uint64,
            offset3=uint64,
            Lq=uint64,
            i=uint64,
            inv=uint64,
            base=uint64,
            digit=uint64,
            sps_val=uint64,
            new_state=uint64,
        ),
    )
    def translation_map_2species(state, N_dummy, sign_ptr, args):
        L_total = args[0]
        offset1 = 1
        offset2 = offset1 + L_total
        offset3 = offset2 + L_total
        Lq = L_total // 2
        new_state = 0

        for i in range(Lq):
            inv = args[offset3 + i]
            base = args[offset2 + inv]
            sps_val = args[offset1 + inv]
            digit = (state // base) % sps_val
            new_state += digit * args[offset2 + i]

        for i in range(Lq, L_total):
            inv = args[offset3 + i]
            base = args[offset2 + inv]
            sps_val = args[offset1 + inv]
            digit = (state // base) % sps_val
            new_state += digit * args[offset2 + i]

        return new_state

    @cfunc(next_state_sig_64, locals=dict(next_state=uint64))
    def next_state(s, counter, N, args):
        next_state = args[counter + 1]
        return next_state

    return {
        "np": np,
        "hamiltonian": hamiltonian,
        "user_basis": user_basis,
        "op_two_species": op_two_species,
        "translation_map_2species": translation_map_2species,
        "next_state": next_state,
    }


def _run_ed_reference_runtime(
    runtime: dict[str, Any],
    *,
    parameters: dict[str, Any],
    policy: dict[str, Any],
    output_path: Path | None,
) -> dict[str, Any]:
    lx = int(parameters["Lx"])
    ly = int(parameters["Ly"])
    lq = lx * ly
    t = float(parameters.get("t", 1.0))
    u1 = float(parameters["U1"])
    u2 = float(parameters["U2"])
    beta = float(parameters["beta"])
    mu = float(parameters["mu"])
    max_total_particles = int(parameters["max_total_particles"])

    partition_total = 0.0
    weighted_total = empty_weighted_observables()
    final_payload: dict[str, Any] | None = None

    for total_ne in range(max_total_particles + 1):
        shell_partition = 0.0
        shell_weighted = empty_weighted_observables()

        if total_ne == 0:
            shell_partition = 1.0
        else:
            for ne_b in range(total_ne + 1):
                ne_c = total_ne - ne_b
                block_weighted, block_partition = _grand_shell_pair(
                    runtime,
                    lx=lx,
                    ly=ly,
                    t=t,
                    u1=u1,
                    u2=u2,
                    beta=beta,
                    mu=mu,
                    ne_b=ne_b,
                    ne_c=ne_c,
                )
                shell_partition += block_partition
                _add_weighted_in_place(shell_weighted, block_weighted)

        partition_total += float(shell_partition)
        _add_weighted_in_place(weighted_total, shell_weighted)

        observables = normalize_weighted_observables(weighted_total, partition_total)
        last_shell = normalize_weighted_observables(shell_weighted, partition_total)
        tails = relative_tail(last_shell, observables)
        completed_shells = int(total_ne)
        status = classify_status(
            parameters,
            tails,
            policy,
            completed_shells=completed_shells,
            observables=observables,
        )
        final_payload = result_payload(
            parameters=parameters,
            observables=observables,
            last_shell_contribution=last_shell,
            relative_last_shell=tails,
            convergence_policy=policy,
            status=status,
            completed_shells=completed_shells,
        )
        if output_path is not None:
            write_results(Path(output_path), final_payload)
        if _should_stop_after_shell(status, policy, completed_shells):
            break

    if final_payload is None:
        raise RuntimeError("ED shell loop did not produce a result payload")
    return final_payload


def _should_stop_after_shell(
    status: str, policy: dict[str, Any], completed_shells: int
) -> bool:
    if not bool(policy.get("stop_when_reliable", True)):
        return False
    min_completed_shells = int(policy.get("min_completed_shells", 0))
    if int(completed_shells) < min_completed_shells:
        return False
    return status in {"converged", "low_density_cutoff_accepted"}


def _grand_shell_pair(
    runtime: dict[str, Any],
    *,
    lx: int,
    ly: int,
    t: float,
    u1: float,
    u2: float,
    beta: float,
    mu: float,
    ne_b: int,
    ne_c: int,
) -> tuple[dict[str, float], float]:
    shell_weighted = empty_weighted_observables()
    shell_partition = 0.0

    for kx in range(lx):
        for ky in range(ly):
            block_weighted, block_partition = _grand_momentum_block(
                runtime,
                lx=lx,
                ly=ly,
                t=t,
                u1=u1,
                u2=u2,
                beta=beta,
                mu=mu,
                ne_b=ne_b,
                ne_c=ne_c,
                kx=kx,
                ky=ky,
            )
            shell_partition += block_partition
            _add_weighted_in_place(shell_weighted, block_weighted)

    return shell_weighted, shell_partition


def assert_dense_memory_budget(dimension: int, dense_memory_cap_gib: float | None = DEFAULT_DENSE_MEMORY_CAP_GIB) -> float:
    """Conservative four-complex-matrix estimate; the cap is not a particle cutoff."""
    estimate = 4 * 16 * int(dimension)**2 / 1024**3
    if dense_memory_cap_gib is not None and float(dense_memory_cap_gib) > 0 and estimate > float(dense_memory_cap_gib):
        raise MemoryError(
            f"Dense ED block dimension {dimension} estimates {estimate:.2f} GiB, exceeding "
            f"--dense-memory-cap-gib={float(dense_memory_cap_gib):g}. The last completed "
            "particle shell remains in the output file. Increase this cap only on a "
            "machine with sufficient RAM; physical cutoffs were not changed.")
    return estimate


def _grand_momentum_block(
    runtime: dict[str, Any],
    *,
    lx: int,
    ly: int,
    t: float,
    u1: float,
    u2: float,
    beta: float,
    mu: float,
    ne_b: int,
    ne_c: int,
    kx: int,
    ky: int,
) -> tuple[dict[str, float], float]:
    np = runtime["np"]
    lq = lx * ly
    total_ne = ne_b + ne_c
    basis = _build_two_species_user_basis(runtime, lx, ly, ne_b, ne_c, kx, ky)
    if int(basis.Ns) == 0:
        return empty_weighted_observables(), 0.0

    hamiltonian = runtime["hamiltonian"]
    static = _hamiltonian_static(lx, ly, t, u1, u2)
    H = hamiltonian(
        static,
        [],
        basis=basis,
        dtype=np.complex128,
        check_symm=False,
        check_herm=False,
        check_pcon=False,
    )
    assert_dense_memory_budget(H.Ns, runtime.get("dense_memory_cap_gib", DEFAULT_DENSE_MEMORY_CAP_GIB))
    E, V = H.eigh()
    weights = np.exp(-beta * (np.real(E) - mu * total_ne))
    partition = float(np.sum(weights))

    weighted = empty_weighted_observables()
    weighted["density_total"] = (total_ne / lq) * partition
    weighted["num_up"] = ne_b * partition
    weighted["num_do"] = ne_c * partition
    weighted["numsquare_up"] = ne_b * ne_b * partition
    weighted["numsquare_do"] = ne_c * ne_c * partition

    operators = {
        "total_kinetic": [["+-", _hopping_two_species(lx, ly, t)]],
        "doubleOcc": [["nn", [[1.0 / lq, i, i + lq] for i in range(lq)]]],
        "onsite_n2_up": [["nn", [[1.0, i, i] for i in range(lq)]]],
        "onsite_n2_do": [["nn", [[1.0, i + lq, i + lq] for i in range(lq)]]],
        "S_SF_K": [["+-", _sf_k_terms(lx, ly)]],
        "S_PSF_Gamma": [["++--", _psf_gamma_terms(lx, ly)]],
        "S_DW_K": [["nn", _dw_k_terms(lx, ly)]],
    }
    for name, op_static in operators.items():
        op = hamiltonian(
            op_static,
            [],
            basis=basis,
            dtype=np.complex128,
            check_symm=False,
            check_herm=False,
            check_pcon=False,
        )
        diagonal = op.matrix_ele(V, V, diagonal=True)
        weighted[name] = _real_scalar(np.sum(weights * diagonal))

    weighted["interaction_energy_density"] = (
        (u1 + u2) * (weighted["onsite_n2_up"] + weighted["onsite_n2_do"]) / lq
        + 2.0 * (u1 - u2) * weighted["doubleOcc"]
    )
    weighted["energy_density"] = (
        weighted["total_kinetic"] / lq + weighted["interaction_energy_density"]
    )
    return weighted, partition


def _generate_two_species_basis(
    np: Any, lx: int, ly: int, ne_b: int, ne_c: int
) -> tuple[Any, Any, Any]:
    lq = lx * ly
    total_sites = 2 * lq
    sps_b = ne_b + 1
    sps_c = ne_c + 1
    sps_arr = np.array([sps_b] * lq + [sps_c] * lq, dtype=np.uint64)
    full_factor = np.empty(total_sites, dtype=np.uint64)
    full_factor[0] = 1
    for i in range(1, total_sites):
        full_factor[i] = full_factor[i - 1] * sps_arr[i - 1]

    def gen_configs(n_sites: int, particles: int, local_dim: int) -> list[list[int]]:
        configs: list[list[int]] = []
        config = [0] * n_sites

        def rec(site: int, remaining: int) -> None:
            if site == n_sites - 1:
                if remaining < local_dim:
                    config[site] = remaining
                    configs.append(config.copy())
                return
            for occupation in range(min(local_dim, remaining + 1)):
                config[site] = occupation
                rec(site + 1, remaining - occupation)

        rec(0, particles)
        return configs

    valid_states: list[int] = []
    b_configs = gen_configs(lq, ne_b, sps_b)
    c_configs = gen_configs(lq, ne_c, sps_c)
    for b_conf in b_configs:
        state_b = 0
        for i, occupation in enumerate(b_conf):
            state_b += int(occupation) * int(full_factor[i])
        for c_conf in c_configs:
            state = state_b
            for i, occupation in enumerate(c_conf):
                state += int(occupation) * int(full_factor[lq + i])
            valid_states.append(state)

    return np.array(valid_states, dtype=np.uint64), sps_arr, full_factor


def _build_two_species_user_basis(
    runtime: dict[str, Any],
    lx: int,
    ly: int,
    ne_b: int,
    ne_c: int,
    kx: int,
    ky: int,
) -> Any:
    np = runtime["np"]
    lq = lx * ly
    total_sites = 2 * lq
    valid_states, sps_arr, full_factor = _generate_two_species_basis(
        np, lx, ly, ne_b, ne_c
    )
    tx, ty = _translation_permutations(np, lx, ly)
    tx_full = np.concatenate([tx, tx + lq])
    ty_full = np.concatenate([ty, ty + lq])
    inv_tx = np.empty(total_sites, dtype=np.uint64)
    inv_ty = np.empty(total_sites, dtype=np.uint64)
    for i in range(total_sites):
        inv_tx[tx_full[i]] = i
        inv_ty[ty_full[i]] = i

    args_x = np.concatenate(
        [np.array([total_sites], dtype=np.uint64), sps_arr, full_factor, inv_tx]
    )
    args_y = np.concatenate(
        [np.array([total_sites], dtype=np.uint64), sps_arr, full_factor, inv_ty]
    )

    class FunctionWrapper:
        def __init__(self, basis_arr: Any) -> None:
            self.basis = basis_arr

        def get_s0_pcon(self, N: int, Np: Any) -> int:
            return int(self.basis[0])

        def get_Ns_pcon(self, N: int, Np: Any) -> int:
            return int(self.basis.size)

    wrapper = FunctionWrapper(valid_states)
    pcon_dict = {
        "Np": (),
        "next_state": runtime["next_state"],
        "next_state_args": valid_states,
        "get_Ns_pcon": wrapper.get_Ns_pcon,
        "get_s0_pcon": wrapper.get_s0_pcon,
    }
    op_args = np.concatenate([sps_arr, full_factor])

    return runtime["user_basis"](
        np.uint64,
        total_sites,
        op_dict={"op": runtime["op_two_species"], "op_args": op_args},
        allowed_ops=set("+-n"),
        sps=int(np.max(sps_arr)),
        pcon_dict=pcon_dict,
        Tx_block=(runtime["translation_map_2species"], lx, kx, args_x),
        Ty_block=(runtime["translation_map_2species"], ly, ky, args_y),
    )


def _translation_permutations(np: Any, lx: int, ly: int) -> tuple[Any, Any]:
    tx = np.empty(lx * ly, dtype=np.uint64)
    ty = np.empty(lx * ly, dtype=np.uint64)
    for y in range(ly):
        for x in range(lx):
            site = _site_index(x, y, lx)
            tx[site] = _site_index((x + 1) % lx, y, lx)
            ty[site] = _site_index(x, (y + 1) % ly, lx)
    return tx, ty


def _hamiltonian_static(
    lx: int, ly: int, t: float, u1: float, u2: float
) -> list[list[Any]]:
    lq = lx * ly
    return [
        ["+-", _hopping_two_species(lx, ly, t)],
        ["nn", [[2.0 * (u1 - u2), i, i + lq] for i in range(lq)]],
        ["nn", [[u1 + u2, i, i] for i in range(2 * lq)]],
    ]


def _hopping_two_species(lx: int, ly: int, t: float) -> list[list[float | int]]:
    lq = lx * ly
    hopping: list[list[float | int]] = []
    directions = ((1, 0), (0, 1), (-1, 1), (-1, 0), (0, -1), (1, -1))
    for offset in (0, lq):
        for y in range(ly):
            for x in range(lx):
                i = _site_index(x, y, lx) + offset
                for dx, dy in directions:
                    j = _site_index((x + dx) % lx, (y + dy) % ly, lx) + offset
                    hopping.append([t, i, j])
    return hopping


def _sf_k_terms(lx: int, ly: int) -> list[list[complex | int]]:
    lq = lx * ly
    norm = 1.0 / float(lq * lq)
    terms: list[list[complex | int]] = []
    for i in range(lq):
        for j in range(lq):
            coeff = norm * _k_phase(lx, ly, i, j)
            terms.append([coeff, i, j])
            terms.append([coeff, i + lq, j + lq])
    return terms


def _psf_gamma_terms(lx: int, ly: int) -> list[list[float | int]]:
    lq = lx * ly
    norm = 1.0 / float(lq * lq)
    return [[norm, i, i + lq, j + lq, j] for i in range(lq) for j in range(lq)]


def _dw_k_terms(lx: int, ly: int) -> list[list[complex | int]]:
    lq = lx * ly
    norm = 1.0 / float(lq * lq)
    terms: list[list[complex | int]] = []
    for i in range(lq):
        for j in range(lq):
            coeff = norm * _k_phase(lx, ly, i, j)
            terms.append([coeff, i, j])
            terms.append([coeff, i + lq, j + lq])
            terms.append([coeff, i, j + lq])
            terms.append([coeff, i + lq, j])
    return terms


def _k_phase(lx: int, ly: int, i: int, j: int) -> complex:
    if lx % 3 != 0 or ly % 3 != 0:
        return 0.0 + 0.0j
    xi, yi = _site_coords(i, lx)
    xj, yj = _site_coords(j, lx)
    dx = (xi - xj) % lx
    dy = (yi - yj) % ly
    kx_cell = 2 * lx // 3
    ky_cell = ly // 3
    angle = 2.0 * math.pi * (
        float(kx_cell * dx) / float(lx) + float(ky_cell * dy) / float(ly)
    )
    return complex(math.cos(angle), math.sin(angle))


def _pbc_1based(value: int, size: int) -> int:
    if value > size:
        return value - size
    if value < 1:
        return value + size
    return value


def _site_index(x: int, y: int, lx: int) -> int:
    return y * lx + x


def _site_coords(site: int, lx: int) -> tuple[int, int]:
    return site % lx, site // lx


def _add_weighted_in_place(total: dict[str, float], increment: dict[str, float]) -> None:
    for name in ED_OBSERVABLES:
        total[name] = float(total.get(name, 0.0)) + float(increment.get(name, 0.0))


def _real_scalar(value: Any) -> float:
    return float(complex(value).real)


def dry_run_payload(
    parameters: dict[str, Any], policy: dict[str, Any]
) -> dict[str, Any]:
    """Build an explicit schema-only payload for CLI validation."""

    observables = empty_weighted_observables()
    last_shell = empty_weighted_observables()
    tails = relative_tail(last_shell, observables)
    return result_payload(
        parameters={**parameters, "dry_run": True},
        observables=observables,
        last_shell_contribution=last_shell,
        relative_last_shell=tails,
        convergence_policy=policy,
        status="dry_run",
        completed_shells=0,
    )


def main(argv: list[str] | None = None) -> int:
    args = _parse_args(argv)
    parameters, policy = _load_inputs(args)
    _validate_parameters(parameters)

    if args.dry_run:
        payload = dry_run_payload(parameters, policy)
    else:
        payload = run_ed_reference(parameters, policy, output_path=args.output,
                                   dense_memory_cap_gib=args.dense_memory_cap_gib)

    write_results(args.output, payload)
    return 0


def _parse_args(argv: list[str] | None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Build or run a 3x3 QuSpin ED reference result payload."
    )
    parser.add_argument(
        "params_json",
        nargs="?",
        type=Path,
        help="JSON parameters file. May contain parameters/convergence_policy keys.",
    )
    parser.add_argument("--params", type=Path, help="JSON parameters file.")
    parser.add_argument("--Lx", type=int)
    parser.add_argument("--Ly", type=int)
    parser.add_argument("--beta", type=float)
    parser.add_argument("--mu", type=float)
    parser.add_argument("--U1", type=float)
    parser.add_argument("--U2", type=float)
    parser.add_argument("--t", type=float, default=None)
    parser.add_argument("--max-total-particles", type=int, dest="max_total_particles")
    parser.add_argument("--output", type=Path, default=Path("results.json"))
    parser.add_argument("--dense-memory-cap-gib", type=float, default=DEFAULT_DENSE_MEMORY_CAP_GIB,
                        help="Refuse dense blocks estimated above this RAM budget; <=0 disables the guard.")
    parser.add_argument("--tail-tolerance", type=float, default=None)
    parser.add_argument("--max-density-total", type=float, default=None)
    parser.add_argument("--min-completed-shells", type=int, default=None)
    parser.add_argument(
        "--no-stop-when-reliable",
        action="store_false",
        dest="stop_when_reliable",
        default=None,
        help="Continue through max_total_particles even after a reliable shell.",
    )
    parser.add_argument(
        "--dry-run",
        action="store_true",
        help="Write an explicit schema-only payload without running ED.",
    )
    return parser.parse_args(argv)


def _load_inputs(args: argparse.Namespace) -> tuple[dict[str, Any], dict[str, Any]]:
    parameter_file = args.params or args.params_json
    parameters: dict[str, Any] = {}
    policy: dict[str, Any] = {}

    if parameter_file is not None:
        data = json.loads(parameter_file.read_text(encoding="utf-8"))
        if not isinstance(data, dict):
            raise ValueError("JSON parameters file must contain an object")
        if "parameters" in data:
            raw_parameters = data["parameters"]
        else:
            raw_parameters = {
                key: value for key, value in data.items() if key != "convergence_policy"
            }
        raw_policy = data.get("convergence_policy", {})
        if not isinstance(raw_parameters, dict):
            raise ValueError("parameters must be a JSON object")
        if not isinstance(raw_policy, dict):
            raise ValueError("convergence_policy must be a JSON object")
        parameters.update(raw_parameters)
        policy.update(raw_policy)

    if "max_total_particles" in policy and "max_total_particles" not in parameters:
        parameters["max_total_particles"] = policy.pop("max_total_particles")

    for name in ("Lx", "Ly", "beta", "mu", "U1", "U2", "t", "max_total_particles"):
        value = getattr(args, name)
        if value is not None:
            parameters[name] = value

    for arg_name, policy_name in (
        ("tail_tolerance", "tail_tolerance"),
        ("max_density_total", "max_density_total"),
        ("min_completed_shells", "min_completed_shells"),
    ):
        value = getattr(args, arg_name)
        if value is not None:
            policy[policy_name] = value
    if args.stop_when_reliable is not None:
        policy["stop_when_reliable"] = bool(args.stop_when_reliable)

    parameters.setdefault("t", 1.0)
    policy.setdefault("tail_tolerance", 1.0e-3)
    policy.setdefault("max_density_total", 1.0)
    policy.setdefault("min_completed_shells", 2)
    policy.setdefault("stop_when_reliable", True)
    return parameters, policy


def _validate_parameters(parameters: dict[str, Any]) -> None:
    missing = [
        name
        for name in ("Lx", "Ly", "beta", "mu", "U1", "U2", "max_total_particles")
        if name not in parameters
    ]
    if missing:
        raise ValueError(f"Missing required parameters: {', '.join(missing)}")

    lx = int(parameters["Lx"])
    ly = int(parameters["Ly"])
    if (lx, ly) != (3, 3):
        raise ValueError("EDtriangle_quspin_3x3.py only supports Lx=Ly=3")
    if float(parameters["beta"]) <= 0.0:
        raise ValueError("beta must be positive")
    if int(parameters["max_total_particles"]) < 0:
        raise ValueError("max_total_particles must be nonnegative")


def _ordered_observable_dict(values: dict[str, float]) -> dict[str, float]:
    return {name: float(values.get(name, 0.0)) for name in ED_OBSERVABLES}


def _utc_timestamp() -> str:
    return datetime.now(timezone.utc).isoformat(timespec="seconds").replace(
        "+00:00", "Z"
    )


def _json_clean(value: Any) -> Any:
    try:
        json.dumps(value)
    except TypeError:
        if isinstance(value, dict):
            return {str(key): _json_clean(item) for key, item in value.items()}
        if isinstance(value, (list, tuple)):
            return [_json_clean(item) for item in value]
        return str(value)
    return value


if __name__ == "__main__":
    raise SystemExit(main())
