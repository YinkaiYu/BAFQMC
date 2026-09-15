#!/usr/bin/env python3
"""Exact diagonalization reference for the 3x2 two-flavor boson benchmark.

The script reads ``params.txt`` in the current directory and writes:

- ``results.txt`` with the historical first four lines:
  total_NE, total_kinetic, last-shell total_NE contribution, last-shell kinetic
  contribution.
- ``results.json`` with all scalar observables used by the benchmark suite.

The implementation uses dense fixed-particle-number blocks. This keeps the ED
tool runnable on the repository's small benchmark cases without optional QuSpin
dependencies.
"""

from __future__ import annotations

import json
import logging
import math
import time
from functools import lru_cache
from pathlib import Path

import numpy as np


SCALAR_OBSERVABLES = (
    "total_NE",
    "total_kinetic",
    "doubleOcc",
    "squareOcc",
    "numsquare_up",
    "numsquare_do",
)

logging.basicConfig(
    filename="calc.log", level=logging.INFO, format="%(asctime)s - %(message)s"
)


def read_params(filename: str = "params.txt") -> dict[str, float]:
    params: dict[str, float] = {}
    with open(filename, "r", encoding="utf-8") as handle:
        for line in handle:
            stripped = line.strip()
            if not stripped or stripped.startswith("#"):
                continue
            key, value = stripped.split("=")
            params[key.strip()] = float(value.strip())
    return params


@lru_cache(maxsize=None)
def occupation_configs(n_sites: int, n_particles: int) -> tuple[tuple[int, ...], ...]:
    configs: list[tuple[int, ...]] = []
    current = [0] * n_sites

    def visit(site: int, remaining: int) -> None:
        if site == n_sites - 1:
            current[site] = remaining
            configs.append(tuple(current))
            return
        for count in range(remaining + 1):
            current[site] = count
            visit(site + 1, remaining - count)

    visit(0, n_particles)
    return tuple(configs)


@lru_cache(maxsize=None)
def directed_triangular_bonds(lx: int, ly: int) -> tuple[tuple[int, int], ...]:
    def site_index(x: int, y: int) -> int:
        return y * lx + x

    bonds: list[tuple[int, int]] = []
    for y in range(ly):
        for x in range(lx):
            source = site_index(x, y)
            neighbors = (
                site_index((x + 1) % lx, y),
                site_index(x, (y + 1) % ly),
                site_index((x - 1) % lx, (y + 1) % ly),
                site_index((x - 1) % lx, y),
                site_index(x, (y - 1) % ly),
                site_index((x + 1) % lx, (y - 1) % ly),
            )
            bonds.extend((target, source) for target in neighbors)
    return tuple(bonds)


def fixed_number_basis(
    lq: int, ne_b: int, ne_c: int
) -> tuple[tuple[tuple[int, ...], tuple[int, ...]], ...]:
    return tuple(
        (b_config, c_config)
        for b_config in occupation_configs(lq, ne_b)
        for c_config in occupation_configs(lq, ne_c)
    )


def add_hopping_terms(
    matrix: np.ndarray,
    states: tuple[tuple[tuple[int, ...], tuple[int, ...]], ...],
    state_index: dict[tuple[tuple[int, ...], tuple[int, ...]], int],
    bonds: tuple[tuple[int, int], ...],
    flavor: str,
    t: float,
) -> None:
    for column, (b_config, c_config) in enumerate(states):
        for target, source in bonds:
            if flavor == "b":
                source_config = list(b_config)
                partner_config = c_config
            else:
                source_config = list(c_config)
                partner_config = b_config

            if source_config[source] == 0:
                continue

            factor = math.sqrt(
                (source_config[target] + 1) * source_config[source]
            )
            source_config[target] += 1
            source_config[source] -= 1
            moved = tuple(source_config)
            if flavor == "b":
                row = state_index[(moved, partner_config)]
            else:
                row = state_index[(partner_config, moved)]
            matrix[row, column] += t * factor


def build_block(
    lx: int, ly: int, t: float, u1: float, u2: float, ne_b: int, ne_c: int
) -> tuple[np.ndarray, np.ndarray, dict[str, np.ndarray]]:
    lq = lx * ly
    states = fixed_number_basis(lq, ne_b, ne_c)
    state_index = {state: index for index, state in enumerate(states)}
    dimension = len(states)

    hamiltonian = np.zeros((dimension, dimension), dtype=float)
    kinetic = np.zeros((dimension, dimension), dtype=float)
    double_occ = np.zeros(dimension, dtype=float)
    square_occ = np.zeros(dimension, dtype=float)

    for index, (b_config, c_config) in enumerate(states):
        b_array = np.asarray(b_config, dtype=float)
        c_array = np.asarray(c_config, dtype=float)
        hamiltonian[index, index] = (u1 + u2) * np.sum(
            b_array * b_array + c_array * c_array
        ) + 2.0 * (u1 - u2) * np.sum(b_array * c_array)
        double_occ[index] = np.sum(b_array * c_array) / float(lq)
        square_occ[index] = (
            0.5
            * np.sum(b_array * (b_array - 1.0) + c_array * (c_array - 1.0))
            / float(lq)
        )

    bonds = directed_triangular_bonds(lx, ly)
    add_hopping_terms(hamiltonian, states, state_index, bonds, "b", t)
    add_hopping_terms(hamiltonian, states, state_index, bonds, "c", t)
    add_hopping_terms(kinetic, states, state_index, bonds, "b", t)
    add_hopping_terms(kinetic, states, state_index, bonds, "c", t)

    diagonal_observables = {
        "doubleOcc": double_occ,
        "squareOcc": square_occ,
    }
    return hamiltonian, kinetic, diagonal_observables


def block_weighted_observables(
    lx: int,
    ly: int,
    t: float,
    u1: float,
    u2: float,
    beta: float,
    mu: float,
    ne_b: int,
    ne_c: int,
) -> tuple[float, dict[str, float]]:
    ne = ne_b + ne_c
    hamiltonian, kinetic, diagonal_observables = build_block(
        lx, ly, t, u1, u2, ne_b, ne_c
    )
    energies, vectors = np.linalg.eigh(hamiltonian)
    weights = np.exp(-(energies - mu * ne) * beta)
    partition = float(np.sum(weights))

    kinetic_vectors = kinetic @ vectors
    kinetic_diagonal = np.einsum("ij,ij->j", vectors, kinetic_vectors)
    probabilities = vectors * vectors

    weighted = {
        "total_NE": partition * float(ne),
        "total_kinetic": float(np.sum(weights * kinetic_diagonal)),
        "doubleOcc": 0.0,
        "squareOcc": 0.0,
        "numsquare_up": partition * float(ne_b * ne_b),
        "numsquare_do": partition * float(ne_c * ne_c),
    }
    for name, diagonal in diagonal_observables.items():
        eigen_diagonal = probabilities.T @ diagonal
        weighted[name] = float(np.sum(weights * eigen_diagonal))
    return partition, weighted


def add_observable_sums(
    target: dict[str, float], source: dict[str, float]
) -> None:
    for name in SCALAR_OBSERVABLES:
        target[name] += source[name]


def normalized_observables(weighted: dict[str, float], partition: float) -> dict[str, float]:
    return {name: weighted[name] / partition for name in SCALAR_OBSERVABLES}


def grand_canonical_observables(params: dict[str, float]) -> dict[str, object]:
    lx = int(params["Lx"])
    ly = int(params["Ly"])
    t = float(params["t"])
    u1 = float(params["U1"])
    u2 = float(params["U2"])
    beta = float(params["beta"])
    mu = float(params["mu"])

    partition_total = 0.0
    weighted_total = {name: 0.0 for name in SCALAR_OBSERVABLES}
    last_shell = {name: 0.0 for name in SCALAR_OBSERVABLES}
    last_ne = 0

    for ne in range(lx * ly * lx * ly + 1):
        shell_partition = 0.0
        shell_weighted = {name: 0.0 for name in SCALAR_OBSERVABLES}
        logging.info("NE = %s", ne)

        for ne_b in range(ne + 1):
            ne_c = ne - ne_b
            logging.info("  NE_b = %s, NE_c = %s", ne_b, ne_c)
            partition, weighted = block_weighted_observables(
                lx, ly, t, u1, u2, beta, mu, ne_b, ne_c
            )
            shell_partition += partition
            add_observable_sums(shell_weighted, weighted)

        partition_total += shell_partition
        add_observable_sums(weighted_total, shell_weighted)
        last_ne = ne
        last_shell = {
            name: shell_weighted[name] / partition_total
            for name in SCALAR_OBSERVABLES
        }

        current = normalized_observables(weighted_total, partition_total)
        logging.info("  partial observables = %s", current)

        if ne > 0 and abs(shell_weighted["total_kinetic"]) < abs(
            0.001 * weighted_total["total_kinetic"]
        ):
            break

    return {
        "parameters": {
            "Lx": lx,
            "Ly": ly,
            "t": t,
            "U1": u1,
            "U2": u2,
            "beta": beta,
            "mu": mu,
        },
        "last_included_NE": last_ne,
        "observables": normalized_observables(weighted_total, partition_total),
        "last_shell_contribution": last_shell,
    }


def write_results(result: dict[str, object]) -> None:
    observables = result["observables"]
    last_shell = result["last_shell_contribution"]
    assert isinstance(observables, dict)
    assert isinstance(last_shell, dict)

    with open("results.txt", "w", encoding="utf-8") as result_file:
        result_file.write(f"{observables['total_NE']}\n")
        result_file.write(f"{observables['total_kinetic']}\n")
        result_file.write(f"{last_shell['total_NE']}\n")
        result_file.write(f"{last_shell['total_kinetic']}\n")
        for name in ("doubleOcc", "squareOcc", "numsquare_up", "numsquare_do"):
            result_file.write(f"{observables[name]}\n")

    Path("results.json").write_text(json.dumps(result, indent=2) + "\n", encoding="utf-8")


def main() -> int:
    logging.info("######################### EDtriangle dense ED begin")
    start_time = time.time()
    params = read_params()
    result = grand_canonical_observables(params)
    write_results(result)
    logging.info("Final results: %s", result["observables"])
    logging.info("Execution completed in %.2f seconds.", time.time() - start_time)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
