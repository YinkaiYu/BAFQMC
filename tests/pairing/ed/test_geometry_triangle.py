#!/usr/bin/env python3
from __future__ import annotations

import sys
import unittest
from pathlib import Path

import numpy as np

REPO_ROOT = Path(__file__).resolve().parents[3] / "src/pairing"
sys.path.insert(0, str(REPO_ROOT / "benchmarks/ed"))

import ed_pairing_triangle_general as general
from geometry_triangle import A1, A2, B1, B2


class TriangleGeometryTests(unittest.TestCase):

    def test_reciprocal_basis_matches_real_space_basis(self) -> None:
        avecs = (np.array(A1), np.array(A2))
        bvecs = (np.array(B1), np.array(B2))
        target = 2.0 * np.pi * np.eye(2)
        actual = np.array([[float(np.dot(b, a)) for a in avecs] for b in bvecs])

        self.assertTrue(np.allclose(actual, target, rtol=1.0e-12, atol=1.0e-12), actual)

    def test_ed_phase_indexing_matches_no_pairing_modulo_convention(self) -> None:
        self.assertEqual(general._npbc(0, 3), 1)
        self.assertEqual(general._npbc(-1, 3), 3)
        self.assertEqual(general._npbc(3, 3), 1)
        self.assertEqual(general._npbc(4, 3), 2)
        self.assertTrue(np.allclose(general._aimj_vector(0, 3), (0.0, 0.0)))

    def test_k_index_fails_closed_for_non_compatible_geometry(self) -> None:
        self.assertEqual(general.structure_k_index(3, 3), 5)
        with self.assertRaises(ValueError):
            general.structure_k_index(2, 2)
        with self.assertRaises(ValueError):
            general.structure_k_index(3, 2)

    def test_k_phase_matrix_matches_direct_fourier_kernel(self) -> None:
        lx = 3
        ly = 3
        k_index = general.structure_k_index(lx, ly)
        pairing_phases = general.triangular_phase_matrix(lx, ly, k_index)
        expected = np.empty_like(pairing_phases)
        kx, ky = general._k_vector(k_index, lx, ly)
        for i in range(lx * ly):
            ix, iy = i % lx + 1, i // lx + 1
            for j in range(lx * ly):
                jx, jy = j % lx + 1, j // lx + 1
                dx = (ix - jx) % lx
                dy = (iy - jy) % ly
                rx = dx * A1[0] + dy * A2[0]
                ry = dx * A1[1] + dy * A2[1]
                expected[i, j] = np.exp(1j * (kx * rx + ky * ry))
        self.assertTrue(
            np.allclose(pairing_phases, expected, rtol=1.0e-12, atol=1.0e-12),
            pairing_phases - expected,
        )

    def test_k_phase_matrix_is_hermitian_kernel(self) -> None:
        lx = 3
        ly = 3
        phases = general.triangular_phase_matrix(
            lx,
            ly,
            general.structure_k_index(lx, ly),
        )

        self.assertTrue(
            np.allclose(np.diag(phases), 1.0, rtol=1.0e-12, atol=1.0e-12),
            np.diag(phases),
        )
        self.assertTrue(
            np.allclose(phases, phases.conj().T, rtol=1.0e-12, atol=1.0e-12),
            phases - phases.conj().T,
        )

    def test_k_structure_factor_operators_are_hermitian_positive(self) -> None:
        params = {"Lx": 3, "Ly": 3}
        basis = general.build_basis(3, 3, 1, 2)

        for operator in (
            general.sf_structure_operator(params, basis),
            general.dw_structure_operator(params, basis),
        ):
            dense = operator.toarray()
            hermitian_part = 0.5 * (dense + dense.conj().T)
            self.assertTrue(
                np.allclose(dense, dense.conj().T, rtol=1.0e-12, atol=1.0e-12)
            )
            self.assertGreaterEqual(float(np.linalg.eigvalsh(hermitian_part).min()), -1.0e-12)


if __name__ == "__main__":
    unittest.main()
