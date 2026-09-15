"""Independent finite-Fock checks of both paired exact-diagonalization solvers.

Run with the QuSpin interpreter. The reference uses NumPy ladder matrices and
a density-matrix trace; it imports no production basis, geometry, operator, or
thermal-average helper. All comparisons concern the same finite Hilbert space.
"""

from itertools import product
import unittest

import numpy as np

from src.pairing.benchmarks.ed import ed_pairing_triangle_general as general
from src.pairing.benchmarks.ed import ed_pairing_triangle_tensor as tensor


PARAMETERS = {
    "case": "independent_two_site_fock_reference",
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
    "basis_cap": 100,
}


def finite_fock_operators(parameters):
    """Assemble H-mu*N and observables on |n_b0,n_b1,n_c0,n_c1>."""
    local_dimension = 3
    occupations = np.array(list(product(range(local_dimension), repeat=4)))
    identity = np.eye(local_dimension)
    local_annihilator = np.diag(np.sqrt([1.0, 2.0]), 1)
    annihilators = []
    for mode in range(4):
        operator = np.ones((1, 1))
        for position in range(4):
            operator = np.kron(
                operator, local_annihilator if position == mode else identity
            )
        annihilators.append(operator)
    numbers = [np.diag(occupations[:, mode].astype(float)) for mode in range(4)]
    total_number = sum(numbers)

    # For a 2 x 1 periodic triangular cell, the three positive bond directions
    # are (0->1, 0->0, 0->1) and (1->0, 1->1, 1->0). Including their Hermitian
    # conjugates gives hopping t * [[2,4],[4,2]] for each flavor. The repeated
    # bonds and self loops are physical images of this small periodic cell.
    kinetic = 2.0 * parameters["t"] * total_number
    for left, right in ((0, 1), (2, 3)):
        hopping = annihilators[left].T @ annihilators[right]
        kinetic = kinetic + 4.0 * parameters["t"] * (hopping + hopping.T)

    interaction = np.zeros_like(total_number)
    pair_annihilator = np.zeros_like(total_number)
    double_occupation = np.zeros_like(total_number)
    for b, c in ((0, 2), (1, 3)):
        total_local = numbers[b] + numbers[c]
        relative_local = numbers[b] - numbers[c]
        interaction += (
            parameters["U1"] * (total_local @ total_local)
            + parameters["U2"] * (relative_local @ relative_local)
        )
        pair_annihilator += annihilators[b] @ annihilators[c]
        double_occupation += numbers[b] @ numbers[c] / 2.0

    pair_source = pair_annihilator + pair_annihilator.T
    physical_energy = kinetic + interaction + parameters["Delta"] * pair_source
    grand_hamiltonian = physical_energy - parameters["mu"] * total_number
    operators = {
        "total_NE": total_number,
        "density_total": total_number / 2.0,
        "kinetic_total": kinetic,
        "interaction_energy_total": interaction,
        "energy_total": physical_energy,
        "energy_density": physical_energy / 2.0,
        "grand_energy_total": grand_hamiltonian,
        "doubleOcc": double_occupation,
        "pair_equal": pair_source / 2.0,
        "S_PSF_Gamma": pair_annihilator.T @ pair_annihilator / 4.0,
        "numsquare_up": (numbers[0] + numbers[1]) @ (numbers[0] + numbers[1]),
        "onsite_n2_up": numbers[0] @ numbers[0] + numbers[1] @ numbers[1],
    }
    return grand_hamiltonian, operators, occupations


def direct_thermal_trace(parameters):
    hamiltonian, operators, occupations = finite_fock_operators(parameters)
    energies, eigenvectors = np.linalg.eigh(hamiltonian)
    # These 81-state test Hamiltonians have modest energies, so the direct
    # Boltzmann factors are representable without production energy shifting.
    weights = np.exp(-parameters["beta"] * energies)
    partition = weights.sum()
    density_matrix = (eigenvectors * (weights / partition)) @ eigenvectors.T
    values = {
        name: float(np.trace(density_matrix @ operator))
        for name, operator in operators.items()
    }
    return {
        "hamiltonian": hamiltonian,
        "occupations": occupations,
        "energies": energies,
        "free_energy": float(-np.log(partition) / parameters["beta"]),
        "observables": values,
    }


class IndependentEDPhysicsTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.references = {}
        cls.results = {}
        cls.spectra = {}
        for sign in (1, -1):
            parameters = {**PARAMETERS, "Delta": sign * PARAMETERS["Delta"]}
            cls.references[sign] = direct_thermal_trace(parameters)
            cls.results[sign] = {
                "general": general.run_ed(parameters),
                "tensor": tensor.run_tensor_ed(parameters),
            }
            general_basis = general.build_basis(2, 1, 2, None)
            tensor_basis = tensor.build_tensor_basis(2, 1, 2)
            cls.spectra[sign] = {
                "general": np.linalg.eigvalsh(
                    general.build_hamiltonian(parameters, general_basis).toarray()
                ),
                "tensor": np.linalg.eigvalsh(
                    tensor.build_tensor_hamiltonian(parameters, tensor_basis).toarray()
                ),
            }

    def test_finite_pairing_spectra_match_independent_fock_hamiltonian(self):
        for sign, spectra in self.spectra.items():
            for implementation, energies in spectra.items():
                with self.subTest(sign=sign, implementation=implementation):
                    self.assertEqual(len(energies), 81)
                    np.testing.assert_allclose(
                        energies, self.references[sign]["energies"],
                        rtol=1.0e-12, atol=1.0e-12,
                    )

    def test_thermal_observables_and_free_energy_match_independent_trace(self):
        for sign, results in self.results.items():
            reference = self.references[sign]
            for implementation, result in results.items():
                with self.subTest(sign=sign, implementation=implementation):
                    self.assertAlmostEqual(
                        result["free_energy"], reference["free_energy"], delta=1.0e-11
                    )
                    for name, expected in reference["observables"].items():
                        with self.subTest(observable=name):
                            self.assertAlmostEqual(
                                result["observables"][name]["value"], expected,
                                delta=1.0e-11,
                            )

    def test_pairing_sign_is_a_flavor_phase_rotation(self):
        positive, negative = self.references[1], self.references[-1]
        phase = (-1.0) ** positive["occupations"][:, 2:].sum(axis=1)
        rotated = phase[:, None] * positive["hamiltonian"] * phase[None, :]
        np.testing.assert_allclose(rotated, negative["hamiltonian"], atol=1.0e-13)
        for implementation in ("general", "tensor"):
            with self.subTest(implementation=implementation):
                plus = self.results[1][implementation]
                minus = self.results[-1][implementation]
                self.assertAlmostEqual(plus["free_energy"], minus["free_energy"], delta=1.0e-12)
                for name in positive["observables"]:
                    sign = -1 if name == "pair_equal" else 1
                    self.assertAlmostEqual(
                        plus["observables"][name]["value"],
                        sign * minus["observables"][name]["value"],
                        delta=1.0e-11,
                    )
                self.assertLess(plus["observables"]["pair_equal"]["value"], -1.0e-3)

    def test_density_and_grand_energy_match_free_energy_derivatives(self):
        step = 1.0e-4
        mu = PARAMETERS["mu"]
        f_minus_mu = direct_thermal_trace({**PARAMETERS, "mu": mu - step})["free_energy"]
        f_plus_mu = direct_thermal_trace({**PARAMETERS, "mu": mu + step})["free_energy"]
        number_from_derivative = -(f_plus_mu - f_minus_mu) / (2.0 * step)
        beta = PARAMETERS["beta"]
        f_minus_beta = direct_thermal_trace({**PARAMETERS, "beta": beta - step})["free_energy"]
        f_plus_beta = direct_thermal_trace({**PARAMETERS, "beta": beta + step})["free_energy"]
        grand_energy_from_derivative = (
            (beta + step) * f_plus_beta - (beta - step) * f_minus_beta
        ) / (2.0 * step)
        for implementation, result in self.results[1].items():
            with self.subTest(implementation=implementation):
                self.assertAlmostEqual(
                    result["observables"]["total_NE"]["value"],
                    number_from_derivative, delta=2.0e-8,
                )
                self.assertAlmostEqual(
                    result["observables"]["grand_energy_total"]["value"],
                    grand_energy_from_derivative, delta=2.0e-8,
                )


if __name__ == "__main__":
    unittest.main()
