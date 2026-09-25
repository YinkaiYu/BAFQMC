"""Check interacting HS propagation against independent dense matrix products.

Run with BAFQMC_RUN_MPI_TESTS=1. The actual number-conserving executable reads
a prescribed nonuniform auxiliary field with U1>0 and U2<0. Zero proposal
displacement and disabled warmup keep that field fixed. Every bin contains a
left and right traversal of imaginary time; three stabilization intervals are
checked against the same direct Gaussian trace at each cyclic time boundary.

This checks configuration weights and fixed-field observables, including both
interaction channels. It does not test Metropolis acceptance decisions or
estimate interacting thermal averages. No stored benchmark measurements are used.
"""

from __future__ import annotations

import os
from pathlib import Path
import shutil
import subprocess
import sys
import tempfile
import unittest

import numpy as np


ROOT = Path(__file__).resolve().parents[2]
SOLVER = ROOT / "src/number_conserving"
LX = LY = 3
NSITE = LX * LY
NSLICE = 6
NBIN = 3
BETA = 0.6
MU = -5.0
U1, U2 = 0.7, -0.1
DT = BETA / NSLICE
TOLERANCE = 2.0e-10


def prescribed_fields():
    """Array order matches the input's time, site, auxiliary-channel order."""
    time = np.arange(1, NSLICE + 1, dtype=float)[:, None]
    site = np.arange(1, NSITE + 1, dtype=float)[None, :]
    density_field = (
        0.23 + 0.31 * np.sin(0.71 * time + 0.37 * site)
        + 0.09 * np.cos(0.53 * time * site)
    )
    relative_field = (
        -0.19 + 0.28 * np.cos(0.43 * time - 0.61 * site)
        + 0.13 * np.sin(0.29 * time * site)
    )
    return np.stack((relative_field, density_field), axis=-1)


def direct_fixed_field_reference(fields):
    """Construct the physical hopping and HS factors without solver helpers."""
    hopping = np.zeros((NSITE, NSITE))
    for y in range(LY):
        for x in range(LX):
            source = y * LX + x
            # Three positive triangular-lattice bonds plus their conjugates.
            for dx, dy in ((1, 0), (0, 1), (-1, 1)):
                destination = ((y + dy) % LY) * LX + (x + dx) % LX
                hopping[source, destination] += 1.0
                hopping[destination, source] += 1.0
    eigenvalues, eigenvectors = np.linalg.eigh(hopping - MU * np.eye(NSITE))
    kinetic_slice = (eigenvectors * np.exp(-DT * eigenvalues)) @ eigenvectors.T
    hs_exponents = (
        1j * np.sqrt(2.0 * U1 * DT) * fields[:, :, 0]
        + np.sqrt(-2.0 * U2 * DT) * fields[:, :, 1]
    )
    slices = [kinetic_slice @ np.diag(np.exp(exponent)) for exponent in hs_exponents]

    # At boundary tau the ordered product is
    # B_tau ... B_1 B_M ... B_(tau+1). Rebuild each product from scratch;
    # no rank-one updates, wrap algorithm, or production Green function is used.
    products, greens = [], []
    for boundary in range(1, NSLICE + 1):
        product = np.eye(NSITE, dtype=complex)
        for index in (*range(boundary, NSLICE), *range(boundary)):
            product = slices[index] @ product
        products.append(product)
        greens.append(np.linalg.inv(np.eye(NSITE) - product))

    site_densities, kinetic_values, double_values, onsite_values = [], [], [], []
    for green in greens:
        occupation = np.diag(green) - 1.0
        site_densities.append(2.0 * occupation.real)
        kinetic_values.append(2.0 * np.trace(hopping @ (green - np.eye(NSITE))).real / NSITE)
        double_values.append(float(np.mean(np.abs(occupation) ** 2)))
        onsite_values.append(float(np.sum(2.0 * occupation**2 + occupation).real))

    density_site = np.mean(site_densities, axis=0)
    onsite_n2 = float(np.mean(onsite_values))
    double_occupation = float(np.mean(double_values))
    kinetic = float(np.mean(kinetic_values))
    interaction = (
        2.0 * (U1 + U2) * onsite_n2 / NSITE
        + 2.0 * (U2 - U1) * double_occupation
    )
    matrix = np.eye(NSITE) - products[-1]
    _, log_abs_determinant = np.linalg.slogdet(matrix)
    pole_eigenvalues = np.linalg.eigvals(matrix)
    pole_distance = float(np.min(np.abs(pole_eigenvalues)))
    return {
        "density_site_total": density_site,
        "density_total": float(density_site.mean()),
        "kinetic": kinetic,
        "doubleOcc": double_occupation,
        "onsite_n2_up": onsite_n2,
        "onsite_n2_do": onsite_n2,
        "interaction_energy_density": interaction,
        "energy_density": kinetic + interaction,
        "log_weight": float(-0.5 * np.sum(fields**2) - 2.0 * log_abs_determinant),
        "pole_distance": pole_distance,
        "pole_x": float(-np.log10(pole_distance)),
        "green_spectral_radius": 1.0 / pole_distance,
        "green_smax": float(np.linalg.svd(greens[-1], compute_uv=False)[0]),
    }


@unittest.skipUnless(
    os.environ.get("BAFQMC_RUN_MPI_TESTS") == "1",
    "set BAFQMC_RUN_MPI_TESTS=1 for the fixed interacting-field solver checks",
)
class FixedFieldSolverTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        if sys.platform != "linux" or shutil.which("mpirun") is None:
            raise RuntimeError("Fixed-field solver checks require Linux/WSL and mpirun")
        temporary = tempfile.TemporaryDirectory(prefix="bafqmc-fixed-field-")
        cls.addClassCleanup(temporary.cleanup)
        cls.work = Path(temporary.name)
        cls.environment = os.environ.copy()
        cls.environment.update({name: "1" for name in (
            "OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS", "NUMBA_NUM_THREADS",
        )})
        cls.fields = prescribed_fields()
        cls.reference = direct_fixed_field_reference(cls.fields)
        cls.run_command(["make", "-C", str(SOLVER), "build"], ROOT)
        cls.outputs = {}
        for nwrap in (1, 2, NSLICE):
            directory = cls.work / f"wrap-{nwrap}"
            directory.mkdir()
            cls.write_inputs(directory, nwrap)
            cls.run_command(["mpirun", "-np", "1", str(SOLVER / "build/bosonDQMC.out")], directory)
            cls.outputs[nwrap] = directory

    @classmethod
    def run_command(cls, command, directory):
        result = subprocess.run(
            command, cwd=directory, env=cls.environment, text=True,
            stdout=subprocess.PIPE, stderr=subprocess.STDOUT, timeout=180,
        )
        if result.returncode:
            raise RuntimeError(f"Command failed: {command}\n{result.stdout[-6000:]}")

    @classmethod
    def write_inputs(cls, directory, nwrap):
        rows = (
            f"{U1} {U2} {MU}",
            f"{LX} {LY} {NSLICE} {BETA}",
            f"{LX} {LY} {NSLICE}",
            f"{nwrap} {NBIN} 1 0.0",
            ".false. 0",
            ".false. 0 0.0 0.0",
            "1 0.0 0.0 0.0",
        )
        (directory / "paramC_sets.txt").write_text("\n".join(rows) + "\n")
        # A nonzero leading seed selects explicit-field input, avoiding random
        # initialization. The second seed file is a required executable input.
        (directory / "confin.txt").write_text(
            "1357911\n" + "".join(f"{value:.17g}\n" for value in cls.fields.ravel())
        )
        (directory / "seeds.txt").write_text("24681357\n")

    def assert_measurements(self, directory, name, expected):
        path = directory / name
        rows = [[float(token.replace("D", "E").replace("d", "e")) for token in line.split()]
                for line in path.read_text().splitlines() if line.strip()]
        expected_row = np.atleast_1d(expected)
        actual = np.array(rows)
        self.assertEqual(actual.shape, (NBIN, expected_row.size), str(path))
        self.assertTrue(np.isfinite(actual).all(), str(path))
        np.testing.assert_allclose(
            actual, np.broadcast_to(expected_row, actual.shape),
            rtol=TOLERANCE, atol=TOLERANCE, err_msg=str(path),
        )

    def test_weight_and_poles_match_direct_determinant(self):
        self.assertGreater(abs(self.reference["log_weight"]), 1.0)
        self.assertGreater(self.reference["green_spectral_radius"], 1.05)
        for nwrap, directory in self.outputs.items():
            for name in ("log_weight", "pole_distance", "pole_x", "green_spectral_radius", "green_smax"):
                with self.subTest(nwrap=nwrap, observable=name):
                    self.assert_measurements(directory, name, self.reference[name])

    def test_cyclic_time_averaged_observables_match_dense_gaussian_traces(self):
        self.assertGreater(self.reference["density_total"], 0.05)
        self.assertGreater(np.ptp(self.reference["density_site_total"]), 1.0e-3)
        self.assertGreater(abs(self.reference["interaction_energy_density"]), 1.0e-3)
        for nwrap, directory in self.outputs.items():
            for name in (
                "density_total", "density_site_total", "kinetic", "doubleOcc",
                "onsite_n2_up", "onsite_n2_do", "interaction_energy_density", "energy_density",
            ):
                with self.subTest(nwrap=nwrap, observable=name):
                    self.assert_measurements(directory, name, self.reference[name])

    def test_auxiliary_configuration_remains_fixed_across_both_sweep_directions(self):
        for nwrap, directory in self.outputs.items():
            with self.subTest(nwrap=nwrap):
                rows = (directory / "confout.txt").read_text().splitlines()
                output_fields = np.array([float(row.replace("D", "E")) for row in rows[1:] if row.strip()])
                np.testing.assert_allclose(output_fields, self.fields.ravel(), rtol=0.0, atol=2.0e-15)


if __name__ == "__main__":
    unittest.main()
