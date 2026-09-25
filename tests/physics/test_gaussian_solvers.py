"""Compare actual Fortran measurements with an independent Gaussian solution.

Run with BAFQMC_RUN_MPI_TESTS=1. No QuSpin installation or stored benchmark data is
needed: these are four-bin, zero-interaction calculations in fresh directories.
The HS fields decouple, so every bin has the same exact expectation value.

For each triangular-lattice momentum, xi=epsilon-mu and
omega=sqrt(xi**2-Delta**2). The two-mode Bogoliubov transformation gives
n_k=(xi/omega*coth(beta*omega/2)-1)/2 and
m_k=<b_k c_-k>=-Delta*coth(beta*omega/2)/(2*omega).
The code uses a positive real pairing coefficient. Thus pair_equal=2*mean(m),
rho=2*mean(n), and E=2*sum(epsilon*n)+Ns*Delta*pair_equal.

At nonzero K, the disconnected density contribution vanishes. Wick's theorem
then gives S_DW(K)=2/Ns**2*sum[n_k*(1+n_(k+K))+m_k*m_(k+K)]. The second term
tests anomalous cross-flavor contractions; the 1 tests the bosonic commutator.
S_SF(K)=2*n_K/Ns. These momentum sums do not call either solver's reference,
geometry, Green-function, or estimator implementation.
"""
from __future__ import annotations

import json
import math
import os
from pathlib import Path
import shutil
import subprocess
import sys
import tempfile
import unittest


ROOT = Path(__file__).resolve().parents[2]
NBIN = 4
BETA = 1.0
MU = -5.0
# This is a floating-point comparison to an exact quadratic solution, not a
# statistical acceptance criterion. It also resolves the small anomalous terms.
ABS_TOL = 2.0e-9
REL_TOL = 2.0e-8
PAPER_OBSERVABLES = ("density_total", "energy_density", "sf_K", "dw_K")


def gaussian_observables(delta: float) -> dict[str, float]:
    """Thermal observables on a periodic 3x3 lattice with positive t=1."""
    length = 3
    nsite = length**2
    occupation, anomalous, dispersion = {}, {}, {}
    for x in range(length):
        for y in range(length):
            # k.a1=2*pi*x/L; k.a2=2*pi*y/L for the three bond directions
            # a1, a2, a1-a2. K=(4*pi/3,0) has indices (2,1).
            k1, k2 = 2 * math.pi * x / length, 2 * math.pi * y / length
            epsilon = 2 * (math.cos(k1) + math.cos(k2) + math.cos(k1-k2))
            xi = epsilon - MU
            if xi <= abs(delta):
                raise ValueError("Gaussian thermal trace requires xi > abs(Delta)")
            omega = math.sqrt(xi*xi-delta*delta)
            coth = 1 / math.tanh(BETA*omega/2)
            occupation[x, y] = (xi/omega*coth-1)/2
            anomalous[x, y] = -delta*coth/(2*omega)
            dispersion[x, y] = epsilon
    pair = 2 * math.fsum(anomalous.values()) / nsite
    energy = 2 * math.fsum(dispersion[k]*occupation[k] for k in occupation)
    energy += nsite*delta*pair
    density_sum = math.fsum(
        occupation[x, y] * (1 + occupation[(x+2) % length, (y+1) % length])
        + anomalous[x, y] * anomalous[(x+2) % length, (y+1) % length]
        for x, y in occupation
    )
    return {
        "density_total": 2 * math.fsum(occupation.values()) / nsite,
        "energy_density": energy / nsite,
        "sf_K": 2 * occupation[2, 1] / nsite,
        "dw_K": 2 * density_sum / nsite**2,
        "pair_equal": pair,
    }


@unittest.skipUnless(
    os.environ.get("BAFQMC_RUN_MPI_TESTS") == "1",
    "set BAFQMC_RUN_MPI_TESTS=1 to build and run the Fortran Gaussian checks",
)
class GaussianSolverTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        if sys.platform != "linux":
            raise RuntimeError("Gaussian solver checks require Linux/WSL")
        if shutil.which("mpirun") is None:
            raise RuntimeError("Gaussian solver checks require an MPI launcher")
        temporary = tempfile.TemporaryDirectory(prefix="bafqmc-gaussian-")
        cls.addClassCleanup(temporary.cleanup)
        cls.work = Path(temporary.name)
        cls.env = os.environ.copy()
        cls.env.update({key: "1" for key in (
            "OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS", "NUMBA_NUM_THREADS",
        )})
        cls.results = {}
        for model in ("number_conserving", "pairing"):
            cls.run_command(["make", "-C", str(ROOT / "src" / model), "build"])
        for nwrap in (10, 2):
            cls.results["number_conserving", 0.0, nwrap] = cls.run_number(nwrap)
        for delta, nwrap in ((0.0, 10), (0.2, 10), (-0.2, 10), (0.2, 2)):
            cls.results["pairing", delta, nwrap] = cls.run_pairing(delta, nwrap)

    @classmethod
    def run_command(cls, command):
        result = subprocess.run(command, cwd=ROOT, env=cls.env, text=True,
                                stdout=subprocess.PIPE, stderr=subprocess.STDOUT,
                                timeout=180)
        if result.returncode:
            logs = "\n".join(
                f"{log.relative_to(cls.work)}:\n{log.read_text(errors='replace')[-2000:]}"
                for log in sorted(cls.work.rglob("*.log"))
            )
            raise RuntimeError(f"Command failed: {command}\n{result.stdout[-6000:]}\n{logs}")

    @classmethod
    def run_number(cls, nwrap):
        work = cls.work / f"number-wrap{nwrap}"
        inputs = work / "inputs"
        source = ROOT / "examples/number_conserving/examples/triangle_3x3_free"
        inputs.mkdir(parents=True)
        for name in ("paramC_sets.txt", "confin.txt", "seeds.txt", "params.json"):
            shutil.copy2(source / name, inputs / name)
        rows = (inputs / "paramC_sets.txt").read_text().splitlines()
        rows[0] = f"0.0 0.0 {MU}"
        rows[1] = f"3 3 100 {BETA}"
        rows[2] = "3 3 100"
        rows[3] = f"{nwrap} {NBIN} 1 1.0"
        rows[5] = ".false. 0 1.0 1.0"
        (inputs / "paramC_sets.txt").write_text("\n".join(rows) + "\n")
        manifest = {"cases": [{"id": "gaussian", "input_dir": "inputs"}]}
        path = work / "manifest.json"
        path.write_text(json.dumps(manifest))
        cls.run_command([sys.executable, str(ROOT / "src/number_conserving/run_paper.py"),
                         "--manifest", str(path), "--output", str(work / "output"),
                         "--mode", "dqmc"])
        return cls.read_measurements(work / "output/gaussian/dqmc", PAPER_OBSERVABLES)

    @classmethod
    def run_pairing(cls, delta, nwrap):
        work = cls.work / f"pairing-delta{delta}-wrap{nwrap}"
        work.mkdir()
        source = ROOT / "src/pairing/benchmarks/campaigns/triangle_pairing_pipeline_smoke.json"
        manifest = json.loads(source.read_text())
        manifest["parameters"].update(U1=0.0, U2=0.0, mu=MU, beta=BETA)
        manifest["deltas"] = [delta]
        manifest["dqmc_defaults"].update(Nbin=NBIN, Nwrap=nwrap, dtau=0.01,
                                          is_warm=False, Nwarm=0)
        path = work / "manifest.json"
        path.write_text(json.dumps(manifest))
        cls.run_command([sys.executable, str(ROOT / "src/pairing/run_paper.py"),
                         "--manifest", str(path), "--output", str(work / "output"),
                         "--mode", "dqmc", "--seed", "24681357"])
        outputs = list((work / "output/inputs").iterdir())
        if len(outputs) != 1:
            raise AssertionError(f"expected one Gaussian case, found {outputs}")
        return cls.read_measurements(outputs[0], (*PAPER_OBSERVABLES, "pair_equal"))

    @staticmethod
    def read_measurements(directory, observables):
        values = {}
        for observable in observables:
            rows = [[float(token.replace("D", "E").replace("d", "e"))
                     for token in line.split()]
                    for line in (directory / observable).read_text().splitlines() if line.strip()]
            ncolumns = 2 if observable in ("sf_K", "dw_K") else 1
            if len(rows) != NBIN or any(len(row) != ncolumns for row in rows):
                raise AssertionError(f"{directory / observable}: expected {NBIN}x{ncolumns} measurements")
            if any(not math.isfinite(value) for row in rows for value in row):
                raise AssertionError(f"{directory / observable}: nonfinite measurement")
            values[observable] = rows
        return values

    def assert_precise(self, actual, expected, label):
        self.assertTrue(math.isclose(actual, expected, abs_tol=ABS_TOL, rel_tol=REL_TOL),
                        f"{label}: measured={actual:.16g}, exact={expected:.16g}")

    def test_all_bins_match_independent_bogoliubov_solution(self):
        for (model, delta, nwrap), measurements in self.results.items():
            exact = gaussian_observables(delta)
            for observable, rows in measurements.items():
                for index, row in enumerate(rows):
                    with self.subTest(model=model, delta=delta, nwrap=nwrap,
                                      observable=observable, bin=index):
                        self.assert_precise(row[0], exact[observable], observable)
                        if len(row) == 2:
                            self.assert_precise(row[1], 0.0, f"Im {observable}")

    def test_zero_pairing_matches_number_conserving_solver(self):
        for observable in PAPER_OBSERVABLES:
            self.assert_precise(self.results["pairing", 0.0, 10][observable][0][0],
                                self.results["number_conserving", 0.0, 10][observable][0][0],
                                observable)

    def test_pairing_phase_changes_only_the_anomalous_sign(self):
        positive, negative = self.results["pairing", 0.2, 10], self.results["pairing", -0.2, 10]
        for observable in PAPER_OBSERVABLES:
            self.assert_precise(positive[observable][0][0], negative[observable][0][0], observable)
        self.assert_precise(positive["pair_equal"][0][0], -negative["pair_equal"][0][0], "pair_equal")
        self.assertLess(positive["pair_equal"][0][0], -0.01)

    def test_stabilization_interval_preserves_observables(self):
        for model, delta in (("number_conserving", 0.0), ("pairing", 0.2)):
            for observable in self.results[model, delta, 10]:
                self.assert_precise(self.results[model, delta, 10][observable][0][0],
                                    self.results[model, delta, 2][observable][0][0],
                                    f"{model} {observable}")


if __name__ == "__main__":
    unittest.main()
