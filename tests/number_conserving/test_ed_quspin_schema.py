#!/usr/bin/env python3
"""Import-light schema tests for the 3x3 QuSpin ED reference helper."""

from __future__ import annotations

import importlib.util
import json
import math
import os
import sys
import subprocess
import tempfile
import unittest
from pathlib import Path


REPO_ROOT = Path(__file__).resolve().parents[2] / "src/number_conserving"
ED_SCRIPT = REPO_ROOT / "benchmarks" / "ed" / "EDtriangle_quspin_3x3.py"
QUSPIN_PYTHON = Path(os.environ.get("ED_PYTHON", sys.executable))


def load_ed_module():
    spec = importlib.util.spec_from_file_location("edtriangle_quspin_3x3", ED_SCRIPT)
    if spec is None or spec.loader is None:
        raise RuntimeError(f"Cannot import {ED_SCRIPT}")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


class EDQuSpinSchemaTests(unittest.TestCase):
    def setUp(self) -> None:
        self.ed = load_ed_module()

    def test_empty_weighted_observables_matches_schema(self) -> None:
        observables = self.ed.empty_weighted_observables()

        self.assertEqual(set(observables), set(self.ed.ED_OBSERVABLES))
        self.assertTrue(all(value == 0.0 for value in observables.values()))

    def test_schema_includes_total_flavor_number_second_moments(self) -> None:
        self.assertIn("numsquare_up", self.ed.ED_OBSERVABLES)
        self.assertIn("numsquare_do", self.ed.ED_OBSERVABLES)

    def test_normalize_weighted_observables_divides_partition(self) -> None:
        weighted = {
            name: float(index + 1)
            for index, name in enumerate(self.ed.ED_OBSERVABLES)
        }

        observables = self.ed.normalize_weighted_observables(weighted, 2.0)

        self.assertEqual(set(observables), set(self.ed.ED_OBSERVABLES))
        for name in self.ed.ED_OBSERVABLES:
            self.assertEqual(observables[name], weighted[name] / 2.0)

    def test_classify_status_converged_for_bounded_repulsive_case(self) -> None:
        relative_tail = {name: 0.01 for name in self.ed.ED_OBSERVABLES}
        policy = {"tail_tolerance": 0.05}

        status = self.ed.classify_status(
            {"U1": 0.0}, relative_tail, policy, completed_shells=1
        )

        self.assertEqual(status, "converged")

    def test_classify_status_waits_for_minimum_shells_before_converging(self) -> None:
        relative_tail = {name: 0.0 for name in self.ed.ED_OBSERVABLES}
        policy = {"tail_tolerance": 0.05, "min_completed_shells": 2}

        status = self.ed.classify_status(
            {"U1": 0.0}, relative_tail, policy, completed_shells=0
        )

        self.assertEqual(status, "incomplete")

    def test_classify_status_accepts_negative_u1_low_density_cutoff(self) -> None:
        relative_tail = {name: 1.0 for name in self.ed.ED_OBSERVABLES}
        policy = {
            "tail_tolerance": 0.05,
            "max_density_total": 0.2,
            "min_completed_shells": 3,
        }
        observables = self.ed.empty_weighted_observables()
        observables["density_total"] = 0.1

        status = self.ed.classify_status(
            {"U1": -0.1},
            relative_tail,
            policy,
            completed_shells=3,
            observables=observables,
        )

        self.assertEqual(status, "low_density_cutoff_accepted")
        self.assertNotEqual(status, "converged")

    def test_classify_status_incomplete_when_policy_not_met(self) -> None:
        relative_tail = {name: 0.10 for name in self.ed.ED_OBSERVABLES}
        policy = {
            "tail_tolerance": 0.05,
            "max_density_total": 0.2,
            "min_completed_shells": 3,
        }

        positive_status = self.ed.classify_status(
            {"U1": 0.0}, relative_tail, policy, completed_shells=4
        )
        negative_status = self.ed.classify_status(
            {"U1": -0.1},
            relative_tail,
            policy,
            completed_shells=2,
            observables={"density_total": 0.1},
        )

        self.assertEqual(positive_status, "incomplete")
        self.assertEqual(negative_status, "incomplete")

    def test_dry_run_payload_is_explicitly_untrusted(self) -> None:
        payload = self.ed.dry_run_payload(
            {
                "Lx": 3,
                "Ly": 3,
                "beta": 1.0,
                "mu": -5.0,
                "U1": 0.0,
                "U2": 1.0,
                "max_total_particles": 2,
            },
            {"tail_tolerance": 0.05},
        )

        self.assertTrue(payload["parameters"]["dry_run"])
        self.assertEqual(payload["status"], "dry_run")

    def test_load_inputs_accepts_manifest_style_max_total_particles(self) -> None:
        args = self.ed._parse_args(
            [
                "--Lx",
                "3",
                "--Ly",
                "3",
                "--beta",
                "1.0",
                "--mu",
                "-5.0",
                "--U1",
                "0.0",
                "--U2",
                "1.0",
                "--max-total-particles",
                "4",
                "--tail-tolerance",
                "0.001",
            ]
        )

        parameters, policy = self.ed._load_inputs(args)

        self.assertEqual(parameters["max_total_particles"], 4)
        self.assertEqual(policy["tail_tolerance"], 0.001)

    def test_load_inputs_migrates_json_policy_particle_cutoff(self) -> None:
        with tempfile.TemporaryDirectory() as tmpdir:
            params = Path(tmpdir) / "params.json"
            params.write_text(
                json.dumps(
                    {
                        "parameters": {
                            "Lx": 3,
                            "Ly": 3,
                            "beta": 1.0,
                            "mu": -5.0,
                            "U1": 0.0,
                            "U2": 1.0,
                        },
                        "convergence_policy": {
                            "tail_tolerance": 1.0e-3,
                            "max_total_particles": 5,
                        },
                    }
                ),
                encoding="utf-8",
            )
            args = self.ed._parse_args([str(params)])

            parameters, policy = self.ed._load_inputs(args)

        self.assertEqual(parameters["max_total_particles"], 5)
        self.assertNotIn("max_total_particles", policy)

    def test_k_phase_matches_three_by_three_fortran_k_point(self) -> None:
        self.assertAlmostEqual(self.ed._k_phase(3, 3, 0, 0).real, 1.0)
        self.assertAlmostEqual(self.ed._k_phase(3, 3, 0, 0).imag, 0.0)
        self.assertAlmostEqual(self.ed._k_phase(3, 3, 1, 0).real, -0.5)
        self.assertAlmostEqual(
            self.ed._k_phase(3, 3, 1, 0).imag,
            -math.sqrt(3.0) / 2.0,
        )
        self.assertAlmostEqual(self.ed._k_phase(3, 3, 0, 1).real, -0.5)
        self.assertAlmostEqual(
            self.ed._k_phase(3, 3, 0, 1).imag,
            math.sqrt(3.0) / 2.0,
        )

    def test_result_payload_schema_and_json_serializable(self) -> None:
        observables = {
            name: float(index)
            for index, name in enumerate(self.ed.ED_OBSERVABLES)
        }
        last_shell = {
            name: float(index) / 10.0
            for index, name in enumerate(self.ed.ED_OBSERVABLES)
        }
        relative_tail = self.ed.relative_tail(last_shell, observables)
        policy = {"tail_tolerance": 0.05}

        payload = self.ed.result_payload(
            parameters={"Lx": 3, "Ly": 3, "U1": 0.0},
            observables=observables,
            last_shell_contribution=last_shell,
            relative_last_shell=relative_tail,
            convergence_policy=policy,
            status="incomplete",
            completed_shells=4,
        )

        self.assertEqual(
            set(payload),
            {
                "schema_version",
                "parameters",
                "observables",
                "last_shell_contribution",
                "relative_last_shell",
                "convergence_policy",
                "status",
                "completed_shells",
                "timestamp_utc",
            },
        )
        self.assertEqual(set(payload["observables"]), set(self.ed.ED_OBSERVABLES))
        json.dumps(payload)

    def test_runtime_stops_after_reliable_shell_when_policy_allows(self) -> None:
        calls = []

        def fake_grand_shell_pair(runtime, **kwargs):
            total_ne = int(kwargs["ne_b"]) + int(kwargs["ne_c"])
            calls.append(total_ne)
            if total_ne == 1:
                weighted = {name: 1.0 for name in self.ed.ED_OBSERVABLES}
                return weighted, 1.0
            weighted = {name: 1.0e-8 for name in self.ed.ED_OBSERVABLES}
            return weighted, 1.0e-8

        original = self.ed._grand_shell_pair
        self.ed._grand_shell_pair = fake_grand_shell_pair
        try:
            with tempfile.TemporaryDirectory() as tmpdir:
                output = Path(tmpdir) / "results.json"
                payload = self.ed._run_ed_reference_runtime(
                    {},
                    parameters={
                        "Lx": 3,
                        "Ly": 3,
                        "beta": 1.0,
                        "mu": -5.0,
                        "U1": 0.0,
                        "U2": 1.0,
                        "max_total_particles": 5,
                    },
                    policy={
                        "tail_tolerance": 1.0e-3,
                        "min_completed_shells": 2,
                        "stop_when_reliable": True,
                    },
                    output_path=output,
                )
                written = json.loads(output.read_text(encoding="utf-8"))
        finally:
            self.ed._grand_shell_pair = original

        self.assertEqual(payload["status"], "converged")
        self.assertEqual(payload["completed_shells"], 2)
        self.assertEqual(written["completed_shells"], 2)
        self.assertEqual(calls, [1, 1, 2, 2, 2])

    @unittest.skipUnless(
        os.environ.get("RUN_ED_SMOKE") == "1",
        "set RUN_ED_SMOKE=1 to run QuSpin ED smoke tests",
    )
    def test_tiny_runtime_ed_smoke_when_quspin_is_available(self) -> None:
        if not QUSPIN_PYTHON.exists():
            self.skipTest(f"QuSpin Python not found: {QUSPIN_PYTHON}")

        with tempfile.TemporaryDirectory() as tmpdir:
            output = Path(tmpdir) / "results.json"
            command = [
                str(QUSPIN_PYTHON),
                str(ED_SCRIPT),
                "--Lx",
                "3",
                "--Ly",
                "3",
                "--beta",
                "0.5",
                "--mu",
                "-5",
                "--U1",
                "0",
                "--U2",
                "1",
                "--max-total-particles",
                "1",
                "--output",
                str(output),
            ]

            completed = subprocess.run(
                command,
                cwd=REPO_ROOT,
                text=True,
                capture_output=True,
                timeout=300,
            )
            if completed.returncode != 0 and "QuSpin is required" in (
                completed.stderr + completed.stdout
            ):
                self.skipTest(completed.stderr or completed.stdout)
            self.assertEqual(
                completed.returncode,
                0,
                msg=f"stdout:\n{completed.stdout}\nstderr:\n{completed.stderr}",
            )

            self.assertTrue(output.exists())
            payload = json.loads(output.read_text(encoding="utf-8"))

        self.assertNotEqual(payload["status"], "dry_run")
        self.assertEqual(payload["completed_shells"], 1)
        self.assertEqual(set(payload["observables"]), set(self.ed.ED_OBSERVABLES))
        self.assertAlmostEqual(
            payload["observables"]["numsquare_up"],
            payload["observables"]["num_up"],
        )
        self.assertAlmostEqual(
            payload["observables"]["numsquare_do"],
            payload["observables"]["num_do"],
        )
        for name in self.ed.ED_OBSERVABLES:
            self.assertTrue(math.isfinite(float(payload["observables"][name])), name)


if __name__ == "__main__":
    unittest.main()
