#!/usr/bin/env python3
"""Tests for the 3x3 observable benchmark campaign helpers."""

from __future__ import annotations

import json
import tempfile
import unittest
from pathlib import Path

from src.number_conserving.benchmarks import campaign_3x3 as campaign


class Campaign3x3Tests(unittest.TestCase):
    def write_matching_outputs(self, run_dir: Path) -> None:
        for name, values in {
            "num_up": ["2.0", "2.0", "2.0", "2.0"],
            "num_do": ["2.0", "2.0", "2.0", "2.0"],
            "density_total": ["4.0", "4.0", "4.0", "4.0"],
            "energy_density": ["-1.0", "-1.0", "-1.0", "-1.0"],
            "doubleOcc": ["0.25", "0.25", "0.25", "0.25"],
            "onsite_n2_up": ["2.0", "2.0", "2.0", "2.0"],
            "onsite_n2_do": ["2.0", "2.0", "2.0", "2.0"],
        }.items():
            (run_dir / name).write_text("\n".join(values) + "\n", encoding="utf-8")
        for name in ("sf_K", "psf_Gamma", "dw_K"):
            (run_dir / name).write_text(
                "1.0 0.0\n1.0 0.0\n1.0 0.0\n1.0 0.0\n",
                encoding="utf-8",
            )
        (run_dir / "density_site_total").write_text(
            "\n".join([" ".join(["4.0"] * 9)] * 4) + "\n",
            encoding="utf-8",
        )

        observables = {
            "density_total": 4.0,
            "total_kinetic": -9.0,
            "doubleOcc": 0.25,
            "onsite_n2_up": 2.0,
            "onsite_n2_do": 2.0,
            "interaction_energy_density": 0.0,
            "energy_density": -1.0,
            "S_SF_K": 1.0,
            "S_PSF_Gamma": 1.0,
            "S_DW_K": 1.0,
            "num_up": 2.0,
            "num_do": 2.0,
        }
        (run_dir / "results.json").write_text(
            json.dumps(
                {
                    "schema_version": 1,
                    "parameters": {
                        "Lx": 3,
                        "Ly": 3,
                        "beta": 4.0,
                        "mu": -3.5,
                        "U1": 0.0,
                        "U2": 0.0,
                        "max_total_particles": 6,
                    },
                    "status": "converged",
                    "completed_shells": 3,
                    "observables": observables,
                    "relative_last_shell": {
                        name: 1.0e-4 for name in observables
                    },
                    "convergence_policy": {
                        "tail_tolerance": 1.0e-3,
                        "max_density_total": 1.0,
                        "min_completed_shells": 2,
                    },
                }
            ),
            encoding="utf-8",
        )


    def test_iter_cases_expands_both_sweeps(self) -> None:
        cases = list(campaign.iter_cases(campaign.default_manifest()))
        labels = {case["label"] for case in cases}

        self.assertIn("U2_0_mu-3.5_b4", labels)
        self.assertIn("U1_0_mu-5_b1", labels)
        self.assertEqual({case["sweep"] for case in cases}, {"U2_sweep", "U1_sweep"})
        self.assertTrue(
            any(case["U2"] == 0.0 and case["mu"] == -3.5 for case in cases)
        )
        self.assertTrue(
            any(case["U1"] == 0.0 and case["U2"] == 1.0 for case in cases)
        )



    def test_write_dqmc_case_creates_fixed_runtime_inputs(self) -> None:
        manifest = campaign.default_manifest()
        case = {
            "sweep": "U2_sweep",
            "Lx": 3,
            "Ly": 3,
            "beta": 4.0,
            "mu": -3.5,
            "U1": 0.0,
            "U2": 1.5,
        }
        with tempfile.TemporaryDirectory() as tmp:
            run_dir = campaign.write_dqmc_case(
                Path(tmp),
                case,
                manifest["dqmc_defaults"],
                seed=13579,
            )
            names = {path.name for path in run_dir.iterdir()}
            param_text = (run_dir / "paramC_sets.txt").read_text(encoding="utf-8")
            confin_text = (run_dir / "confin.txt").read_text(encoding="utf-8")
            seeds_text = (run_dir / "seeds.txt").read_text(encoding="utf-8")

        self.assertGreaterEqual(names, {"paramC_sets.txt", "confin.txt", "seeds.txt"})
        self.assertEqual(confin_text, "0\n")
        self.assertEqual(seeds_text, "13579\n")
        rows = [line.split() for line in param_text.splitlines()[:7]]
        self.assertEqual([float(x) for x in rows[0]], [0.0, 1.5, -3.5])
        self.assertEqual([float(x) for x in rows[1]], [3, 3, 400, 4.0])
        self.assertEqual([float(x) for x in rows[3]], [10, 100000, 1, 1.5])
        self.assertEqual(rows[4][0].lower(), ".false.")
        self.assertEqual(int(rows[4][1]), 0)
        self.assertEqual(rows[5][0].lower(), ".true.")
        self.assertEqual([float(x) for x in rows[5][1:]], [500, 1.0, 1.0])
        self.assertEqual([float(x) for x in rows[6]], [2, 0.1, 0.0, 0.0])

    def test_write_dqmc_case_removes_stale_output_files_on_regeneration(self) -> None:
        manifest = campaign.default_manifest()
        case = next(iter(campaign.iter_cases(manifest)))
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            run_dir = campaign.write_dqmc_case(
                root, case, manifest["dqmc_defaults"], seed=13579
            )
            (run_dir / "density_total").write_text("old\n", encoding="utf-8")
            (run_dir / "density_site_total").write_text("old\n", encoding="utf-8")
            (run_dir / "sf_K").write_text("old\n", encoding="utf-8")
            (run_dir / "notes.txt").write_text("keep\n", encoding="utf-8")

            rerun_dir = campaign.write_dqmc_case(
                root, case, manifest["dqmc_defaults"], seed=24680
            )

            self.assertEqual(rerun_dir, run_dir)
            self.assertFalse((run_dir / "density_total").exists())
            self.assertFalse((run_dir / "density_site_total").exists())
            self.assertFalse((run_dir / "sf_K").exists())
            self.assertEqual(
                (run_dir / "notes.txt").read_text(encoding="utf-8"), "keep\n"
            )
            self.assertEqual(
                (run_dir / "seeds.txt").read_text(encoding="utf-8"), "24680\n"
            )

    def test_write_ed_case_creates_script_compatible_json(self) -> None:
        manifest = campaign.default_manifest()
        case = next(iter(campaign.iter_cases(manifest)))
        with tempfile.TemporaryDirectory() as tmp:
            run_dir = campaign.write_ed_case(Path(tmp), case, manifest["ed_policy"])
            params = json.loads((run_dir / "params.json").read_text(encoding="utf-8"))

        self.assertEqual(run_dir.name, campaign.case_name(case))
        self.assertEqual(params["parameters"]["Lx"], 3)
        self.assertEqual(params["parameters"]["Ly"], 3)
        self.assertEqual(params["parameters"]["max_total_particles"], 6)
        self.assertEqual(params["convergence_policy"]["tail_tolerance"], 1.0e-3)
        self.assertNotIn("max_total_particles", params["convergence_policy"])

    def test_init_local_writes_dqmc_and_ed_inputs(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            output_dir = Path(tmp) / "inputs"
            status = campaign.main(
                ["init-local", "--output-dir", str(output_dir), "--seed", "13579"]
            )
            case = next(iter(campaign.iter_cases(campaign.default_manifest())))
            run_dir = output_dir / campaign.case_name(case)

            names = {path.name for path in run_dir.iterdir()}

        self.assertEqual(status, 0)
        self.assertGreaterEqual(
            names, {"paramC_sets.txt", "confin.txt", "seeds.txt", "params.json"}
        )


    def test_results_summary_retains_valid_measurements(self) -> None:
        manifest = campaign.default_manifest()
        manifest["sweeps"] = [{**manifest["sweeps"][0], "U2": [0.0]}]
        manifest["dqmc_defaults"]["block_size"] = 2
        case = next(iter(campaign.iter_cases(manifest)))
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            input_dir = root / "inputs"
            output_dir = root / "summary"
            run_dir = campaign.write_dqmc_case(
                input_dir,
                case,
                manifest["dqmc_defaults"],
                seed=13579,
            )
            campaign.write_ed_case(input_dir, case, manifest["ed_policy"])
            self.write_matching_outputs(run_dir)
            paths = campaign.write_results_summary(input_dir, output_dir, manifest)

            payload = json.loads(paths["json"].read_text(encoding="utf-8"))
            self.assertTrue(payload["cases"][0]["trusted"])
            self.assertEqual(payload["counts"], {"trusted": 1})


if __name__ == "__main__":
    unittest.main()
