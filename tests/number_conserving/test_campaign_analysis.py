#!/usr/bin/env python3
"""Tests for 3x3 BAFQMC-vs-ED campaign analysis helpers."""

from __future__ import annotations

import json
import math
import tempfile
import unittest
from pathlib import Path

from src.number_conserving.benchmarks import campaign_analysis as ca


class CampaignAnalysisTests(unittest.TestCase):
    def ed_payload(
        self,
        *,
        u1: float,
        status: str,
        density_total: float = 0.1,
        completed_shells: int = 3,
        tail: float = 1.0e-4,
    ) -> dict:
        observables = {name: 0.1 for name in ca.ED_REFERENCE_OBSERVABLES}
        observables["density_total"] = density_total
        return {
            "schema_version": 1,
            "parameters": {"U1": u1},
            "status": status,
            "completed_shells": completed_shells,
            "observables": observables,
            "relative_last_shell": {
                name: tail for name in ca.ED_REFERENCE_OBSERVABLES
            },
            "convergence_policy": {
                "tail_tolerance": 1.0e-3,
                "max_density_total": 1.0,
                "min_completed_shells": 2,
            },
        }

    def matching_ed_payload(self, *, status: str = "converged") -> dict:
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
        return {
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
            "status": status,
            "completed_shells": 3,
            "observables": observables,
            "relative_last_shell": {
                name: 1.0e-4 for name in ca.ED_REFERENCE_OBSERVABLES
            },
            "convergence_policy": {
                "tail_tolerance": 1.0e-3,
                "max_density_total": 1.0,
                "min_completed_shells": 2,
            },
        }

    def write_matching_dqmc_outputs(self, root: Path) -> None:
        for name, values in {
            "num_up": ["2.0", "2.0", "2.0", "2.0"],
            "num_do": ["2.0", "2.0", "2.0", "2.0"],
            "density_total": ["4.0", "4.0", "4.0", "4.0"],
            "energy_density": ["-1.0", "-1.0", "-1.0", "-1.0"],
            "doubleOcc": ["0.25", "0.25", "0.25", "0.25"],
            "onsite_n2_up": ["2.0", "2.0", "2.0", "2.0"],
            "onsite_n2_do": ["2.0", "2.0", "2.0", "2.0"],
        }.items():
            (root / name).write_text("\n".join(values) + "\n", encoding="utf-8")
        for name in ("sf_K", "psf_Gamma", "dw_K"):
            (root / name).write_text(
                "1.0 0.0\n1.0 0.0\n1.0 0.0\n1.0 0.0\n",
                encoding="utf-8",
            )
        self.write_density_site_total(root, rows=4)

    def write_density_site_total(
        self, root: Path, *, rows: int, profile: list[float] | None = None
    ) -> None:
        values = profile or [4.0] * 9
        line = " ".join(str(value) for value in values)
        (root / "density_site_total").write_text(
            "\n".join([line] * rows) + "\n",
            encoding="utf-8",
        )

    def test_read_scalar_and_complex_series(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            (root / "density_total").write_text("1.0\n2.5\n", encoding="utf-8")
            (root / "sf_K").write_text("1.0 0.1\n2.0 -0.2\n", encoding="utf-8")

            self.assertEqual(ca.read_scalar_series(root / "density_total"), [1.0, 2.5])
            complex_values = ca.read_complex_series(root / "sf_K")
            self.assertEqual(complex_values, [complex(1.0, 0.1), complex(2.0, -0.2)])

            (root / "malformed_complex").write_text("1.0\n", encoding="utf-8")
            with self.assertRaises(ValueError):
                ca.read_complex_series(root / "malformed_complex")

    def test_read_scalar_series_rejects_non_finite_values(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            (root / "density_total").write_text("nan\n", encoding="utf-8")

            with self.assertRaisesRegex(ValueError, "non-finite.*line 1"):
                ca.read_scalar_series(root / "density_total")

    def test_read_complex_series_rejects_non_finite_values(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            (root / "sf_K_real").write_text("nan 0.0\n", encoding="utf-8")
            (root / "sf_K_imag").write_text("1.0 inf\n", encoding="utf-8")

            with self.assertRaisesRegex(ValueError, "non-finite.*line 1"):
                ca.read_complex_series(root / "sf_K_real")
            with self.assertRaisesRegex(ValueError, "non-finite.*line 1"):
                ca.read_complex_series(root / "sf_K_imag")

    def test_block_statistics_uses_block_means(self) -> None:
        stats = ca.block_statistics([1.0, 3.0, 5.0, 7.0], block_size=2)
        self.assertEqual(stats["blocks"], 2)
        self.assertEqual(stats["samples_used"], 4)
        self.assertAlmostEqual(stats["actual"], 4.0)
        self.assertAlmostEqual(stats["stderr"], math.sqrt(8.0 / 2.0))

    def test_block_statistics_for_complex_series_reports_real_and_imag(self) -> None:
        stats = ca.complex_block_statistics(
            [1 + 1j, 3 + 3j, 5 + 1j, 7 + 3j], block_size=2
        )
        self.assertAlmostEqual(stats["real"]["actual"], 4.0)
        self.assertAlmostEqual(stats["imag"]["actual"], 2.0)
        self.assertEqual(stats["real"]["blocks"], 2)
        self.assertEqual(stats["imag"]["blocks"], 2)

    def test_k_point_index_requires_lattice_multiple_of_three(self) -> None:
        self.assertEqual(ca.k_point_cell_index(3, 3), (3, 2))
        self.assertEqual(ca.k_point_cell_index(6, 6), (5, 3))
        with self.assertRaisesRegex(ValueError, "multiples of 3"):
            ca.k_point_cell_index(3, 2)


    def test_ipr_uses_normalized_density_profile(self) -> None:
        value = ca.compute_ipr([1.0, 1.0, 1.0, 1.0], min_number=1.0e-8)
        self.assertAlmostEqual(value.value, 0.25)
        self.assertTrue(value.reliable)

        localized = ca.compute_ipr([0.0, 0.0, 2.0, 0.0], min_number=1.0e-8)
        self.assertAlmostEqual(localized.value, 1.0)
        self.assertTrue(localized.reliable)

        unstable = ca.compute_ipr([1.0e-12, 0.0, 0.0], min_number=1.0e-8)
        self.assertFalse(unstable.reliable)
        self.assertTrue(math.isnan(unstable.value))

    def test_ed_reliability_accepts_bounded_and_negative_u1_statuses(self) -> None:
        bounded = self.ed_payload(u1=0.0, status="converged")
        self.assertTrue(ca.ed_reference_reliable(bounded).reliable)

        negative = self.ed_payload(
            u1=-0.1,
            status="low_density_cutoff_accepted",
            density_total=0.01,
            completed_shells=6,
            tail=0.5,
        )
        decision = ca.ed_reference_reliable(negative)
        self.assertTrue(decision.reliable)
        self.assertEqual(decision.kind, "finite_window")

        incomplete = self.ed_payload(u1=0.0, status="incomplete", tail=0.2)
        incomplete_decision = ca.ed_reference_reliable(incomplete)
        self.assertFalse(incomplete_decision.reliable)
        self.assertEqual(
            incomplete_decision.reason,
            "unsupported_ed_status: status=incomplete, U1=0.0",
        )

    def test_ed_reliability_rejects_minimal_or_stale_payloads(self) -> None:
        minimal = {"parameters": {"U1": 0.0}, "status": "converged"}

        decision = ca.ed_reference_reliable(minimal)

        self.assertFalse(decision.reliable)
        self.assertEqual(decision.reason, "unsupported_or_missing_schema_version")

    def test_ed_reliability_cross_checks_tail_and_cutoff_policy(self) -> None:
        bad_tail = self.ed_payload(u1=0.0, status="converged", tail=0.2)
        too_few_converged_shells = self.ed_payload(
            u1=0.0,
            status="converged",
            completed_shells=1,
            tail=1.0e-5,
        )
        high_density = self.ed_payload(
            u1=-0.1,
            status="low_density_cutoff_accepted",
            density_total=2.0,
            tail=0.5,
        )
        too_few_shells = self.ed_payload(
            u1=-0.1,
            status="low_density_cutoff_accepted",
            completed_shells=1,
            tail=0.5,
        )
        too_shallow_for_trusted_window = self.ed_payload(
            u1=-0.2,
            status="low_density_cutoff_accepted",
            density_total=0.032,
            completed_shells=2,
            tail=1.0,
        )

        self.assertEqual(ca.ed_reference_reliable(bad_tail).reason, "tail_not_converged")
        self.assertEqual(
            ca.ed_reference_reliable(too_few_converged_shells).reason,
            "insufficient_completed_shells",
        )
        self.assertEqual(
            ca.ed_reference_reliable(high_density).reason,
            "density_exceeds_cutoff",
        )
        self.assertEqual(
            ca.ed_reference_reliable(too_few_shells).reason,
            "insufficient_completed_shells",
        )
        self.assertEqual(
            ca.ed_reference_reliable(too_shallow_for_trusted_window).reason,
            "insufficient_low_density_window_shells",
        )

    def test_ed_reliability_rejects_schema_only_dry_run_payloads(self) -> None:
        dry_run = self.ed_payload(u1=0.0, status="dry_run")
        dry_run["parameters"]["dry_run"] = True

        decision = ca.ed_reference_reliable(dry_run)

        self.assertFalse(decision.reliable)
        self.assertEqual(decision.reason, "dry_run_payload")

    def test_ed_target_observables_compute_postprocessed_ipr(self) -> None:
        targets = ca.ed_target_observables(self.matching_ed_payload())

        self.assertAlmostEqual(targets["IPR"]["actual"], 1.0 / 9.0)
        self.assertTrue(targets["IPR"]["reliable"])
        self.assertAlmostEqual(targets["S_SF_K"]["actual"], 1.0)

    def test_compare_observable_handles_unreliable_nan_fail_closed(self) -> None:
        comparison = ca.compare_observable(
            {"actual": math.nan, "stderr": math.nan, "reliable": False, "reason": "bad"},
            {"actual": 1.0, "reliable": True},
        )

        self.assertFalse(comparison["reliable"])
        self.assertEqual(comparison["reason"], "bad")
        self.assertFalse(comparison["pass_uncertainty"])

    def test_compare_observable_rejects_large_structure_factor_imaginary_part(self) -> None:
        comparison = ca.compare_observable(
            {
                "actual": 1.0,
                "stderr": 0.1,
                "imag_actual": 1.0,
                "imag_stderr": 0.1,
            },
            {"actual": 1.0, "reliable": True},
        )

        self.assertFalse(comparison["reliable"])
        self.assertFalse(comparison["pass_imaginary"])
        self.assertEqual(comparison["reason"], "imaginary_part_not_zero")

    def test_compare_observable_accepts_roundoff_imaginary_part(self) -> None:
        comparison = ca.compare_observable(
            {
                "actual": 1.0,
                "stderr": 0.1,
                "imag_actual": 1.0e-21,
                "imag_stderr": 1.0e-24,
            },
            {"actual": 1.0, "reliable": True},
        )

        self.assertTrue(comparison["reliable"])
        self.assertTrue(comparison["pass_imaginary"])
        self.assertEqual(comparison["reason"], "ok")

    def test_density_comparison_rejects_non_positive_dqmc_density(self) -> None:
        dqmc_observable = ca._dqmc_observable_for_comparison(
            "density_total",
            {
                "actual": -1.0,
                "stderr": 100.0,
            },
        )
        comparison = ca.compare_observable(
            dqmc_observable,
            {"actual": 0.1, "reliable": True},
        )

        self.assertFalse(comparison["reliable"])
        self.assertEqual(comparison["reason"], "non_positive_density")

    def test_ed_result_matches_case_rejects_stale_parameters(self) -> None:
        payload = self.matching_ed_payload()
        payload["parameters"] = {
            "Lx": 3,
            "Ly": 3,
            "beta": 4.0,
            "mu": -3.0,
            "U1": 0.0,
            "U2": 0.0,
            "max_total_particles": 6,
        }

        decision = ca.ed_result_matches_case(
            payload,
            {
                "Lx": 3,
                "Ly": 3,
                "beta": 4.0,
                "mu": -3.5,
                "U1": 0.0,
                "U2": 0.0,
            },
            {"max_total_particles": 6},
        )

        self.assertFalse(decision.reliable)
        self.assertIn("mu", decision.reason)

    def test_compare_dqmc_ed_case_accepts_matching_reliable_payload(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            self.write_matching_dqmc_outputs(root)
            (root / "results.json").write_text(
                json.dumps(self.matching_ed_payload()),
                encoding="utf-8",
            )

            summary = ca.compare_dqmc_ed_case(root, block_size=2, lq=9)

        self.assertTrue(summary["trusted"])
        self.assertEqual(summary["status"], "trusted")
        self.assertTrue(summary["observables"]["IPR"]["reliable"])

    def test_compare_dqmc_ed_case_rejects_incomplete_ed_payload(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            self.write_matching_dqmc_outputs(root)
            (root / "results.json").write_text(
                json.dumps(self.matching_ed_payload(status="incomplete")),
                encoding="utf-8",
            )

            summary = ca.compare_dqmc_ed_case(root, block_size=2, lq=9)

        self.assertFalse(summary["trusted"])
        self.assertEqual(summary["status"], "ed_unreliable")
        self.assertIn("unsupported_ed_status", summary["reason"])

    def test_compare_dqmc_ed_case_uses_exact_reference_for_free_boson_case(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            exact = ca.free_boson_reference_observables(
                {"Lx": 3, "Ly": 3, "beta": 4.0, "mu": -3.5, "U1": 0.0, "U2": 0.0}
            )
            self.write_matching_dqmc_outputs(root)
            for name in ("density_total", "energy_density", "doubleOcc"):
                value = exact[name]
                (root / name).write_text(
                    "\n".join([repr(value)] * 4) + "\n",
                    encoding="utf-8",
                )
            for filename, observable in {
                "sf_K": "S_SF_K",
                "psf_Gamma": "S_PSF_Gamma",
                "dw_K": "S_DW_K",
            }.items():
                value = exact[observable]
                (root / filename).write_text(
                    "\n".join([f"{repr(value)} 0.0"] * 4) + "\n",
                    encoding="utf-8",
                )
            site_value = exact["density_total"]
            self.write_density_site_total(root, rows=4, profile=[site_value] * 9)
            (root / "results.json").write_text(
                json.dumps(self.matching_ed_payload(status="incomplete")),
                encoding="utf-8",
            )

            summary = ca.compare_dqmc_ed_case(
                root,
                block_size=2,
                lq=9,
                expected_parameters={
                    "Lx": 3,
                    "Ly": 3,
                    "beta": 4.0,
                    "mu": -3.5,
                    "U1": 0.0,
                    "U2": 0.0,
                },
            )

        self.assertTrue(summary["trusted"])
        self.assertEqual(summary["status"], "trusted")
        self.assertEqual(summary["ed_reliability"]["kind"], "exact_free_boson")
        self.assertEqual(summary["ed"]["status"], "exact_free_boson")
        self.assertAlmostEqual(
            summary["observables"]["density_total"]["ed"],
            0.06956450548352416,
        )

    def test_free_boson_reference_reports_flavor_number_second_moments(self) -> None:
        exact = ca.free_boson_reference_observables(
            {"Lx": 3, "Ly": 3, "beta": 4.0, "mu": -3.5, "U1": 0.0, "U2": 0.0}
        )

        self.assertIn("numsquare_up", exact)
        self.assertIn("numsquare_do", exact)
        self.assertGreater(exact["numsquare_up"], exact["num_up"] ** 2)
        self.assertAlmostEqual(exact["numsquare_up"], exact["numsquare_do"])

    def test_summarize_dqmc_run_computes_target_observables(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            for name, values in {
                "num_up": ["2.0", "2.0"],
                "num_do": ["2.0", "2.0"],
                "density_total": ["4.0", "4.0"],
                "energy_density": ["-1.0", "-1.0"],
                "doubleOcc": ["0.25", "0.25"],
                "onsite_n2_up": ["2.0", "2.0"],
                "onsite_n2_do": ["2.0", "2.0"],
            }.items():
                (root / name).write_text("\n".join(values) + "\n", encoding="utf-8")
            for name in ("sf_K", "psf_Gamma", "dw_K"):
                (root / name).write_text("1.0 0.0\n1.0 0.0\n", encoding="utf-8")
            self.write_density_site_total(root, rows=2)

            summary = ca.summarize_dqmc_run(root, block_size=1, lq=9)

        self.assertEqual(set(summary["observables"]), ca.TARGET_OBSERVABLES)
        self.assertAlmostEqual(summary["observables"]["IPR"]["actual"], 1.0 / 9.0)
        self.assertAlmostEqual(summary["observables"]["S_SF_K"]["actual"], 1.0)

    def test_summarize_dqmc_run_ipr_uses_blocked_sample_window(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            for name, values in {
                "num_up": ["2.0", "2.0", "2.0", "2.0", "100.0"],
                "num_do": ["2.0", "2.0", "2.0", "2.0", "100.0"],
                "density_total": ["4.0", "4.0", "4.0", "4.0", "200.0"],
                "energy_density": ["-1.0", "-1.0", "-1.0", "-1.0", "100.0"],
                "doubleOcc": ["0.25", "0.25", "0.25", "0.25", "100.0"],
                "onsite_n2_up": ["2.0", "2.0", "2.0", "2.0", "100.0"],
                "onsite_n2_do": ["2.0", "2.0", "2.0", "2.0", "100.0"],
            }.items():
                (root / name).write_text("\n".join(values) + "\n", encoding="utf-8")
            for name in ("sf_K", "psf_Gamma", "dw_K"):
                (root / name).write_text(
                    "1.0 0.0\n1.0 0.0\n1.0 0.0\n1.0 0.0\n100.0 0.0\n",
                    encoding="utf-8",
                )
            self.write_density_site_total(root, rows=4)
            with (root / "density_site_total").open("a", encoding="utf-8") as handle:
                handle.write(" ".join(["100.0"] * 9) + "\n")

            summary = ca.summarize_dqmc_run(root, block_size=2, lq=9)

        self.assertAlmostEqual(summary["observables"]["IPR"]["actual"], 1.0 / 9.0)
        self.assertEqual(summary["observables"]["IPR"]["samples_used"], 4)

    def test_summarize_dqmc_run_ipr_reports_block_stderr(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            for name, values in {
                "num_up": ["2.0", "2.0", "4.0", "4.0"],
                "num_do": ["4.0", "4.0", "6.0", "6.0"],
                "density_total": ["6.0", "6.0", "10.0", "10.0"],
                "energy_density": ["-1.0", "-1.0", "-2.0", "-2.0"],
                "doubleOcc": ["0.25", "0.25", "0.5", "0.5"],
                "onsite_n2_up": ["2.0", "2.0", "16.0", "16.0"],
                "onsite_n2_do": ["8.0", "8.0", "18.0", "18.0"],
            }.items():
                (root / name).write_text("\n".join(values) + "\n", encoding="utf-8")
            for name in ("sf_K", "psf_Gamma", "dw_K"):
                (root / name).write_text(
                    "1.0 0.0\n1.0 0.0\n2.0 0.0\n2.0 0.0\n",
                    encoding="utf-8",
                )
            (root / "density_site_total").write_text(
                "1.0 1.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0\n"
                "1.0 1.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0\n"
                "2.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0\n"
                "2.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0\n",
                encoding="utf-8",
            )

            summary = ca.summarize_dqmc_run(root, block_size=2, lq=9)

        ipr = summary["observables"]["IPR"]
        self.assertAlmostEqual(ipr["actual"], 0.625)
        self.assertAlmostEqual(ipr["stderr"], 0.25)
        self.assertTrue(math.isfinite(ipr["stderr"]))

    def test_summarize_dqmc_run_ipr_rejects_small_block_denominator(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            for name, values in {
                "num_up": ["1.0e-12", "1.0e-12", "2.0", "2.0"],
                "num_do": ["2.0", "2.0", "2.0", "2.0"],
                "density_total": ["2.0", "2.0", "4.0", "4.0"],
                "energy_density": ["-1.0", "-1.0", "-1.0", "-1.0"],
                "doubleOcc": ["0.25", "0.25", "0.25", "0.25"],
                "onsite_n2_up": ["2.0", "2.0", "2.0", "2.0"],
                "onsite_n2_do": ["2.0", "2.0", "2.0", "2.0"],
            }.items():
                (root / name).write_text("\n".join(values) + "\n", encoding="utf-8")
            for name in ("sf_K", "psf_Gamma", "dw_K"):
                (root / name).write_text(
                    "1.0 0.0\n1.0 0.0\n1.0 0.0\n1.0 0.0\n",
                    encoding="utf-8",
                )
            (root / "density_site_total").write_text(
                "1.0e-12 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0\n"
                "1.0e-12 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0\n"
                "4.0 4.0 4.0 4.0 4.0 4.0 4.0 4.0 4.0\n"
                "4.0 4.0 4.0 4.0 4.0 4.0 4.0 4.0 4.0\n",
                encoding="utf-8",
            )

            summary = ca.summarize_dqmc_run(root, block_size=2, lq=9)

        ipr = summary["observables"]["IPR"]
        self.assertFalse(ipr["reliable"])
        self.assertTrue(math.isnan(ipr["actual"]))
        self.assertTrue(math.isnan(ipr["stderr"]))
        self.assertIn("block", ipr["reason"])

    def test_summarize_dqmc_run_rejects_mismatched_output_lengths(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            for name, values in {
                "num_up": ["2.0", "2.0"],
                "num_do": ["2.0"],
                "density_total": ["4.0", "4.0"],
                "energy_density": ["-1.0", "-1.0"],
                "doubleOcc": ["0.25", "0.25"],
                "onsite_n2_up": ["2.0", "2.0"],
                "onsite_n2_do": ["2.0", "2.0"],
            }.items():
                (root / name).write_text("\n".join(values) + "\n", encoding="utf-8")
            for name in ("sf_K", "psf_Gamma", "dw_K"):
                (root / name).write_text("1.0 0.0\n1.0 0.0\n", encoding="utf-8")
            self.write_density_site_total(root, rows=2)

            with self.assertRaisesRegex(ValueError, "num_up=2.*num_do=1"):
                ca.summarize_dqmc_run(root, block_size=1, lq=9)

    def test_summarize_dqmc_run_rejects_complex_output_length_mismatch(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            for name, values in {
                "num_up": ["2.0", "2.0"],
                "num_do": ["2.0", "2.0"],
                "density_total": ["4.0", "4.0"],
                "energy_density": ["-1.0", "-1.0"],
                "doubleOcc": ["0.25", "0.25"],
                "onsite_n2_up": ["2.0", "2.0"],
                "onsite_n2_do": ["2.0", "2.0"],
            }.items():
                (root / name).write_text("\n".join(values) + "\n", encoding="utf-8")
            (root / "sf_K").write_text("1.0 0.0\n", encoding="utf-8")
            for name in ("psf_Gamma", "dw_K"):
                (root / name).write_text("1.0 0.0\n1.0 0.0\n", encoding="utf-8")
            self.write_density_site_total(root, rows=2)

            with self.assertRaisesRegex(ValueError, "sf_K=1.*psf_Gamma=2"):
                ca.summarize_dqmc_run(root, block_size=1, lq=9)

    def test_summarize_dqmc_run_rejects_non_positive_lq(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            for name, values in {
                "num_up": ["2.0", "2.0"],
                "num_do": ["2.0", "2.0"],
                "density_total": ["4.0", "4.0"],
                "energy_density": ["-1.0", "-1.0"],
                "doubleOcc": ["0.25", "0.25"],
                "onsite_n2_up": ["2.0", "2.0"],
                "onsite_n2_do": ["2.0", "2.0"],
            }.items():
                (root / name).write_text("\n".join(values) + "\n", encoding="utf-8")
            for name in ("sf_K", "psf_Gamma", "dw_K"):
                (root / name).write_text("1.0 0.0\n1.0 0.0\n", encoding="utf-8")

            with self.assertRaisesRegex(ValueError, "lq must be positive"):
                ca.summarize_dqmc_run(root, block_size=1, lq=0)


if __name__ == "__main__":
    unittest.main()
