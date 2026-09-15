#!/usr/bin/env python3
from __future__ import annotations

import importlib
import json
import math
import tempfile
import unittest
from collections.abc import Mapping
from pathlib import Path
from typing import Any


def _sample_sem(values: list[float]) -> float:
    average = sum(values) / float(len(values))
    variance = sum((value - average) ** 2 for value in values) / float(len(values) - 1)
    return math.sqrt(variance / float(len(values)))


def _field(container: Any, *names: str) -> Any:
    if isinstance(container, Mapping):
        for name in names:
            if name in container:
                return container[name]
    for name in names:
        if hasattr(container, name):
            return getattr(container, name)
    raise AssertionError(f"missing any of fields {names!r} in {container!r}")


def _stat_mean(stats: Any) -> float:
    return float(_field(stats, "mean", "actual", "value"))


def _stat_sem(stats: Any) -> float:
    return float(_field(stats, "sem", "stderr"))


def _ipr_value_and_reliability(result: Any) -> tuple[float, bool]:
    if isinstance(result, (float, int)):
        return float(result), True
    return float(_field(result, "ipr", "value")), bool(
        _field(result, "reliable", "is_reliable")
    )


def _issues_text(container: Any) -> str:
    parts: list[str] = []
    names = ("issues", "errors", "warnings", "reason", "message", "notes", "status")
    if isinstance(container, Mapping):
        parts.extend(str(container[name]) for name in names if name in container)
    for name in names:
        if hasattr(container, name):
            parts.append(str(getattr(container, name)))
    return " ".join(parts)


def _is_unreliable_or_unavailable(container: Any) -> bool:
    if isinstance(container, Mapping):
        for name in ("reliable", "is_reliable", "ok", "success", "available", "is_available", "valid"):
            if name in container:
                return not bool(container[name])
        status = str(container.get("status", "")).lower()
        if any(word in status for word in ("unavailable", "unreliable", "invalid")):
            return True
    for name in ("reliable", "is_reliable", "ok", "success", "available", "is_available", "valid"):
        if hasattr(container, name):
            return not bool(getattr(container, name))
    status = str(getattr(container, "status", "")).lower()
    return any(word in status for word in ("unavailable", "unreliable", "invalid"))


class PairingDeltaAnalysisTests(unittest.TestCase):
    def setUp(self) -> None:
        self.analysis = importlib.import_module("src.pairing.benchmarks.pairing_delta_analysis")

    def test_block_statistics_reports_sem_from_block_means(self) -> None:
        samples = [1.0, 3.0, 5.0, 7.0, 10.0, 14.0, 18.0, 22.0]
        block_means = [2.0, 6.0, 12.0, 20.0]

        stats = self.analysis.block_statistics(samples, block_size=2)

        self.assertAlmostEqual(_stat_mean(stats), sum(block_means) / len(block_means))
        self.assertAlmostEqual(_stat_sem(stats), _sample_sem(block_means))
        self.assertNotAlmostEqual(_stat_sem(stats), _sample_sem(samples))
        self.assertEqual(int(_field(stats, "blocks", "n_blocks")), len(block_means))

    def test_complex_block_statistics_separates_real_and_imaginary_columns(self) -> None:
        rows = [
            (1.0, 0.0),
            (3.0, 2.0),
            (5.0, 4.0),
            (7.0, 6.0),
            (10.0, 8.0),
            (14.0, 10.0),
            (18.0, 12.0),
            (22.0, 14.0),
        ]

        stats = self.analysis.complex_block_statistics(rows, block_size=2)

        real_stats = _field(stats, "real", "real_stats")
        imag_stats = _field(stats, "imag", "imaginary", "imag_stats")
        self.assertAlmostEqual(_stat_mean(real_stats), 10.0)
        self.assertAlmostEqual(_stat_sem(real_stats), _sample_sem([2.0, 6.0, 12.0, 20.0]))
        self.assertAlmostEqual(_stat_mean(imag_stats), 7.0)
        self.assertAlmostEqual(_stat_sem(imag_stats), _sample_sem([1.0, 5.0, 9.0, 13.0]))

    def test_summarize_dqmc_run_reads_required_complex_observable_files(self) -> None:
        with tempfile.TemporaryDirectory() as tmpdir:
            run_dir = Path(tmpdir)
            for filename in (
                "density_total",
                "energy_density",
                "interaction_energy_density",
                "pairing_energy_density",
                "chemical_energy_density",
                "grand_energy_density",
                "onsite_n2_up",
                "onsite_n2_do",
            ):
                (run_dir / filename).write_text("1.0\n3.0\n5.0\n7.0\n", encoding="utf-8")
            (run_dir / "density_site_total").write_text(
                "\n".join(["1 1 1 1 1 1 1 1 1"] * 4) + "\n",
                encoding="utf-8",
            )
            for filename, offset in (("sf_K", 0.0), ("dw_K", 10.0), ("psf_Gamma", 20.0)):
                (run_dir / filename).write_text(
                    "\n".join(
                        f"{offset + value:.1f} {0.5 * value:.1f}"
                        for value in (1.0, 3.0, 5.0, 7.0)
                    )
                    + "\n",
                    encoding="utf-8",
                )

            summary = self.analysis.summarize_dqmc_run(run_dir, block_size=2)

            for filename in ("sf_K", "dw_K", "psf_Gamma"):
                complex_summary = _field(summary, filename)
                real_stats = _field(complex_summary, "real", "real_stats")
                imag_stats = _field(complex_summary, "imag", "imaginary", "imag_stats")
                self.assertAlmostEqual(_stat_sem(real_stats), _sample_sem([2.0, 6.0]))
                self.assertAlmostEqual(_stat_sem(imag_stats), _sample_sem([1.0, 3.0]))

            for filename in (
                "energy_density",
                "interaction_energy_density",
                "pairing_energy_density",
                "chemical_energy_density",
                "grand_energy_density",
            ):
                scalar_stats = _field(summary, filename)
                self.assertAlmostEqual(_stat_mean(scalar_stats), 4.0)
                self.assertAlmostEqual(_stat_sem(scalar_stats), _sample_sem([2.0, 6.0]))

    def test_malformed_density_site_total_fails_closed_for_lq9(self) -> None:
        malformed_rows = {
            "wrong_column_count": "1 1 1 1 1 1 1 1\n",
            "overlong_column_count": "1 1 1 1 1 1 1 1 1 1\n",
            "nan_value": "1 1 nan 1 1 1 1 1 1\n",
            "inf_value": "1 1 inf 1 1 1 1 1 1\n",
        }
        for name, density_site_total in malformed_rows.items():
            with self.subTest(name=name), tempfile.TemporaryDirectory() as tmpdir:
                run_dir = Path(tmpdir)
                self._write_minimal_summary_run(run_dir, lq=9)
                (run_dir / "density_site_total").write_text(
                    density_site_total,
                    encoding="utf-8",
                )

                try:
                    summary = self.analysis.summarize_dqmc_run(
                        run_dir,
                        Lq=9,
                        block_size=2,
                    )
                except (FileNotFoundError, ValueError) as exc:
                    self.assertIn("density_site_total", str(exc))
                    continue

                self.assertFalse(
                    bool(_field(summary, "reliable", "is_reliable", "ok", "success"))
                )
                self.assertIn("density_site_total", _issues_text(summary))

    def test_incompatible_k_geometry_marks_k_observables_unavailable(self) -> None:
        with tempfile.TemporaryDirectory() as tmpdir:
            run_dir = Path(tmpdir)
            self._write_minimal_summary_run(run_dir, lq=4)
            for filename in ("sf_K", "dw_K"):
                (run_dir / filename).write_text(
                    "\n".join(["0.0 0.0"] * 4) + "\n",
                    encoding="utf-8",
                )

            summary = self.analysis.summarize_dqmc_run(
                run_dir,
                Lx=2,
                Ly=2,
                block_size=2,
            )

            for filename in ("sf_K", "dw_K"):
                observable_summary = _field(summary, filename)
                self.assertTrue(
                    _is_unreliable_or_unavailable(observable_summary),
                    f"{filename} should be unavailable for non-3-compatible geometry",
                )

    def test_summarize_dqmc_run_marks_single_k_imaginary_part_unreliable(self) -> None:
        with tempfile.TemporaryDirectory() as tmpdir:
            run_dir = Path(tmpdir)
            self._write_minimal_summary_run(run_dir, lq=9)
            for filename in ("sf_K", "dw_K"):
                (run_dir / filename).write_text(
                    "\n".join(["1.0 10.0"] * 4) + "\n",
                    encoding="utf-8",
                )

            summary = self.analysis.summarize_dqmc_run(
                run_dir,
                Lx=3,
                Ly=3,
                block_size=2,
            )

            self.assertFalse(bool(_field(summary, "reliable", "is_reliable", "ok", "success")))
            for filename in ("S_SF_K", "S_DW_K"):
                observable = _field(_field(summary, "observables"), filename)
                self.assertFalse(bool(_field(observable, "reliable", "is_reliable")))
                self.assertFalse(bool(_field(observable, "pass_imaginary")))
                self.assertFalse(bool(_field(observable, "imaginary_diagnostic_only")))
                self.assertIn("imaginary_part_not_zero", _issues_text(observable))

    def test_summarize_dqmc_run_marks_negative_k_structure_factor_unreliable(self) -> None:
        with tempfile.TemporaryDirectory() as tmpdir:
            run_dir = Path(tmpdir)
            self._write_minimal_summary_run(run_dir, lq=9)
            for filename in ("sf_K", "dw_K"):
                (run_dir / filename).write_text(
                    "\n".join(["-1.0 0.0"] * 4) + "\n",
                    encoding="utf-8",
                )

            summary = self.analysis.summarize_dqmc_run(
                run_dir,
                Lx=3,
                Ly=3,
                block_size=2,
            )

            self.assertFalse(bool(_field(summary, "reliable", "is_reliable", "ok", "success")))
            for filename in ("S_SF_K", "S_DW_K"):
                observable = _field(_field(summary, "observables"), filename)
                self.assertFalse(bool(_field(observable, "reliable", "is_reliable")))
                self.assertFalse(bool(_field(observable, "pass_nonnegative")))
                self.assertIn("negative_structure_factor", _issues_text(observable))

    def test_explicit_required_files_mark_single_k_imaginary_unreliable(self) -> None:
        with tempfile.TemporaryDirectory() as tmpdir:
            run_dir = Path(tmpdir)
            for filename in ("sf_K", "dw_K"):
                (run_dir / filename).write_text(
                    "\n".join(["1.0 10.0"] * 4) + "\n",
                    encoding="utf-8",
                )
            (run_dir / "psf_Gamma").write_text(
                "\n".join(["1.0 0.0"] * 4) + "\n",
                encoding="utf-8",
            )

            summary = self.analysis.summarize_dqmc_run(
                run_dir,
                required_files=("sf_K", "dw_K", "psf_Gamma"),
                block_size=2,
            )

            self.assertFalse(bool(_field(summary, "reliable", "is_reliable", "ok", "success")))
            for filename in ("S_SF_K", "S_DW_K"):
                observable = _field(_field(summary, "observables"), filename)
                self.assertFalse(bool(_field(observable, "reliable", "is_reliable")))
                self.assertFalse(bool(_field(observable, "pass_imaginary")))
                self.assertFalse(bool(_field(observable, "imaginary_diagnostic_only")))

    def test_summarize_dqmc_run_marks_bad_physics_diagnostics_unreliable(self) -> None:
        with tempfile.TemporaryDirectory() as tmpdir:
            run_dir = Path(tmpdir)
            self._write_minimal_summary_run(run_dir, lq=9)
            (run_dir / "density_total").write_text(
                "\n".join(["-1.0"] * 4) + "\n",
                encoding="utf-8",
            )
            (run_dir / "psf_Gamma").write_text(
                "\n".join(["1.0 10.0"] * 4) + "\n",
                encoding="utf-8",
            )

            summary = self.analysis.summarize_dqmc_run(
                run_dir,
                Lx=3,
                Ly=3,
                block_size=2,
            )

            self.assertFalse(bool(_field(summary, "reliable", "is_reliable", "ok", "success")))
            self.assertIn("density_total", _issues_text(summary))
            self.assertIn("imag", _issues_text(summary).lower())
            self.assertTrue(_is_unreliable_or_unavailable(_field(summary, "density_total")))
            self.assertTrue(_is_unreliable_or_unavailable(_field(summary, "psf_Gamma")))

    def test_compute_ipr_uniform_3x3_and_near_zero_reliability(self) -> None:
        ipr, reliable = _ipr_value_and_reliability(self.analysis.compute_ipr([1.0] * 9))
        self.assertTrue(reliable)
        self.assertAlmostEqual(ipr, 1.0 / 9.0)

        _, reliable = _ipr_value_and_reliability(self.analysis.compute_ipr([1.0e-16] * 9))
        self.assertFalse(reliable)

        _, reliable = _ipr_value_and_reliability(self.analysis.compute_ipr([-1.0] * 9))
        self.assertFalse(reliable)

    def test_ed_target_observables_exposes_campaign_observable_set(self) -> None:
        observables = {
            "density_total": {"value": 1.0},
            "energy_density": {"value": -2.0},
            "doubleOcc": {"value": 0.1},
            "interaction_energy_density": {"value": 0.2},
            "pairing_energy_density": {"value": 0.3},
            "chemical_energy_density": {"value": 0.4},
            "grand_energy_density": {"value": 0.5},
            "onsite_n2_up": {"value": 0.6},
            "onsite_n2_do": {"value": 0.7},
            "S_SF_K": {"value": 0.8},
            "S_DW_K": {"value": 0.9},
            "S_PSF_Gamma": {"value": 1.0},
        }

        targets = self.analysis.ed_target_observables({"observables": observables})

        for name in observables:
            self.assertIn(name, targets)
            self.assertAlmostEqual(float(_field(targets[name], "actual", "value")), float(observables[name]["value"]))
        self.assertNotIn(
            "IPR",
            self.analysis.ed_target_observables(
                {
                    "parameters": {"Lx": 3, "Ly": 3},
                    "observables": {"density_total": {"value": 1.0}},
                }
            ),
        )

    def test_compare_dqmc_ed_case_marks_unreliable_ed_payload_untrusted(self) -> None:
        with tempfile.TemporaryDirectory() as tmpdir:
            run_dir = Path(tmpdir)
            self._write_minimal_summary_run(run_dir, lq=9)
            observables = {
                "density_total": {"value": 1.0},
                "energy_density": {"value": 1.0},
                "doubleOcc": {"value": 1.0},
                "S_SF_K": {"value": 1.0},
                "S_DW_K": {"value": 1.0},
                "S_PSF_Gamma": {"value": 1.0},
            }
            (run_dir / "results.json").write_text(
                '{"status":"cutoff_unreliable","reliable":false,"observables":'
                + json.dumps(observables)
                + "}\n",
                encoding="utf-8",
            )

            comparison = self.analysis.compare_dqmc_ed_case(
                run_dir,
                Lx=3,
                Ly=3,
                block_size=2,
            )

            self.assertFalse(bool(_field(comparison, "trusted", "reliable", "ok", "success")))
            self.assertIn("ed", _issues_text(comparison).lower())

    def test_compare_observable_reports_z_score_for_readable_tables(self) -> None:
        comparison = self.analysis._compare_observable(
            {"actual": 2.0, "stderr": 0.25},
            {"actual": 1.0},
            stderr_tolerance=3.0,
            atol=0.0,
            rtol=0.0,
        )

        self.assertAlmostEqual(comparison["difference"], 1.0)
        self.assertAlmostEqual(comparison["z_score"], 4.0)

    def test_summarize_dqmc_run_fails_closed_when_required_file_is_missing(self) -> None:
        with tempfile.TemporaryDirectory() as tmpdir:
            run_dir = Path(tmpdir)
            (run_dir / "density_total").write_text("1.0\n1.1\n", encoding="utf-8")

            try:
                summary = self.analysis.summarize_dqmc_run(
                    run_dir,
                    required_files=("density_total", "energy_density"),
                )
            except (FileNotFoundError, ValueError) as exc:
                self.assertIn("energy_density", str(exc))
                return

            self.assertFalse(bool(_field(summary, "reliable", "is_reliable", "ok", "success")))
            self.assertIn("energy_density", str(_field(summary, "issues", "errors", "warnings")))

    def _write_minimal_summary_run(self, run_dir: Path, lq: int) -> None:
        for filename in (
            "density_total",
            "energy_density",
            "interaction_energy_density",
            "pairing_energy_density",
            "chemical_energy_density",
            "grand_energy_density",
            "onsite_n2_up",
            "onsite_n2_do",
        ):
            (run_dir / filename).write_text("1.0\n3.0\n5.0\n7.0\n", encoding="utf-8")
        (run_dir / "density_site_total").write_text(
            "\n".join([" ".join(["1.0"] * lq)] * 4) + "\n",
            encoding="utf-8",
        )
        for filename in ("sf_K", "dw_K", "psf_Gamma"):
            (run_dir / filename).write_text(
                "\n".join(["1.0 0.0"] * 4) + "\n",
                encoding="utf-8",
            )


if __name__ == "__main__":
    unittest.main()
