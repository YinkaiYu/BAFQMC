#!/usr/bin/env python3
from __future__ import annotations

import importlib
import json
import tempfile
import unittest
from collections.abc import Mapping
from pathlib import Path
from typing import Any


EXPECTED_DELTAS = [0.0, 0.1, 0.2, 0.3, 0.4]
EXPECTED_SHARED_PARAMS = {
    "Lx": 3,
    "Ly": 3,
    "U1": 0.0,
    "U2": 1.0,
    "beta": 4.0,
    "mu": -3.5,
}
REQUIRED_INPUTS = ("paramC_sets.txt", "confin.txt", "seeds.txt")
PARAMC_DATA_ROWS = 7
KNOWN_APPEND_OUTPUTS = (
    "density",
    "density_up",
    "density_do",
    "density_total",
    "density_site_total",
    "num_up",
    "num_do",
    "kinetic",
    "doubleOcc",
    "squareOcc",
    "local_numsquare",
    "numsquare_up",
    "numsquare_do",
    "pair_equal",
    "den_upup_sub11",
    "den_dodo_sub11",
    "den_updo",
    "energy_density",
    "interaction_energy_density",
    "pairing_energy_density",
    "chemical_energy_density",
    "grand_energy_density",
    "onsite_n2_up",
    "onsite_n2_do",
    "sf_K",
    "dw_K",
    "psf_Gamma",
)


def _case_parameters(case: Any) -> Mapping[str, Any]:
    if isinstance(case, Mapping):
        for key in ("parameters", "params", "dqmc_parameters"):
            if key in case and isinstance(case[key], Mapping):
                return case[key]
        return case
    for key in ("parameters", "params", "dqmc_parameters"):
        value = getattr(case, key, None)
        if isinstance(value, Mapping):
            return value
    raise AssertionError(f"cannot find parameters in case {case!r}")


def _param(case: Any, key: str) -> Any:
    params = _case_parameters(case)
    aliases = {
        "Delta": ("Delta", "delta", "RDelta"),
        "Lx": ("Lx", "Nlx"),
        "Ly": ("Ly", "Nly"),
        "beta": ("beta", "Beta"),
    }
    for candidate in aliases.get(key, (key,)):
        if candidate in params:
            return params[candidate]
    raise AssertionError(f"missing parameter {key!r} in {params!r}")


def _has_label_text(line: str) -> bool:
    allowed_logicals = {"t", "f", "true", "false"}
    for token in line.replace("#", " ").split():
        try:
            float(token.replace("D", "E").replace("d", "e"))
            continue
        except ValueError:
            pass
        normalized = token.strip(".").lower()
        if normalized in allowed_logicals:
            continue
        if any(ch.isalpha() for ch in token):
            return True
    return False


class PairingCampaignTests(unittest.TestCase):
    def setUp(self) -> None:
        self.campaign = importlib.import_module("src.pairing.benchmarks.campaign_pairing_delta")

    def _default_cases(self) -> list[Any]:
        manifest = self.campaign.default_manifest()
        if hasattr(self.campaign, "iter_cases"):
            return list(self.campaign.iter_cases(manifest))
        if isinstance(manifest, Mapping) and isinstance(manifest.get("cases"), list):
            return list(manifest["cases"])
        raise AssertionError(f"default_manifest did not expose cases: {manifest!r}")


    def test_write_dqmc_case_paramc_first_data_row_has_ru1_ru2_mu_rdelta(self) -> None:
        case = self._default_cases()[0]

        with tempfile.TemporaryDirectory() as tmpdir:
            run_dir = Path(tmpdir) / "case"
            self.campaign.write_dqmc_case(case, run_dir)

            lines = (run_dir / "paramC_sets.txt").read_text(encoding="utf-8").splitlines()
            nonempty_lines = [line.strip() for line in lines if line.strip()]
            non_comment_lines = [
                line for line in nonempty_lines if not line.lstrip().startswith("#")
            ]
            first_data_row = non_comment_lines[0]
            tokens = first_data_row.split()
            self.assertEqual(len(tokens), 4)
            try:
                values = [float(value) for value in tokens]
            except ValueError as exc:
                raise AssertionError(
                    "first nonempty non-comment paramC_sets.txt line must be numeric"
                ) from exc

            self.assertEqual(len(values), 4)
            self.assertEqual(
                values,
                [
                    float(_param(case, "U1")),
                    float(_param(case, "U2")),
                    float(_param(case, "mu")),
                    float(_param(case, "Delta")),
                ],
            )
            data_rows_seen = 0
            data_block_last_index = -1
            for index, line in enumerate(nonempty_lines):
                if line.lstrip().startswith("#"):
                    if data_rows_seen < PARAMC_DATA_ROWS:
                        self.assertFalse(_has_label_text(line), line)
                    continue
                data_rows_seen += 1
                if data_rows_seen <= PARAMC_DATA_ROWS:
                    self.assertFalse(_has_label_text(line), line)
                    data_block_last_index = index

            self.assertGreaterEqual(data_rows_seen, PARAMC_DATA_ROWS)
            for line in nonempty_lines[: data_block_last_index + 1]:
                self.assertFalse(_has_label_text(line), line)

    def test_write_dqmc_case_refreshes_inputs_and_removes_known_outputs_only(self) -> None:
        case = self._default_cases()[0]

        with tempfile.TemporaryDirectory() as tmpdir:
            run_dir = Path(tmpdir) / "case"
            run_dir.mkdir()
            for filename in KNOWN_APPEND_OUTPUTS:
                (run_dir / filename).write_text("stale\n", encoding="utf-8")
            user_file = run_dir / "user_notes.txt"
            user_file.write_text("keep this file\n", encoding="utf-8")

            self.campaign.write_dqmc_case(case, run_dir)

            for filename in REQUIRED_INPUTS:
                self.assertTrue((run_dir / filename).is_file(), filename)
            for filename in KNOWN_APPEND_OUTPUTS:
                self.assertFalse((run_dir / filename).exists(), filename)
            self.assertEqual(user_file.read_text(encoding="utf-8"), "keep this file\n")

    def test_case_local_defaults_override_manifest_defaults(self) -> None:
        manifest = self.campaign.default_manifest(Nbin=100000, nmax=2, ncut=6)
        case = self._default_cases()[0]
        case["dqmc"] = {"Nbin": 7}
        case["ed"] = {"nmax": 3, "ncut": None, "basis_cap": 12345}

        with tempfile.TemporaryDirectory() as tmpdir:
            run_dir = Path(tmpdir) / "case"
            self.campaign.write_dqmc_case(case, run_dir, manifest["dqmc_defaults"])
            self.campaign.write_ed_case(case, run_dir, manifest["ed_defaults"])

            paramc_rows = [
                line.strip().split()
                for line in (run_dir / "paramC_sets.txt").read_text(encoding="utf-8").splitlines()
                if line.strip()
            ]
            self.assertEqual(paramc_rows[3][1], "7")
            ed_params = json.loads((run_dir / "params.json").read_text(encoding="utf-8"))
            self.assertEqual(ed_params["nmax"], 3)
            self.assertIsNone(ed_params["ncut"])
            self.assertEqual(ed_params["basis_cap"], 12345)

    def test_cli_reports_malformed_delta_list_without_traceback(self) -> None:
        with tempfile.TemporaryDirectory() as tmpdir:
            with self.assertRaises(SystemExit) as raised:
                self.campaign.main(
                    [
                        "write-manifest",
                        "--output",
                        str(Path(tmpdir) / "manifest.json"),
                        "--deltas",
                        "",
                    ]
                )
            self.assertEqual(raised.exception.code, 2)

    def test_results_summary_marks_single_k_imaginary_failure(self) -> None:
        manifest = self.campaign.default_manifest(
            deltas=[0.0],
            block_size=2,
            Nbin=4,
        )
        case = next(iter(self.campaign.iter_cases(manifest)))

        with tempfile.TemporaryDirectory() as tmpdir:
            input_dir = Path(tmpdir) / "inputs"
            output_dir = Path(tmpdir) / "summary"
            run_dir = input_dir / self.campaign.case_name(case)
            run_dir.mkdir(parents=True)
            self._write_minimal_dqmc_summary_run(run_dir, lq=9)
            for filename in ("sf_K", "dw_K"):
                (run_dir / filename).write_text(
                    "\n".join(["1.0 10.0"] * 4) + "\n",
                    encoding="utf-8",
                )

            outputs = self.campaign.write_results_summary(
                input_dir,
                output_dir,
                manifest,
                block_size=2,
            )

            payload = json.loads(outputs["json"].read_text(encoding="utf-8"))
            record = payload["cases"][0]
            self.assertEqual(record["status"], "unreliable")
            self.assertFalse(record["reliable"])
            self.assertEqual(record["S_SF_K_imag_mean"], 10.0)
            self.assertEqual(record["S_SF_K_pass_imaginary"], False)
            self.assertEqual(record["S_SF_K_imaginary_diagnostic_only"], False)

    def _write_minimal_dqmc_summary_run(self, run_dir: Path, lq: int) -> None:
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
