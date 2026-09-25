"""Check that manuscript scopes preserve the published inputs and figure meaning."""
import copy
import csv
import json
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest
from unittest.mock import patch

from benchmarks.paper.analysis import load_index, processed_cases, save_tables
from benchmarks.paper.manuscript import select_index


ROOT = Path(__file__).resolve().parents[2]
DATA = ROOT / "benchmarks/paper/data"


def case_keys(cases):
    return {(case["model"], case["id"]) for case in cases}


class ManuscriptTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.index = load_index(DATA)

    def test_published_scans_and_notation_map_to_the_same_hamiltonians(self):
        groups = {}
        for case in self.index["cases"]:
            groups.setdefault(case["manuscript"]["benchmark"], []).append(case)
        self.assertEqual(set(groups), {"relative_density_u1", "pairing_delta", "total_density_u2"})
        expected = {
            "relative_density_u1": ([0, .25, .5, .75, 1, 1.25, 1.5, 2], "U1", "combined", "benchmark_combined", 0, "a-d", -3.5, 4),
            "pairing_delta": ([0, .05, .1, .15, .2, .25, .3], "Delta", "combined", "benchmark_combined", 1, "e-h", -5, 4),
            "total_density_u2": ([-.6, -.5, -.4, -.3, -.2, -.1, 0], "U2", "total_density", "benchmark_total_density", 0, "a-d", -7, 1),
        }
        for group, (grid, scan_parameter, section, figure, row, panels, mu, beta) in expected.items():
            with self.subTest(benchmark=group):
                cases = groups[group]
                self.assertEqual([case["parameters"][scan_parameter] for case in cases], grid)
                for case in cases:
                    params, mapping = case["parameters"], case["manuscript"]
                    self.assertEqual((params["Lx"], params["Ly"], params["t"], params["mu"], params["beta"]), (3, 3, 1, mu, beta))
                    self.assertEqual((mapping["section"], mapping["figure"], mapping["figure_row"], mapping["panels"]), (section, figure, row, panels))
                    self.assertEqual(mapping["scan_parameter"], scan_parameter)
                    self.assertEqual(case["x"], params[scan_parameter])
                    if section == "combined":
                        self.assertEqual(params["U2"], 0)
                        self.assertLess(params["mu"], -3 * params["t"] - abs(params.get("Delta", 0)))
                    else:
                        self.assertEqual(params["U1"], 1)
                    if group == "relative_density_u1":
                        self.assertEqual(params["U2"], 0)
                    if group == "pairing_delta":
                        self.assertEqual((params["nmax"], params["ncut"]), (3, 4))
                    else:
                        self.assertEqual(params.get("Delta", 0), 0)

    def test_combined_and_total_density_partition_the_full_default_campaign(self):
        all_cases = select_index(self.index)["cases"]
        combined = select_index(self.index, scope="combined")["cases"]
        total_density = select_index(self.index, scope="total_density")["cases"]
        self.assertEqual((len(all_cases), len(combined), len(total_density)), (22, 15, 7))
        self.assertEqual(case_keys(combined) | case_keys(total_density), case_keys(all_cases))
        self.assertFalse(case_keys(combined) & case_keys(total_density))
        self.assertEqual(len(select_index(self.index, models=("number_conserving",), scope="combined")["cases"]), 8)
        self.assertEqual(len(select_index(self.index, models=("pairing",), scope="combined")["cases"]), 7)
        with self.assertRaises(ValueError):
            select_index(self.index, models=("pairing",), scope="total_density")

    def test_selection_preserves_source_inputs_seeds_and_numeric_records(self):
        original = json.loads((DATA / "index.json").read_text())
        before = copy.deepcopy(original)
        selected = select_index(original, scope="combined")
        source_by_key = {(case["model"], case["id"]): case for case in original["cases"]}
        for case in selected["cases"]:
            source = source_by_key[(case["model"], case["id"])]
            self.assertEqual({key: value for key, value in case.items() if key != "manuscript"}, source)
            if case["model"] == "pairing":
                seed_file = DATA / case["input_dir"] / "seeds.txt"
                self.assertEqual(seed_file, DATA / source["input_dir"] / "seeds.txt")
                self.assertEqual(len(seed_file.read_text().splitlines()), 128)
        selected["cases"][0]["parameters"]["U2"] = 99
        selected["cases"][0]["observables"]["density_total"]["dqmc"] = 99
        self.assertEqual(original, before)

    def test_processed_blocks_follow_case_selection_within_each_solver(self):
        for scope, count, block_count in (("combined", 15, 600), ("total_density", 7, 280)):
            with self.subTest(scope=scope):
                index = select_index(self.index, scope=scope)
                cases, blocks = processed_cases(index, DATA)
                self.assertEqual((len(cases), len(blocks)), (count, block_count))
                self.assertEqual({(block["model"], block["case"]) for block in blocks}, case_keys(cases))

    def test_plan_validates_scope_before_creating_output(self):
        with tempfile.TemporaryDirectory() as temporary:
            output = Path(temporary) / "not-created"
            command = [sys.executable, str(ROOT / "reproduce.py"), "--plan", "--output", str(output)]
            combined = subprocess.run(command + ["--scope", "combined"], capture_output=True, text=True)
            self.assertEqual(combined.returncode, 0, combined.stderr)
            self.assertIn("Mode: full", combined.stdout)
            self.assertIn("15", combined.stdout)
            self.assertIn("combined", combined.stdout)
            paired = subprocess.run(command + ["--scope", "combined", "--model", "pairing"], capture_output=True, text=True)
            self.assertEqual(paired.returncode, 0, paired.stderr)
            self.assertIn("benchmark_combined_pairing.pdf (e-h)", paired.stdout)
            invalid = subprocess.run(command + ["--scope", "total_density", "--model", "pairing"], capture_output=True, text=True)
            self.assertNotEqual(invalid.returncode, 0)
            self.assertFalse(output.exists())

    def test_exported_tables_carry_the_current_manuscript_coordinates(self):
        cases, blocks = processed_cases(self.index, DATA)
        with tempfile.TemporaryDirectory() as temporary:
            output = Path(temporary)
            save_tables(cases, blocks, output)
            with (output / "observables.csv").open(newline="") as handle:
                rows = list(csv.DictReader(handle))
            self.assertEqual(len(rows), 88)
            by_key = {(case["model"], case["id"]): case for case in cases}
            fields = ("benchmark", "section", "figure", "figure_row", "panels", "scan_parameter")
            for row in rows:
                mapping = by_key[(row["model"], row["case"])]["manuscript"]
                for field in fields:
                    expected = "" if mapping[field] is None else str(mapping[field])
                    self.assertEqual(row[field], expected)
            exported = json.loads((output / "records.json").read_text())
            self.assertEqual([case["manuscript"] for case in exported], [case["manuscript"] for case in cases])



if __name__ == "__main__":
    unittest.main()
