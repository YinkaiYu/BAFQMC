"""Integration checks for the small reader-facing data package."""
import math
from pathlib import Path
import shutil
import subprocess
import sys
import tempfile
import unittest

from benchmarks.paper.analysis import block_stats, load_index, processed_cases, verify_checksums

ROOT = Path(__file__).resolve().parents[2]
DATA = ROOT / "benchmarks/paper/data"


class ReproductionTests(unittest.TestCase):
    def test_default_is_full_and_plan_does_not_launch_or_create_output(self):
        with tempfile.TemporaryDirectory() as temp:
            output = Path(temp) / "not-created"
            result = subprocess.run([sys.executable, str(ROOT / "reproduce.py"), "--plan", "--output", str(output)], check=True, capture_output=True, text=True)
            self.assertIn("Mode: full", result.stdout)
            self.assertIn("12-24 hours", result.stdout)
            self.assertFalse(output.exists())

    def test_clean_data_package_without_raw_measurements(self):
        with tempfile.TemporaryDirectory() as temp:
            data = Path(temp) / "data"
            shutil.copytree(DATA, data, ignore=shutil.ignore_patterns("raw"))
            self.assertFalse(list(data.rglob("*.tar.gz")))
            verify_checksums(data)
            cases, blocks = processed_cases(load_index(data), data)
            self.assertEqual(len(cases), 22)
            self.assertEqual(len(blocks), 880)
            self.assertEqual(sum(c["model"] == "number_conserving" for c in cases), 15)
            self.assertEqual(sum(c["reference_kind"] == "exact_free_boson" for c in cases), 1)

    def test_standard_error_uses_independent_blocks(self):
        mean, sem, blocks = block_stats(b"1\n3\n5\n7\n", columns=1, samples=4, block_size=2)
        self.assertEqual(blocks, [2, 6])
        self.assertEqual(mean, 4)
        self.assertEqual(sem, 2)
        self.assertNotEqual(sem, math.sqrt(20 / 3))  # raw-sample standard deviation

    def test_partial_or_nonfinite_measurements_are_rejected(self):
        for raw in (b"1\n3\n5\n", b"1\n3\nnan\n7\n"):
            with self.assertRaises(ValueError):
                block_stats(raw, columns=1, samples=4, block_size=2)

    def test_modified_processed_data_is_detected(self):
        with tempfile.TemporaryDirectory() as temp:
            data = Path(temp) / "data"
            shutil.copytree(DATA, data, ignore=shutil.ignore_patterns("raw"))
            path = data / "block_means.csv"
            path.write_text(path.read_text() + "\n")
            with self.assertRaisesRegex(ValueError, "checksum mismatch"):
                verify_checksums(data)


if __name__ == "__main__":
    unittest.main()
