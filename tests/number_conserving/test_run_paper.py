"""Regression checks for reference routing and ED resource guards."""
import json
import math
import os
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest

from src.number_conserving.benchmarks.ed import EDtriangle_quspin_3x3 as ed

ROOT = Path(__file__).resolve().parents[2] / "src/number_conserving"


class PaperRunnerTests(unittest.TestCase):
    def test_free_reference_bypasses_cutoff_ed(self):
        with tempfile.TemporaryDirectory() as temp:
            root = Path(temp)
            inputs = root / 'inputs'
            inputs.mkdir()
            params = {'Lx': 3, 'Ly': 3, 'U1': 0., 'U2': 0., 'beta': 4.,
                      'mu': -3.5, 't': 1., 'max_total_particles': 10}
            (inputs / 'params.json').write_text(json.dumps({'parameters': params}))
            manifest = root / 'manifest.json'
            manifest.write_text(json.dumps({'cases': [{'id': 'free', 'input_dir': 'inputs',
                                                       'ed_params': 'inputs/params.json'}]}))
            subprocess.run([sys.executable, str(ROOT / 'run_paper.py'), '--manifest', str(manifest),
                            '--output', str(root / 'output'), '--mode', 'ed',
                            '--python', '/nonexistent/quspin/python'], check=True, capture_output=True)
            result = json.loads((root / 'output/free/ed/results.json').read_text())
            self.assertEqual(result['status'], 'exact_free_boson')
            expected = 2/9*(2/math.expm1(2)+6/math.expm1(14)+1/math.expm1(38))
            self.assertAlmostEqual(result['observables']['density_total'], expected, places=14)

    def test_dense_memory_budget_preserves_paper_shell_seven(self):
        self.assertLess(ed.assert_dense_memory_budget(9075, 12.), 5.)
        with self.assertRaisesRegex(MemoryError, '27225'):
            ed.assert_dense_memory_budget(27225, 12.)

    @unittest.skipUnless(os.environ.get('RUN_ED_SMOKE') == '1', 'set RUN_ED_SMOKE=1 for QuSpin guard runtime')
    def test_guard_preserves_completed_shell_checkpoint(self):
        python = os.environ.get('ED_PYTHON', sys.executable)
        with tempfile.TemporaryDirectory() as temp:
            output = Path(temp) / 'results.json'
            command = [python, str(ROOT / 'benchmarks/ed/EDtriangle_quspin_3x3.py'),
                       '--Lx', '3', '--Ly', '3', '--U1', '0', '--U2', '1',
                       '--beta', '.5', '--mu', '-5', '--max-total-particles', '1',
                       '--output', str(output)]
            rejected = subprocess.run(command + ['--dense-memory-cap-gib', '1e-9'],
                                      capture_output=True, text=True, timeout=60)
            self.assertNotEqual(rejected.returncode, 0)
            self.assertIn('last completed particle shell', rejected.stderr)
            self.assertEqual(json.loads(output.read_text())['completed_shells'], 0)
            accepted = subprocess.run(command + ['--dense-memory-cap-gib', '12'],
                                      capture_output=True, text=True, timeout=60)
            self.assertEqual(accepted.returncode, 0, accepted.stderr)
            self.assertEqual(json.loads(output.read_text())['completed_shells'], 1)
