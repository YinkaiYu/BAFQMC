"""Opt-in real MPI/ED checks of the resumable production driver.

Run with BAFQMC_RUN_MPI_TESTS=1 and optionally BAFQMC_PYTHON_ED pointing to
the QuSpin interpreter. These tests use two eight-bin smoke cases only.
"""
import copy
import json
import math
import os
from pathlib import Path
import shutil
import subprocess
import sys
import tempfile
import unittest
from unittest.mock import patch

from benchmarks.paper.analysis import recompute
from benchmarks.paper.production import produce, sha, stage_path
from src.pairing.benchmarks import campaign_pairing_delta as pairing


ROOT = Path(__file__).resolve().parents[2]
MODELS = ("number_conserving", "pairing")


def prepare_smoke_data(data):
    """Keep all smoke inputs independent of the 22-point paper campaign."""
    nc_source = ROOT / "src/number_conserving/benchmarks/campaigns/pipeline_smoke.json"
    nc = json.loads(nc_source.read_text())
    for case in nc["cases"]:
        for key in ("input_dir", "ed_params"):
            case[key] = str((nc_source.parent / case[key]).resolve())
    pair_source = ROOT / "src/pairing/benchmarks/campaigns/triangle_pairing_pipeline_smoke.json"
    paired = json.loads(pair_source.read_text())
    for model, manifest in zip(MODELS, (nc, paired)):
        directory = data / model
        directory.mkdir(parents=True)
        (directory / "manifest.json").write_text(json.dumps(manifest))
    pairing.init_local(paired, data / "pairing/inputs", seed=paired["dqmc_defaults"]["seeds"][0])
    cases = []
    for model, case_id, params, kind in (
        ("number_conserving", nc["cases"][0]["id"], nc["cases"][0]["parameters"], "exact_free_boson"),
        ("pairing", pairing.case_name(next(pairing.iter_cases(paired))),
         {**next(pairing.iter_cases(paired))["parameters"], **paired["ed_defaults"]}, "fixed_cutoff_dense_trace"),
    ):
        cases.append({
            "model": model, "id": case_id, "parameters": params,
            "reference_kind": kind, "samples": 8, "block_size": 2,
            "skip_samples": 0, "observables": {}, "row": 0, "x": params.get("Delta", 0),
            "raw_archive": "smoke_has_no_raw_chain", "ed_status": "not_run",
            "ed_file": "smoke", "input_dir": "smoke",
        })
    index = {"schema_version": 1, "cases": cases}
    (data / "index.json").write_text(json.dumps(index))
    return index


@unittest.skipUnless(os.environ.get("BAFQMC_RUN_MPI_TESTS") == "1", "set BAFQMC_RUN_MPI_TESTS=1 for real MPI/ED smoke runs")
class ProductionIntegrationTests(unittest.TestCase):
    def test_both_solvers_copy_resume_and_restart_interrupted_steps(self):
        if sys.platform != "linux" or shutil.which("mpirun") is None:
            self.skipTest("Linux and mpirun are required")
        python_ed = os.environ.get("BAFQMC_PYTHON_ED", sys.executable)
        with tempfile.TemporaryDirectory(prefix="bafqmc-driver-test-") as temporary:
            base = Path(temporary)
            data, output, scratch = base / "data", base / "output", base / "scratch"
            output.mkdir()
            index = prepare_smoke_data(data)
            arguments = (ROOT, data, index, output, MODELS, ("dqmc", "ed"), python_ed, 1)
            def invoke(**kwargs):
                try:
                    produce(*arguments, work_dir=scratch, **kwargs)
                except subprocess.CalledProcessError:
                    for log in sorted(base.rglob("*.log")):
                        print(f"\n{log}:\n{log.read_text(errors='replace')[-3000:]}", file=sys.stderr)
                    raise
            with patch.dict(os.environ, {key: "1" for key in ("OMP_NUM_THREADS", "MKL_NUM_THREADS", "OPENBLAS_NUM_THREADS", "NUMBA_NUM_THREADS")}):
                invoke()
                state_file = output / "progress.json"
                original_state = json.loads(state_file.read_text())
                self.assertEqual(len(original_state["steps"]), 4)
                for case in index["cases"]:
                    for stage in ("dqmc", "ed"):
                        key = f"{case['model']}/{case['id']}/{stage}"
                        self.assertEqual(original_state["steps"][key]["status"], "complete")
                        source = stage_path(scratch, case["model"], case["id"], stage)
                        destination = stage_path(output / "runs", case["model"], case["id"], stage)
                        self.assertTrue(destination.exists())
                        if source.is_file():
                            self.assertEqual(sha(source), sha(destination))
                        for name, expected in original_state["steps"][key]["files"].items():
                            self.assertEqual(sha(output / name), expected)

                records, blocks = recompute(index, data, fresh_dir=output / "runs")
                self.assertEqual(len(records), 2)
                self.assertEqual(len(blocks), 32)
                for record in records:
                    self.assertIn("stored_data", record)
                    for values in record["observables"].values():
                        self.assertTrue(all(math.isfinite(values[key]) for key in ("dqmc", "stderr", "ed")))
                nc_result = stage_path(output / "runs", "number_conserving", index["cases"][0]["id"], "ed") / "results.json"
                self.assertEqual(json.loads(nc_result.read_text())["status"], "exact_free_boson")
                pair_result = stage_path(output / "runs", "pairing", index["cases"][1]["id"], "ed")
                self.assertEqual(json.loads(pair_result.read_text())["cutoffs"]["ncut"], 2)

                # Completed stages retain their exact records and are never relaunched.
                invoke(resume=True)
                self.assertEqual(json.loads(state_file.read_text())["steps"], original_state["steps"])
                with self.assertRaisesRegex(ValueError, "resume"):
                    invoke()

                # The same solvers with a different manuscript case selection
                # must not silently continue a different campaign in this output.
                subset = {**index, "cases": index["cases"][:1]}
                with self.assertRaisesRegex(ValueError, "run settings differ"):
                    produce(ROOT, data, subset, output, MODELS, ("dqmc", "ed"), python_ed, 1, resume=True)

                # A selected scan launches only its own case with the original
                # full manifest, retaining the campaign's input/seed ordering.
                selected_output = base / "selected-output"
                produce(ROOT, data, subset, selected_output, ("number_conserving",), ("ed",), python_ed, 1, work_dir=base / "selected-scratch")
                selected_state = json.loads((selected_output / "progress.json").read_text())
                self.assertEqual(list(selected_state["steps"]), [f"number_conserving/{subset['cases'][0]['id']}/ed"])
                self.assertFalse((selected_output / "runs/pairing").exists())

                # A changed result is rejected before it can be treated as completed.
                density = stage_path(output / "runs", "number_conserving", index["cases"][0]["id"], "dqmc") / "density_total"
                correct_bytes = density.read_bytes()
                density.write_bytes(correct_bytes + b"999\n")
                with self.assertRaisesRegex(ValueError, "completed output changed"):
                    invoke(resume=True)
                density.write_bytes(correct_bytes)

                # Simulate interruption before either selected stage was copied back.
                # The paired chain contains a partial append that must be preserved
                # in scratch/incomplete and excluded from the restarted chain.
                interrupted = copy.deepcopy(original_state)
                for model, case_id, stage in (
                    ("number_conserving", index["cases"][0]["id"], "ed"),
                    ("pairing", index["cases"][1]["id"], "dqmc"),
                ):
                    interrupted["steps"][f"{model}/{case_id}/{stage}"] = {"status": "running"}
                    destination = stage_path(output / "runs", model, case_id, stage)
                    shutil.rmtree(destination) if destination.is_dir() else destination.unlink()
                partial = stage_path(scratch, "pairing", index["cases"][1]["id"], "dqmc") / "density_total"
                partial.write_bytes(partial.read_bytes() + b"999\n")
                state_file.write_text(json.dumps(interrupted))
                invoke(resume=True)
                recovered = json.loads(state_file.read_text())
                self.assertTrue(all(step["status"] == "complete" for step in recovered["steps"].values()))
                restarted = stage_path(output / "runs", "pairing", index["cases"][1]["id"], "dqmc") / "density_total"
                self.assertEqual(len(restarted.read_text().splitlines()), 8)
                backups = list((scratch / "incomplete/pairing" / index["cases"][1]["id"]).glob("dqmc-*/density_total"))
                self.assertEqual(len(backups), 1)
                self.assertEqual(len(backups[0].read_text().splitlines()), 9)
                records, blocks = recompute(index, data, fresh_dir=output / "runs")
                self.assertEqual((len(records), len(blocks)), (2, 32))


if __name__ == "__main__":
    unittest.main()
