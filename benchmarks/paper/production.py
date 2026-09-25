"""Production campaign orchestration with per-case checkpoints and Linux scratch."""
from __future__ import annotations
import hashlib
import json
import os
from pathlib import Path
import shutil
import subprocess
import sys
import tempfile
import time


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def write_state(path, state):
    pending = path.with_suffix(".tmp")
    pending.write_text(json.dumps(state, indent=2) + "\n")
    pending.replace(path)


def stage_path(base, model, case_id, stage):
    if model == "number_conserving":
        return base / model / case_id / ("dqmc" if stage == "dqmc" else "ed")
    return base / model / "inputs" / case_id if stage == "dqmc" else base / model / "ed_results" / (case_id + ".json")


def hashes(path, root):
    files = [path] if path.is_file() else [p for p in path.rglob("*") if p.is_file()]
    return {str(p.relative_to(root)): sha(p) for p in sorted(files)}


def produce(root, data, index, output, models, stages, python_ed, threads, *, resume=False, work_dir=None):
    output.mkdir(parents=True, exist_ok=True)
    state_path = output / "progress.json"
    source_files = [p for p in (root / "src").rglob("*") if p.is_file() and "build" not in p.parts and (p.suffix in {".f90", ".py", ".mk"} or p.name == "Makefile")]
    source_hash = hashlib.sha256("\n".join(f"{p.relative_to(root)} {sha(p)}" for p in sorted(source_files)).encode()).hexdigest()
    package_files = [p for p in data.rglob("*") if p.is_file() and "raw" not in p.relative_to(data).parts]
    package_hash = hashlib.sha256("\n".join(f"{p.relative_to(data)} {sha(p)}" for p in sorted(package_files)).encode()).hexdigest()
    cases = [c for c in index["cases"] if c["model"] in models]
    signature = {"data_package_sha256": package_hash, "solver_source_sha256": source_hash, "models": list(models), "cases": [f"{c['model']}/{c['id']}" for c in cases], "python_ed": str(python_ed), "threads": threads, "build_environment": {k: os.environ.get(k) for k in ("FC", "FFLAGS", "LDFLAGS", "LDLIBS", "MPIEXEC")}}
    if state_path.exists():
        state = json.loads(state_path.read_text())
        if state["signature"] != signature:
            raise ValueError("run settings differ from progress.json; choose a fresh --output")
        # Separate --mode dqmc and --mode ed can use the same output without --resume.
        if not resume and any(key.rsplit("/", 1)[-1] in stages for key in state["steps"]):
            raise ValueError("this campaign already contains these stages; use --resume with the same --output")
        scratch = Path(state["scratch"])
        if work_dir is not None and scratch != work_dir.resolve():
            raise ValueError("--work-dir differs from the recorded scratch directory")
        scratch.mkdir(parents=True, exist_ok=True)
    else:
        if resume:
            raise ValueError("--resume needs an existing progress.json in --output")
        if (output / "runs").exists() and any((output / "runs").iterdir()):
            raise ValueError("existing runs without progress.json; choose a fresh --output")
        scratch = work_dir.resolve() if work_dir else Path(tempfile.mkdtemp(prefix="bafqmc-"))
        if scratch.exists() and any(scratch.iterdir()):
            raise ValueError("--work-dir must initially be empty")
        scratch.mkdir(parents=True, exist_ok=True)
        state = {"signature": signature, "scratch": str(scratch), "steps": {}}
        write_state(state_path, state)
    if "ed" in stages:
        subprocess.run([str(python_ed), "-c", "import quspin, numpy, scipy"], check=True)
    print(f"Working directory: {scratch}\nCompleted raw results are copied to: {output / 'runs'}", flush=True)
    total = len(cases) * len(stages)
    current = 0
    for model in models:
        if "dqmc" in stages:
            subprocess.run(["make", "-C", str(root / "src" / model), "build"], check=True)
        for case in [c for c in cases if c["model"] == model]:
            case_id = case["id"]
            for stage in stages:
                current += 1
                key = f"{model}/{case_id}/{stage}"
                previous = state["steps"].get(key)
                if previous and previous["status"] == "complete":
                    for name, expected in previous["files"].items():
                        if sha(output / name) != expected:
                            raise ValueError(f"completed output changed: {name}")
                    print(f"[{current}/{total}] Already complete: {key}", flush=True)
                    continue
                source = stage_path(scratch, model, case_id, stage)
                if previous and source.exists():
                    # A failed chain restarts from its manifest input, never appends.
                    backup = scratch / "incomplete" / model / case_id / f"{stage}-{time.time_ns()}"
                    backup.parent.mkdir(parents=True, exist_ok=True)
                    shutil.move(str(source), str(backup))
                    if model == "pairing" and stage == "dqmc":
                        source.mkdir(parents=True)
                        for name in ("params.json", "paramC_sets.txt", "confin.txt", "seeds.txt"):
                            shutil.copyfile(data / model / "inputs" / case_id / name, source / name)
                command = [sys.executable, str(root / "src" / model / "run_paper.py"), "--manifest", str(data / model / "manifest.json"), "--output", str(scratch / model), "--mode", stage, "--case", case_id, "--python", str(python_ed)]
                log = output / "logs" / f"{model}.{case_id}.{stage}.log"
                log.parent.mkdir(parents=True, exist_ok=True)
                if log.exists():
                    log.rename(log.with_name(log.name + f".previous-{time.time_ns()}"))
                if model == "pairing" and (scratch / model / "run_manifest.json").exists():
                    # Recover missing initial inputs without touching measured outputs.
                    for initial in (data / model / "inputs").iterdir():
                        target = scratch / model / "inputs" / initial.name
                        target.mkdir(parents=True, exist_ok=True)
                        for name in ("params.json", "paramC_sets.txt", "confin.txt", "seeds.txt"):
                            if not (target / name).exists():
                                shutil.copyfile(initial / name, target / name)
                state["steps"][key] = {"status": "running", "command": command, "log": str(log.relative_to(output))}
                write_state(state_path, state)
                print(f"[{current}/{total}] {stage.upper()}: {model}/{case_id} (log: {log})", flush=True)
                start = time.monotonic()
                with log.open("w") as handle:
                    subprocess.run(command, stdout=handle, stderr=subprocess.STDOUT, check=True)
                destination = stage_path(output / "runs", model, case_id, stage)
                destination.parent.mkdir(parents=True, exist_ok=True)
                if source.is_file():
                    shutil.copy2(source, destination)
                else:
                    shutil.copytree(source, destination, dirs_exist_ok=True)
                state["steps"][key].update(status="complete", elapsed_seconds=time.monotonic() - start, files=hashes(destination, output))
                write_state(state_path, state)
                print(f"  Completed in {state['steps'][key]['elapsed_seconds']:.1f} s", flush=True)
    state["completed_at"] = time.strftime("%Y-%m-%dT%H:%M:%S%z")
    write_state(state_path, state)
