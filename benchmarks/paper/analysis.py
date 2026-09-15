"""Reblock the archived measurements and recover the four paper observables."""
from __future__ import annotations

import copy
import csv
import io
import json
import math
from pathlib import Path
import tarfile

import numpy as np

RAW_OBSERVABLES = {
    "density_total": "density_total", "energy_density": "energy_density",
    "S_SF_K": "sf_K", "S_DW_K": "dw_K",
}


def load_index(data_dir):
    from .manuscript import select_index
    index = json.loads((data_dir / "index.json").read_text())
    if index["schema_version"] != 1 or len(index["cases"]) != 22:
        raise ValueError("expected the version-1, 22-point paper data package")
    return select_index(index)


def processed_cases(index, data_dir, *, models=("number_conserving", "pairing")):
    """Check the compact, tracked block means against the plotted means/SEM."""
    with (data_dir / "block_means.csv").open(newline="") as handle:
        blocks = list(csv.DictReader(handle))
    selected = [copy.deepcopy(c) for c in index["cases"] if c["model"] in models]
    for c in selected:
        ed = reference_values(c, data_dir / c["ed_file"])
        for name in RAW_OBSERVABLES:
            rows = [b for b in blocks if b["model"] == c["model"] and b["case"] == c["id"] and b["observable"] == name]
            n = c["samples"] // c["block_size"]
            if sorted(int(b["block"]) for b in rows) != list(range(n)):
                raise ValueError(f"missing or duplicate blocks: {c['id']}/{name}")
            means = [float(b["block_mean"]) for b in sorted(rows, key=lambda b: int(b["block"]))]
            if not all(math.isfinite(v) for v in means):
                raise ValueError(f"nonfinite blocks: {c['id']}/{name}")
            mean = sum(means) / n
            sem = math.sqrt(sum((v - mean) ** 2 for v in means) / (n - 1) / n)
            for key, value in (("dqmc", mean), ("stderr", sem), ("ed", ed[name])):
                if not math.isclose(value, c["observables"][name][key], rel_tol=5e-11, abs_tol=5e-14):
                    raise ValueError(f"processed {key} mismatch: {c['id']}/{name}")
    selected_keys = {(c["model"], c["id"]) for c in selected}
    return selected, [b for b in blocks if (b["model"], b["case"]) in selected_keys]


def block_stats(raw: bytes, *, columns: int, samples: int, block_size: int, skip_samples: int = 0):
    values = np.loadtxt(io.BytesIO(raw), ndmin=2)
    if values.shape != (samples, columns) or not np.isfinite(values).all():
        raise ValueError(f"expected {samples} finite rows x {columns} columns, got {values.shape}")
    real = values[skip_samples:, 0]
    if block_size <= 0 or len(real) % block_size or len(real) < 2 * block_size:
        raise ValueError("measurements must contain at least two complete blocks")
    # Match the original scripts' ordered Python sums, including the free point.
    means = [sum(block.tolist()) / block_size for block in real.reshape(-1, block_size)]
    mean = sum(means) / len(means)
    stderr = math.sqrt(sum((v - mean) ** 2 for v in means) / (len(means) - 1) / len(means))
    return mean, stderr, means


def reference_values(case, path):
    if case["reference_kind"] == "exact_free_boson":
        from src.number_conserving.benchmarks.campaign_analysis import free_boson_reference_observables
        return free_boson_reference_observables(case["parameters"])
    result = json.loads(path.read_text())
    for key in ("Lx", "Ly", "U1", "U2", "mu", "beta", "t", "Delta"):
        if key in case["parameters"]:
            if result.get("parameters", {}).get(key) != case["parameters"][key]:
                raise ValueError(f"ED parameter mismatch: {case['id']}/{key}")
    if case["model"] == "pairing":
        for key in ("nmax", "ncut"):
            if result.get("cutoffs", {}).get(key) != case["parameters"][key]:
                raise ValueError(f"ED cutoff mismatch: {case['id']}/{key}")
    values = result["observables"]
    values = {name: float(values[name]["value"] if isinstance(values[name], dict) else values[name]) for name in RAW_OBSERVABLES}
    if not all(math.isfinite(value) for value in values.values()):
        raise ValueError(f"nonfinite ED reference: {case['id']}")
    return values


def recompute(index, data_dir, *, fresh_dir=None, models=("number_conserving", "pairing")):
    cases, block_records = [], []
    for original in index["cases"]:
        if original["model"] not in models:
            continue
        case = copy.deepcopy(original)
        model, case_id = case["model"], case["id"]
        if fresh_dir is None:
            with tarfile.open(data_dir / case["raw_archive"], "r:gz") as archive:
                members = archive.getmembers()
                if len(members) != 4 or {m.name for m in members} != set(RAW_OBSERVABLES.values()) or any(not m.isfile() for m in members):
                    raise ValueError(f"unexpected archive members: {case_id}")
                raw = {m.name: archive.extractfile(m).read() for m in members}
            ed_path = data_dir / case["ed_file"]
        else:
            base = fresh_dir / model
            if model == "number_conserving":
                run_dir, ed_path = base / case_id / "dqmc", base / case_id / "ed/results.json"
            else:
                run_dir, ed_path = base / "inputs" / case_id, base / "ed_results" / (case_id + ".json")
            raw = {name: (run_dir / name).read_bytes() for name in RAW_OBSERVABLES.values()}
            archive_keys = ("raw_archive", "original_run_dir", "ed_status", "ed_file", "input_dir")
            case["archived_source"] = {key: case.pop(key) for key in archive_keys}
            case["run_dir"] = str(run_dir)
            case["ed_file"] = str(ed_path)
            case["ed_status"] = "exact_free_boson" if case["reference_kind"] == "exact_free_boson" else json.loads(ed_path.read_text()).get("status", "fixed_cutoff_dense_trace")
        ed = reference_values(case, ed_path)
        for name, filename in RAW_OBSERVABLES.items():
            mean, stderr, blocks = block_stats(raw[filename], columns=2 if filename in {"sf_K", "dw_K"} else 1, samples=case["samples"], block_size=case["block_size"], skip_samples=case["skip_samples"])
            if fresh_dir is None:
                for key, value in (("dqmc", mean), ("stderr", stderr), ("ed", ed[name])):
                    if not math.isclose(value, case["observables"][name][key], rel_tol=5e-11, abs_tol=5e-14):
                        raise ValueError(f"archived {key} mismatch: {model}/{case_id}/{name}: {value} vs {case['observables'][name][key]}")
            case["observables"][name] = {"dqmc": mean, "stderr": stderr, "ed": ed[name]}
            for i, value in enumerate(blocks):
                block_records.append({"model": model, "case": case_id, "observable": name, "block": i, "block_mean": value})
        cases.append(case)
    return cases, block_records


def save_tables(cases, blocks, output):
    output.mkdir(parents=True, exist_ok=True)
    rows = []
    for c in cases:
        for name, value in c["observables"].items():
            rows.append({"model": c["model"], "case": c["id"], **c.get("manuscript", {}), "row": c["row"], "x": c["x"], "observable": name, "dqmc": value["dqmc"], "stderr": value["stderr"], "ed": value["ed"], "difference": value["dqmc"] - value["ed"]})
    for filename, records in (("observables.csv", rows), ("block_means.csv", blocks)):
        if not records:
            continue
        with (output / filename).open("w", newline="") as handle:
            writer = csv.DictWriter(handle, fieldnames=records[0], lineterminator="\n")
            writer.writeheader()
            writer.writerows(records)
    (output / "records.json").write_text(json.dumps(cases, indent=2, allow_nan=False) + "\n")


def plot(cases, output):
    import matplotlib
    matplotlib.use("Agg")
    from .plot_manuscript import plot_figures
    plot_figures(cases, output)
