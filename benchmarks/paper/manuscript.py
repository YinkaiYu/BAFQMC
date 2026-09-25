"""Describe the current benchmark layout in the paper's notation.

The records already use the paper convention directly: ``U1`` is the
repulsive relative-density coupling and ``U2`` is the attractive total-density
coupling.
"""
from __future__ import annotations

import copy


MODELS = ("number_conserving", "pairing")
SCOPES = ("all", "main", "supplement")
BENCHMARKS = {
    ("number_conserving", 0): {
        "benchmark": "repulsive_u", "section": "main",
        "figure": "benchmark_combined", "figure_label": "fig:benchmark",
        "figure_row": 0, "panels": "a-d",
        "scan_parameter": "U1",
        "fixed_parameters": {"U2": 0.0},
    },
    ("pairing", 0): {
        "benchmark": "pairing_delta", "section": "main",
        "figure": "benchmark_combined", "figure_label": "fig:benchmark",
        "figure_row": 1, "panels": "e-h",
        "scan_parameter": "Delta",
        "fixed_parameters": {"U1": 1.0, "U2": 0.0},
    },
    ("number_conserving", 1): {
        "benchmark": "attractive_u2", "section": "supplement",
        "figure": "benchmark_attractive", "figure_label": "fig:benchmark_attractive",
        "figure_row": 0, "panels": "a-d",
        "scan_parameter": "U2",
        "fixed_parameters": {"U1": 1.0},
    },
}


def describe_case(case):
    """Return the canonical manuscript notation for one benchmark case."""
    key = (case["model"], case["row"])
    if key not in BENCHMARKS:
        raise ValueError(f"unknown manuscript benchmark: {key}")
    description = dict(BENCHMARKS[key])
    params = case["parameters"]
    fixed_parameters = description.pop("fixed_parameters")
    for name, expected in fixed_parameters.items():
        if params[name] != expected:
            raise ValueError(f"{case['id']} has {name}={params[name]}, expected {expected}")
    if case["x"] != params[description["scan_parameter"]]:
        raise ValueError(f"scan coordinate differs from parameter: {case['id']}")
    return description


def select_index(index, *, models=MODELS, scope="all"):
    """Select paper/SM cases without modifying manifests or seed ordering."""
    if scope not in SCOPES:
        raise ValueError(f"unknown benchmark scope: {scope}")
    selected = []
    for original in index["cases"]:
        if original["model"] not in models:
            continue
        description = describe_case(original)
        if scope != "all" and description["section"] != scope:
            continue
        case = copy.deepcopy(original)
        case["manuscript"] = description
        selected.append(case)
    if not selected:
        raise ValueError(f"no benchmark points for scope={scope}, models={','.join(models)}")
    return {**index, "cases": selected}


def output_stem(cases, figure):
    """Name individual main-text rows explicitly when one solver is selected."""
    models = {c["model"] for c in cases if c["manuscript"]["figure"] == figure}
    if figure == "benchmark_combined" and len(models) == 1:
        return figure + "_" + next(iter(models))
    return figure


def selection_summary(cases):
    """Human-readable counts in manuscript reading order."""
    lines = []
    for description in BENCHMARKS.values():
        count = sum(c["manuscript"]["benchmark"] == description["benchmark"] for c in cases)
        if count:
            lines.append(
                f"  {description['section']}: {description['scan_parameter']} scan, "
                f"{count} points -> {output_stem(cases, description['figure'])}.pdf ({description['panels']})"
            )
    return "\n".join(lines)
