"""Map the archived solver inputs to the manuscript's three benchmark scans.

The archive's ``row`` identifies a row of the original solver-specific figures.
It remains part of the provenance; ``figure_row`` identifies the current layout.
Neither the solver couplings nor the original case order/seeds are rewritten.
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
        "scan_parameter": "U", "solver_parameter": "U2",
    },
    ("pairing", 0): {
        "benchmark": "pairing_delta", "section": "main",
        "figure": "benchmark_combined", "figure_label": "fig:benchmark",
        "figure_row": 1, "panels": "e-h",
        "scan_parameter": "Delta", "solver_parameter": "Delta",
    },
    ("number_conserving", 1): {
        "benchmark": "attractive_u1", "section": "supplement",
        "figure": "benchmark_attractive", "figure_label": "fig:benchmark_attractive",
        "figure_row": 0, "panels": "a-d",
        "scan_parameter": "U1", "solver_parameter": "U1",
    },
}


def describe_case(case):
    """Return paper notation while retaining the solver's U1/U2 convention."""
    key = (case["model"], case["row"])
    if key not in BENCHMARKS:
        raise ValueError(f"unknown manuscript benchmark: {key}")
    description = dict(BENCHMARKS[key])
    params = case["parameters"]
    if description["section"] == "main" and params["U1"] != 0:
        raise ValueError(f"main-text U requires U1=0: {case['id']}")
    if case["x"] != params[description["solver_parameter"]]:
        raise ValueError(f"scan coordinate differs from solver input: {case['id']}")
    description["U"] = params["U2"] if description["section"] == "main" else None
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
