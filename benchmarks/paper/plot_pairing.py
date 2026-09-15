"""Original paired paper layout, with an explicit output directory.
Source: code_bosonDQMC_paring make_production_figures.py; hashes in data/index.json.
"""
from __future__ import annotations
from pathlib import Path
from typing import Any, Callable
import numpy as np
from . import plot_style
NSITE = 9

def comp(case: dict[str, Any], observable: str) -> dict[str, Any]:
    return case["observables"][observable]

def values(
    cases: list[dict[str, Any]],
    observable: str,
    transform: Callable[[float], float] = lambda value: value,
    stderr_scale: float = 1.0,
) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    x = np.array([float(case["Delta"]) for case in cases], dtype=float)
    ed = np.array([transform(float(comp(case, observable)["ed"])) for case in cases], dtype=float)
    dqmc = np.array([transform(float(comp(case, observable)["dqmc"])) for case in cases], dtype=float)
    stderr = np.array([abs(stderr_scale) * float(comp(case, observable)["stderr"]) for case in cases], dtype=float)
    return x, ed, dqmc, stderr

def y_limits(*series: np.ndarray, log: bool = False) -> tuple[float, float] | None:
    finite = np.concatenate([np.asarray(values, dtype=float).ravel() for values in series])
    finite = finite[np.isfinite(finite)]
    if finite.size == 0:
        return None
    if log:
        positive = finite[finite > 0.0]
        if positive.size == 0:
            return None
        return float(positive.min() * 0.72), float(positive.max() * 1.38)
    ymin = float(finite.min())
    ymax = float(finite.max())
    if ymin == ymax:
        span = abs(ymax) if ymax else 1.0
    else:
        span = ymax - ymin
    return ymin - 0.12 * span, ymax + 0.16 * span

def should_log(ed: np.ndarray, dqmc: np.ndarray) -> bool:
    finite = np.concatenate([ed, dqmc])
    finite = finite[np.isfinite(finite)]
    if finite.size == 0 or np.any(finite <= 0.0):
        return False
    return bool(finite.max() / finite.min() > 50.0)

def setup_matplotlib():
    import matplotlib.pyplot as plt

    plot_style.apply_revtex_style(plt, label_size=7.5, tick_size=6.0, legend_size=6.5)
    return plt

def figure_main(cases: list[dict[str, Any]], output_dir: Path, *, log_positive: bool = False, suffix: str = "linear") -> None:
    plt = setup_matplotlib()
    width, height = plot_style.revtex_figure_size(
        width_pt=plot_style.REVTEX_DOUBLE_COLUMN_WIDTH_PT,
        ncols=4,
        nrows=1,
        panel_aspect=1.03,
        min_height_in=1.45,
    )
    fig, axes = plt.subplots(1, 4, figsize=(width, height))
    specs = [
        ("density_total", r"$\rho$", lambda value: value, 1.0),
        ("energy_density", r"$-E$", lambda value: -NSITE * value, NSITE),
        ("S_SF_K", r"$S_{\rm SF}(K)$", lambda value: value, 1.0),
        ("S_DW_K", r"$S_{\rm DW}(K)$", lambda value: value, 1.0),
    ]
    for ax, (observable, label, transform, stderr_scale) in zip(axes, specs):
        x, ed, dqmc, stderr = values(cases, observable, transform, stderr_scale)
        use_log = log_positive and observable != "energy_density" and should_log(ed, dqmc)
        ax.plot(x, ed, **plot_style.ed_line_kwargs(label="ED"))
        ax.errorbar(
            x,
            dqmc,
            yerr=stderr,
            **plot_style.bafqmc_errorbar_kwargs(markersize=5.4, mew=0.9, capsize=1.7),
        )
        ax.set_xlabel(r"$\Delta$")
        ax.set_ylabel(label)
        ax.set_xlim(float(x.min()) - 0.015, float(x.max()) + 0.015)
        if use_log:
            ax.set_yscale("log")
        limits = y_limits(ed, dqmc, log=use_log)
        if limits is not None:
            ax.set_ylim(*limits)
        plot_style.style_axis(ax)
        exponent = plot_style.compact_linear_y_tick_labels(ax)
        plot_style.add_y_scale_text(ax, exponent, x=0.0, y=1.012)
    handles: list[Any] = []
    labels: list[str] = []
    for ax in axes:
        ax_handles, ax_labels = ax.get_legend_handles_labels()
        for handle, label in zip(ax_handles, ax_labels):
            if label not in labels:
                handles.append(handle)
                labels.append(label)
    if handles:
        axes[1].legend(
            handles,
            labels,
            loc="upper left",
            frameon=False,
            handlelength=1.5,
            borderaxespad=0.15,
            handletextpad=0.45,
            labelspacing=0.25,
        )
    plot_style.add_panel_labels(axes, x=-0.12, y=1.0)
    fig.subplots_adjust(left=0.066, right=0.992, bottom=0.22, top=0.970, wspace=0.24)
    out_dir = output_dir
    out_dir.mkdir(parents=True, exist_ok=True)
    for ext in ("pdf", "png"):
        fig.savefig(out_dir / f"benchmark_pairing_delta.{ext}")
    plt.close(fig)
