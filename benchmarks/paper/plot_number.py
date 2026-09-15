"""Original number-conserving paper layout, extracted without style changes.
Source: code_bosonDQMC campaign make_publication_report.py; hashes in data/index.json.
"""
from __future__ import annotations
import math
from pathlib import Path
from typing import Any
from . import plot_style

def _plot_candidate_grid(
    plt: Any,
    figure_dir: Path,
    row_specs: list[tuple[str, list[dict[str, Any]]]],
    observables: tuple[str, ...],
    stem: str,
    *,
    panel_aspect: float,
    plot_checkpoints: bool,
) -> dict[str, str] | None:
    if not any(records for _, records in row_specs):
        return None
    old_params = dict(plt.rcParams)
    plot_style.apply_revtex_style(plt, label_size=7.5, tick_size=6.0, legend_size=6.0)
    fig, axes = plt.subplots(
        len(row_specs),
        len(observables),
        figsize=plot_style.revtex_figure_size(
            width_pt=plot_style.REVTEX_DOUBLE_COLUMN_WIDTH_PT,
            ncols=len(observables),
            nrows=len(row_specs),
            panel_aspect=panel_aspect,
        ),
        squeeze=False,
    )
    fig.subplots_adjust(
        left=0.066,
        right=0.992,
        bottom=0.112,
        top=0.970,
        wspace=0.24,
        hspace=0.12,
    )
    any_plotted = False
    for row, (row_label, row_records) in enumerate(row_specs):
        for col, observable in enumerate(observables):
            ax = axes[row][col]
            plotted = _plot_observable_axis(
                ax,
                row_records,
                observable,
                plot_checkpoints=plot_checkpoints,
            )
            plot_style.style_axis(ax)
            any_plotted = any_plotted or plotted
            scale_exponent = plot_style.compact_linear_y_tick_labels(ax)
            ax.set_ylabel(_observable_label(observable))
            plot_style.add_y_scale_text(ax, scale_exponent)
            ax.set_xlabel(_x_label_for_records(row_records))
            if not plotted:
                ax.text(
                    0.5,
                    0.5,
                    "no data",
                    ha="center",
                    va="center",
                    transform=ax.transAxes,
                    color="#667085",
                )
    if not any_plotted:
        plt.close(fig)
        plt.rcParams.update(old_params)
        return None
    handles: list[Any] = []
    labels: list[str] = []
    for ax in axes.flat:
        ax_handles, ax_labels = ax.get_legend_handles_labels()
        for handle, label in zip(ax_handles, ax_labels):
            if label not in labels:
                handles.append(handle)
                labels.append(label)
    if handles:
        axes[0][1].legend(
            handles,
            labels,
            loc="upper right",
            frameon=False,
            handlelength=1.5,
            borderaxespad=0.15,
            handletextpad=0.45,
            labelspacing=0.25,
        )
    plot_style.add_panel_labels(axes, x=-0.12, y=1.0, fontsize=8.0)
    png_path = figure_dir / f"{stem}.png"
    pdf_path = figure_dir / f"{stem}.pdf"
    fig.savefig(png_path, dpi=300)
    fig.savefig(pdf_path)
    plt.close(fig)
    plt.rcParams.update(old_params)
    return {
        "relative_path": f"figures/{stem}.png",
        "sweep": stem,
        "observable": "combined",
    }

def _plot_observable_axis(
    ax: Any,
    records: list[dict[str, Any]],
    observable: str,
    *,
    plot_checkpoints: bool = True,
) -> bool:
    dqmc_points: list[tuple[float, float, float]] = []
    ed_points: list[tuple[float, float]] = []
    ed_checkpoint_points: list[tuple[float, float]] = []
    for record in records:
        values = record["observables"].get(observable, {})
        dqmc = _finite_or_none(values.get("dqmc"))
        ed = _finite_or_none(values.get("ed"))
        stderr = _finite_or_none(values.get("stderr"))
        if dqmc is not None:
            dqmc_points.append((record["x"], dqmc, 0.0 if stderr is None else stderr))
        if ed is not None:
            target = ed_points if bool(record.get("ed_reliable", False)) else ed_checkpoint_points
            target.append((record["x"], ed))
    if not dqmc_points and not ed_points and not ed_checkpoint_points:
        return False

    dqmc_points.sort(key=lambda item: item[0])
    ed_points.sort(key=lambda item: item[0])
    ed_checkpoint_points.sort(key=lambda item: item[0])
    plotted_values: list[float] = []
    if ed_points:
        xs = [item[0] for item in ed_points]
        ys = [item[1] for item in ed_points]
        plotted_values.extend(ys)
        ax.plot(xs, ys, **plot_style.ed_line_kwargs(label="ED"))
    if ed_checkpoint_points:
        xs = [item[0] for item in ed_checkpoint_points]
        ys = [item[1] for item in ed_checkpoint_points]
        plotted_values.extend(ys)
        if plot_checkpoints:
            ax.plot(xs, ys, **plot_style.checkpoint_kwargs(label="ED checkpoint"))
    if dqmc_points:
        xs = [item[0] for item in dqmc_points]
        ys = [item[1] for item in dqmc_points]
        err = [max(0.0, item[2]) for item in dqmc_points]
        plotted_values.extend(ys)
        ax.errorbar(
            xs,
            ys,
            yerr=err,
            **plot_style.bafqmc_errorbar_kwargs(
                label="BAFQMC",
                markersize=5.4,
                mew=0.9,
                capsize=1.7,
            ),
        )
    if _should_log_y(plotted_values, observable):
        ax.set_yscale("log")
    ax.margins(y=0.08)
    return True

def _x_label_for_records(records: list[dict[str, Any]]) -> str:
    return "$U_2$" if any(record.get("vary") == "U2" for record in records) else "$U_1$"

def _should_log_y(values: list[float], observable: str) -> bool:
    if observable in {"energy_density", "IPR"} or observable.endswith("_over_density"):
        return False
    positive = [value for value in values if math.isfinite(value) and value > 0.0]
    if len(positive) < 2 or len(positive) != len([v for v in values if math.isfinite(v)]):
        return False
    return max(positive) / min(positive) > 50.0 or observable.endswith("_over_density")

def _observable_label(name: str) -> str:
    labels = {
        "density_total": r"$\rho$",
        "energy_density": r"$e$",
        "total_energy": r"$E$",
        "minus_total_energy": r"$-E$",
        "doubleOcc": r"$D$",
        "squareOcc": r"$P_{\rm same}$",
        "num_up": r"$N_b$",
        "num_do": r"$N_c$",
        "numsquare_up": r"$N_b^2$",
        "numsquare_do": r"$N_c^2$",
        "onsite_n2_up": r"$\langle n_b^2\rangle$",
        "onsite_n2_do": r"$\langle n_c^2\rangle$",
        "IPR": r"$\mathrm{IPR}_{\rho}$",
        "S_SF_K": r"$S_{\rm SF}(K)$",
        "S_PSF_Gamma": r"$S_{\rm PSF}(\Gamma)$",
        "S_DW_K": r"$S_{\rm DW}(K)$",
        "S_SF_K_over_density": r"$S_{\rm SF}(K)/\rho$",
        "S_DW_K_over_density": r"$S_{\rm DW}(K)/\rho$",
    }
    return labels.get(name, name)

def _finite_or_none(value: Any) -> float | None:
    try:
        number = float(value)
    except (TypeError, ValueError):
        return None
    return number if math.isfinite(number) else None
