"""Matplotlib style helpers for REVTEX-width benchmark figures."""

from __future__ import annotations

import math
import string
from typing import Any


PT_PER_INCH = 72.27
REVTEX_LINEWIDTH_PT = 246.0
REVTEX_DOUBLE_COLUMN_WIDTH_PT = 2.0 * REVTEX_LINEWIDTH_PT
DEFAULT_FONT = "DejaVu Sans"
DEFAULT_FONT_SIZE_PT = 10.0
MIN_FONT_SIZE_PT = 6.0
LINE_WIDTH_PT = 0.75
ED_MARKER_SIZE_PT = 2.4
BAFQMC_MARKER_SIZE_PT = 4.0


def pt_to_inch(value_pt: float) -> float:
    """Convert TeX points to inches."""

    return float(value_pt) / PT_PER_INCH


def revtex_figure_size(
    *,
    width_pt: float = REVTEX_LINEWIDTH_PT,
    ncols: int = 4,
    nrows: int = 1,
    panel_aspect: float = 1.05,
    min_height_in: float = 1.1,
) -> tuple[float, float]:
    """Return a REVTEX-width figure size.

    The width is fixed by the manuscript line width.  The height is derived
    from the per-panel width, because a common paper layout packs four panels
    into one manuscript line including labels and panel margins.
    """

    if ncols <= 0 or nrows <= 0:
        raise ValueError("ncols and nrows must be positive")
    width_in = pt_to_inch(width_pt)
    panel_width_in = width_in / float(ncols)
    height_in = max(float(min_height_in), panel_width_in * float(panel_aspect) * nrows)
    return (width_in, height_in)


def rcparams(
    *,
    font_size: float = DEFAULT_FONT_SIZE_PT,
    label_size: float = 7.5,
    tick_size: float = 6.0,
    legend_size: float = 6.5,
) -> dict[str, Any]:
    """Return Matplotlib rcParams for publication figures."""

    for name, value in {
        "font_size": font_size,
        "label_size": label_size,
        "tick_size": tick_size,
        "legend_size": legend_size,
    }.items():
        if float(value) < MIN_FONT_SIZE_PT:
            raise ValueError(f"{name} must be at least {MIN_FONT_SIZE_PT:g} pt")
    return {
        "font.family": "sans-serif",
        "font.sans-serif": [DEFAULT_FONT],
        "font.size": float(font_size),
        "axes.labelsize": float(label_size),
        "axes.titlesize": float(label_size),
        "axes.linewidth": LINE_WIDTH_PT,
        "axes.grid": False,
        "axes.spines.left": True,
        "axes.spines.right": True,
        "axes.spines.top": True,
        "axes.spines.bottom": True,
        "xtick.direction": "in",
        "ytick.direction": "in",
        "xtick.top": False,
        "ytick.right": False,
        "xtick.labelsize": float(tick_size),
        "ytick.labelsize": float(tick_size),
        "xtick.major.width": LINE_WIDTH_PT,
        "ytick.major.width": LINE_WIDTH_PT,
        "xtick.minor.width": LINE_WIDTH_PT,
        "ytick.minor.width": LINE_WIDTH_PT,
        "xtick.major.pad": 1.0,
        "ytick.major.pad": 4.0,
        "axes.labelpad": 3.0,
        "lines.linewidth": LINE_WIDTH_PT,
        "lines.markersize": ED_MARKER_SIZE_PT,
        "legend.fontsize": float(legend_size),
        "figure.dpi": 300,
        "savefig.dpi": 300,
        "savefig.bbox": "tight",
        "savefig.pad_inches": 0.015,
        "pdf.fonttype": 42,
        "ps.fonttype": 42,
    }


def apply_revtex_style(plt: Any, **kwargs: Any) -> None:
    """Apply the project REVTEX paper style to a Matplotlib pyplot module."""

    plt.rcParams.update(rcparams(**kwargs))


def style_axis(ax: Any) -> None:
    """Apply axis-level style that rcParams cannot fully enforce."""

    from matplotlib.ticker import FuncFormatter, LogFormatterMathtext, LogLocator
    from matplotlib.ticker import MaxNLocator, NullFormatter

    ax.grid(False)
    ax.set_box_aspect(1.0)
    for spine in ax.spines.values():
        spine.set_visible(True)
        spine.set_linewidth(LINE_WIDTH_PT)
    ax.xaxis.set_ticks_position("bottom")
    ax.yaxis.set_ticks_position("left")
    ax.tick_params(axis="x", which="both", direction="in", top=False, bottom=True)
    ax.tick_params(
        axis="y",
        which="both",
        direction="in",
        right=False,
        left=True,
        pad=4.0,
    )
    if ax.get_yscale() == "log":
        ax.yaxis.set_major_locator(LogLocator(base=10.0, subs=(1.0,)))
        ax.yaxis.set_major_formatter(LogFormatterMathtext(base=10.0))
        ax.yaxis.set_minor_locator(LogLocator(base=10.0, subs=tuple(range(2, 10))))
        ax.yaxis.set_minor_formatter(NullFormatter())
    else:
        ax.yaxis.set_major_locator(MaxNLocator(nbins=3))
        ax.yaxis.set_major_formatter(FuncFormatter(_compact_numeric_tick))
    for tick in ax.get_yticklabels():
        tick.set_rotation(90)
        tick.set_rotation_mode("anchor")
        tick.set_va("center")
        tick.set_ha("center")


def compact_linear_y_tick_labels(
    ax: Any,
    *,
    small_threshold: float = 1.0e-2,
    large_threshold: float = 1.0e3,
) -> int | None:
    """Factor a common power of ten out of compact linear y tick labels.

    Returns the extracted base-10 exponent, or ``None`` when no scaling was
    applied.  Use :func:`add_y_scale_text` to display the returned exponent as
    a small axis offset text.
    """

    if ax.get_yscale() != "linear":
        return None

    from matplotlib.ticker import FuncFormatter

    ymin, ymax = ax.get_ylim()
    finite = [abs(float(value)) for value in (ymin, ymax) if math.isfinite(float(value))]
    if not finite:
        return None
    max_abs = max(finite)
    if max_abs == 0.0 or small_threshold <= max_abs < large_threshold:
        return None

    exponent = int(math.floor(math.log10(max_abs)))
    scale = 10.0**exponent
    ax.yaxis.set_major_formatter(
        FuncFormatter(lambda value, _position=None: _plain_numeric_tick(value / scale))
    )
    return exponent


def add_y_scale_text(
    ax: Any,
    exponent: int | None,
    *,
    x: float = 0.0,
    y: float = 1.015,
    fontsize: float = MIN_FONT_SIZE_PT,
) -> Any | None:
    """Draw a Matplotlib-style common y-axis multiplier near the y ticks."""

    if exponent is None:
        return None
    return ax.text(
        x,
        y,
        rf"$\times10^{{{exponent}}}$",
        transform=ax.transAxes,
        ha="left",
        va="bottom",
        fontsize=float(fontsize),
        clip_on=False,
    )


def _compact_numeric_tick(value: float, _position: int | None = None) -> str:
    """Return short numeric tick labels for narrow REVTEX panels."""

    try:
        numeric = float(value)
    except (TypeError, ValueError):
        return ""
    if not numeric:
        return "0"
    magnitude = abs(numeric)
    if magnitude < 1.0e-2 or magnitude >= 1.0e3:
        mantissa, exponent = f"{magnitude:.1e}".split("e")
        mantissa = mantissa.rstrip("0").rstrip(".")
        exponent_int = int(exponent)
        sign = "-" if numeric < 0.0 else ""
        return rf"${sign}{mantissa}\times10^{{{exponent_int}}}$"
    return f"{numeric:.3g}"


def _plain_numeric_tick(value: float, _position: int | None = None) -> str:
    """Return a compact non-scientific numeric tick label."""

    try:
        numeric = float(value)
    except (TypeError, ValueError):
        return ""
    if abs(numeric) < 1.0e-12:
        return "0"
    return f"{numeric:.3g}"


def add_panel_labels(
    axes: Any,
    *,
    labels: list[str] | None = None,
    x: float = -0.12,
    y: float = 1.0,
    fontsize: float = 8.0,
    ha: str = "right",
    va: str = "top",
) -> None:
    """Add bold lowercase panel labels at the upper-left frame."""

    if float(fontsize) < MIN_FONT_SIZE_PT:
        raise ValueError(f"panel label fontsize must be at least {MIN_FONT_SIZE_PT:g} pt")
    flat_axes = list(getattr(axes, "flat", axes))
    default_labels = list(string.ascii_lowercase)
    panel_labels = labels if labels is not None else default_labels
    if len(panel_labels) < len(flat_axes):
        raise ValueError("not enough panel labels for axes")
    for ax, label in zip(flat_axes, panel_labels):
        ax.text(
            x,
            y,
            label,
            transform=ax.transAxes,
            ha=ha,
            va=va,
            fontsize=float(fontsize),
            fontweight="bold",
            clip_on=False,
        )


def ed_line_kwargs(**overrides: Any) -> dict[str, Any]:
    """Keyword defaults for ED curves in paper comparison figures."""

    kwargs: dict[str, Any] = {
        "marker": "o",
        "markersize": ED_MARKER_SIZE_PT,
        "linewidth": LINE_WIDTH_PT,
        "color": "#c23b22",
        "markerfacecolor": "#c23b22",
        "markeredgecolor": "#c23b22",
        "label": "ED / exact",
        "zorder": 3,
    }
    kwargs.update(overrides)
    return kwargs


def bafqmc_errorbar_kwargs(**overrides: Any) -> dict[str, Any]:
    """Keyword defaults for BAFQMC open-marker data in paper figures."""

    kwargs: dict[str, Any] = {
        "fmt": "o",
        "markersize": BAFQMC_MARKER_SIZE_PT,
        "markerfacecolor": "none",
        "mec": "#1f6f8b",
        "mew": LINE_WIDTH_PT,
        "elinewidth": LINE_WIDTH_PT,
        "capsize": 1.6,
        "fillstyle": "none",
        "linestyle": "none",
        "color": "#1f6f8b",
        "label": "BAFQMC",
        "zorder": 4,
    }
    kwargs.update(overrides)
    return kwargs


def checkpoint_kwargs(**overrides: Any) -> dict[str, Any]:
    """Keyword defaults for incomplete ED checkpoints."""

    kwargs: dict[str, Any] = {
        "linestyle": "none",
        "marker": "D",
        "markersize": ED_MARKER_SIZE_PT,
        "markerfacecolor": "white",
        "markeredgecolor": "#8a6f2a",
        "markeredgewidth": LINE_WIDTH_PT,
        "color": "#8a6f2a",
        "label": "ED checkpoint",
        "zorder": 2,
    }
    kwargs.update(overrides)
    return kwargs
