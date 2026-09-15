# Plotting Style

Implementation paths and commands in this guide are relative to
`src/pairing/` from the repository root.

This project uses a compact REVTEX paper style for publication figures.  The
reusable Matplotlib helpers live in `benchmarks/plot_style.py`.

## REVTEX Size

The manuscript uses a REVTEX document class with a single-column line width

```text
\linewidth = 246.0 pt
```

A single-column figure should therefore use width `246.0 pt`, or
`246.0 / 72.27 = 3.405 in`.  A two-column paper figure should use
`492.0 pt`, or `2 * 246.0 / 72.27 = 6.811 in`.  The 3x3 benchmark candidate
figures are four panels per row and should normally use the two-column width,
so each data frame has enough room for axis labels and tick labels.

Use `plot_style.revtex_figure_size(...)` to derive the Matplotlib figure size:
pass `width_pt=plot_style.REVTEX_LINEWIDTH_PT` for single-column figures and
`width_pt=plot_style.REVTEX_DOUBLE_COLUMN_WIDTH_PT` for two-column figures.

## Matplotlib Rules

- Font: DejaVu Sans.
- Base font size: 10 pt, matching the manuscript body text.
- Smaller tick and legend labels are allowed, but never below 6 pt.
- Line width: 0.75 pt, matching the recent REVTEX benchmark-figure style.
- ED/reference curves: point-line plot with small filled markers by default.
- BAFQMC data: transparent open markers by default, slightly larger than ED
  markers, so nearly overlapping benchmark points remain distinguishable
  without covering the ED marker underneath.  Do not use white-filled markers
  for this overlay, because they hide the reference point.
- Axes: four-sided frame, with square panel boxes for paper figures.
- Ticks: inward on the left and bottom axes.  Keep the top and right frame
  lines, but do not draw tick marks on the top or right frame lines.
- Y-axis tick labels in compact multi-panel figures may be rotated by 90
  degrees; when rotated, center-anchor them on the tick mark rather than
  right-aligning them.
- Scientific tick labels: use math notation such as
  `$4.5\times10^{-3}$`, not `4.5e-3`.
- For compact linear axes where every tick would need a long scientific
  label, factor the common power of ten into a small y-axis offset text near
  the upper-left of the frame, for example `\times10^{-4}`, while keeping the
  y-axis label as the physical quantity and tick labels as short mantissas
  such as `4.5`, `6`, and `7.5`.
- Grid: off.
- Panel titles: omit titles inside paper figures.
- Panel labels: bold lowercase `a,b,c,d,...`, placed at the upper-left frame
  without adding extra top whitespace.  The default helper aligns labels to
  the top edge (`va="top"`), similar to the reference paper style.  Labels
  should not float above the row and should not overlap ticks or axis labels.
- Final paper output: vector PDF, with `dpi=300` used for raster previews.

## Example

```python
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt

from benchmarks import plot_style

plot_style.apply_revtex_style(plt)
fig, axes = plt.subplots(
    2,
    4,
    figsize=plot_style.revtex_figure_size(
        width_pt=plot_style.REVTEX_DOUBLE_COLUMN_WIDTH_PT,
        ncols=4,
        nrows=2,
    ),
    constrained_layout=True,
)

for ax in axes.flat:
    plot_style.style_axis(ax)
    ax.plot(x_ed, y_ed, **plot_style.ed_line_kwargs())
    ax.errorbar(x_dqmc, y_dqmc, yerr=err, **plot_style.bafqmc_errorbar_kwargs())

plot_style.add_panel_labels(axes)
fig.savefig("figure.pdf")
fig.savefig("figure.png", dpi=300)
```

For one-off diagnostic plots it is acceptable to use larger figures, but any
figure proposed for the paper should follow this style.
