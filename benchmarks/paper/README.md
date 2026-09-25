# Benchmark reproduction

The default `python3 reproduce.py` reruns all 22 BAFQMC/reference calculations,
processes the new measurements, and draws the main-text and supplemental
benchmark figures of the accompanying study.
See [RESOURCES.md](RESOURCES.md) for the 12–24 hour production budget, memory,
scratch storage, and continuation commands. The optional `--mode plot` redraws
the small stored dataset. Both paths draw the combined main-text benchmark and
the separate attractive total-density benchmark in the Supplemental Material.

## Data and reference conventions

| Figure / row | Varying parameter | Fixed parameters | Reference |
|---|---|---|---|
| Main Fig. 2(a–d) | U1 = 0, 0.25, 0.5, 0.75, 1, 1.25, 1.5, 2 | U2 = 0, mu = -3.5, beta = 4 | Particle-shell ED; analytic free-boson reference at U1 = 0 |
| Main Fig. 2(e–h) | Delta = 0, 0.05, 0.1, 0.15, 0.2, 0.25, 0.3 | U1 = 1, U2 = 0, mu = -5, beta = 4 | Dense ED with nmax = 3 and ncut = 4 |
| Supplemental Fig. S1(a–d) | U2 = -0.6, -0.5, -0.4, -0.3, -0.2, -0.1, 0 | U1 = 1, mu = -7, beta = 1 | Low-density finite-occupation-window reference for U2 < 0 |

All inputs, ED records, processed tables, and plot metadata use the manuscript
notation directly: `U1` is the repulsive relative-density coefficient and `U2`
is the attractive total-density coefficient.


`--scope main` selects the 15 points in the main-text figure, and
`--scope supplement` selects the seven attractive-density points. The default
`--scope all` includes both. Scope and solver selection are independent:

| Selection | Points | Scan |
|---|---|---|
| `--scope main` | 15 | U1 and Delta |
| `--scope main --model number_conserving` | 8 | U1 |
| `--scope main --model pairing` | 7 | Delta |
| `--scope supplement` | 7 | U2 |
| `--model number_conserving` | 15 | U1 and U2 |

The main-text model satisfies `mu < -3*t - abs(Delta)` at every plotted point.
As proved in the study's supplemental subsection "Convergence throughout the auxiliary-field
domain", its quadratic propagation and unitary relative-density HS factors
give a finite trace throughout the auxiliary-field domain. This result applies
to the main scan with `U2=0`; the supplemental scan uses the separate
two-channel model and its stated finite-occupation references.

The pairing implementation writes the real pair term with a positive input
coefficient, whereas the manuscript writes it as `-Delta (b^+ c^+ + b c)`.
The phase convention `c_code = -c_paper` relates these expressions at the same
positive scan value. Density, total physical energy, and the two plotted
structure factors are invariant under this transformation. The solver's
anomalous pair amplitude `pair_equal` changes sign when expressed in manuscript
operators; see [the pairing conventions](../../docs/solvers/pairing/physics.md).

All points have t = 1, a periodic 3 x 3 triangular lattice, Delta tau = 0.01,
100000 measurement bins, and one MPI rank. The error bars are standard errors
from ten consecutive blocks of 10000 bins. Warmup and stabilization settings,
proposal widths, and exact initial seeds are in each `inputs/<case>/` directory.
The paired `seeds.txt` contains a list of rank seeds; the paper uses its first
entry. It is not a multiword state for one random-number stream.

The four measurement files are `density_total`, `energy_density`, `sf_K`, and
`dw_K`. Structure-factor files store real and imaginary components. The plots
use their real parts and show the physical total energy as -E = -9 e, with the
energy standard error multiplied by 9. The chemical-potential term is excluded
from E. No statistical acceptance criterion removes points during reproduction.

The supplemental attractive-U2 reference uses the low-density occupation
window: six completed particle shells, with a configured maximum of eight.
This is the finite-window comparison encoded by `low_density_cutoff_accepted`
in the ED records. The U1 = 0 reference in the main-text interaction scan is
computed analytically; all other reference values are read from the stored ED
results.

## Layout

```text
data/index.json                  selected cases, parameters, plotted values, and metadata
data/observables.csv             compact processed means, SEM, and ED values
data/block_means.csv             ten block means per case and observable
data/<model>/manifest.json       executable paper campaign
data/<model>/inputs/<case>/      exact initial input files and ED parameters
data/<model>/ed/<case>.json      compact ED reference results
data/<model>/raw/<case>.tar.gz   optional local raw chains, ignored by Git
analysis.py                     raw-data statistics, references, tables, plot dispatch
plot_manuscript.py               main U1/Delta grid and supplemental U2 scan
plot_number.py                  2 x 4 number-conserving layout
plot_pairing.py                 1 x 4 paired layout
plot_style.py                   shared figure style
```

The tracked data package is approximately 250 kB. Large raw measurement files
are generated by BAFQMC and excluded from Git. The small ED JSON files contain
reference expectation values and cutoff metadata. Stored block means permit
independent reconstruction of the plotted mean and standard error without
shipping the full Monte Carlo chains. `--mode check` checks finite records,
ED parameters, and reconstruction of the plotted means and SEM without making
figures. `--mode raw` reconstructs the statistics from optional local raw
chains when they are present locally.

The current outputs are `figures/benchmark_combined.pdf` (two rows) and
`figures/benchmark_attractive.pdf` (one row), with PNG previews alongside them.
Only the figures and rows covered by the selected scope and solver are drawn.
An individual main-text row is written as
`figures/benchmark_combined_number_conserving.pdf` or
`figures/benchmark_combined_pairing.pdf`, preserving panel labels (a–d) or
(e–h), respectively.
The layout modules remain available as reusable plotting helpers.

Generated tables include `benchmark`, `section`, `figure`, `figure_row`, and
`panels` to locate each point in the current manuscript. `scan_parameter` is
the input key and plotted symbol. The JSON records collect these layout fields
under `manuscript`. Case IDs follow the same `U1`/`U2` notation; `row`
identifies the benchmark family and `figure_row` identifies the current figure
arrangement.

## Fresh calculations

`--mode full` (the default) builds each solver, executes the paper campaign, runs ED,
reblocks the new measurements, and draws both figures from the new results.
The tracked processed data stay fixed in `data/`. Output goes to the directory selected
by `--output` and includes the numerical differences from ED.

Fresh Monte Carlo trajectories can differ with compiler and numerical library:
the tracked data remain the exact source for regenerating the paper curves.
PDF metadata and font-library versions can also change PDF byte hashes even
when the plotted data and layout are identical.

For separate stages, `--mode dqmc` and `--mode ed` write into `output/runs/`.
After both stages, `--mode analyze` reblocks their outputs and draws the figures.
Use the same `--scope` and `--model` for each stage and for `--resume`.
Each solver's `run_paper.py --help` also exposes case selection or campaign
configuration and an analysis entry point. Copy the supplied inputs into a
separate configuration directory when creating a new model campaign. All source
files, campaign inputs, processed data, and reference implementations required
for reproduction are included here. Production generates fresh raw measurements
directly from the included inputs.
