# Working in BAFQMC

This is the public computational repository for bosonic auxiliary-field quantum
Monte Carlo (BAFQMC). Help researchers run calculations, reproduce the published
benchmarks, and develop the method. The human README gives the broad picture;
use this file and the linked technical guides for implementation work.

## Agent skills

Task-specific instructions are available in `.agents/skills/`:

| Task | Skill |
| --- | --- |
| Reproduce the paper's data and figures | [bafqmc-reproduce](.agents/skills/bafqmc-reproduce/SKILL.md) |
| Implement a different lattice or Hamiltonian | [bafqmc-new-model](.agents/skills/bafqmc-new-model/SKILL.md) |
| Prepare and run a new research campaign | [bafqmc-new-calculation](.agents/skills/bafqmc-new-calculation/SKILL.md) |
| Add and verify a physical observable | [bafqmc-add-observable](.agents/skills/bafqmc-add-observable/SKILL.md) |

Read the matching skill when the task calls for it. Any coding agent can
follow these Markdown instructions; no particular agent service is required.
`CLAUDE.md` and `.github/copilot-instructions.md` also point to this file.

## Working style

- Carry the user's authorized task through implementation and relevant checks.
  Make routine reversible choices independently; keep the user informed of
  findings and results. Do not introduce extra approval steps for ordinary
  inspection, setup, local edits, or calculations already requested.
- Establish the requested Hamiltonian, observable, numerical precision, and
  computing budget from the conversation. Ask only when missing information
  would change the physical calculation or exceed the authorized resources.
- Distinguish benchmark reproduction from a new research campaign. Preserve
  published inputs and use separate files and output directories for new work.
- Report what actually ran, its parameters, and the output paths. A smoke run
  is an installation check; a full reproduction runs all selected production
  points with their original statistics and reference settings.
- Use the paper's notation in reader-facing formulas: main coupling `U`,
  creation operators `b^+,c^+`, and the paper's pairing sign. Explain code-variable
  and phase mappings explicitly. Use fenced `math` blocks and protected inline math (`$` + backticks)
  for GitHub-rendered LaTeX; put copyable agent requests in separate fenced `text` blocks.
- Use plain scientific language. Numerical comparisons should retain every
  valid point and its uncertainty; a deviation larger than three standard
  errors is not a reason to reject a reproduction or suppress a result.

For documentation changes, update the corresponding Chinese pages in `docs/zh/`
and build both languages with `python3 scripts/build_docs.py`. Verify formulas,
links, search, and copy controls in a browser; see [the documentation workflow](docs/development.md#bilingual-documentation).

## Platform and environment

Linux is the supported computational environment; Windows users run through
WSL. A Linux or WSL agent runs commands directly. A native Windows agent uses
the PowerShell launchers in `scripts/` to invoke Linux commands in this checkout.
Do not create a second checkout merely to cross the Windows/WSL boundary.

Follow [docs/getting-started.md](docs/getting-started.md) for dependencies and
the first run. The build uses MPI Fortran and BLAS/LAPACK; Python analysis uses
NumPy/Matplotlib and ED uses QuSpin. The default compiler is `mpifort`; active
Intel MPI/MKL environments are also supported through
`src/common/compiler.mk`. Inspect `make -C src/pairing print-config`
when diagnosing build configuration.

Use fast Linux storage for scratch and raw measurements, especially in WSL.
The root production runner creates scratch on the system temporary filesystem
by default and copies completed stages into the chosen output directory.
Set one numerical-library thread for baseline comparisons; the root runner
does this automatically. An alternative ED environment can be supplied with
`--python-ed /path/to/python`.

## Repository map

| Location | Responsibility |
| --- | --- |
| `reproduce.py` | Complete paper workflow; defaults to fresh full computation |
| `benchmarks/paper/production.py` | Stage execution, provenance, and continuation |
| `benchmarks/paper/analysis.py` | Means, blocking, references, and output tables |
| `benchmarks/paper/manuscript.py` | Main/supplement selection and figure notation |
| `benchmarks/paper/plot_manuscript.py` | Main two-row figure and supplemental figure |
| `benchmarks/paper/data/` | Paper inputs and compact processed benchmark data |
| `src/number_conserving/src/` | Number-conserving Fortran solver |
| `src/pairing/src/` | Full Nambu Fortran solver with onsite pairing |
| `src/common/` | Shared BLAS/LAPACK adapter, random generator, compiler setup |
| `src/*/run_paper.py` | Solver-specific execution of custom or paper manifests |
| `src/*/benchmarks/ed/` | Exact diagonalization and reference calculations |
| `docs/solvers/*/physics.md` | Hamiltonian, HS, Green functions, and estimators |
| `examples/<solver>/` | Curated example and regression input directories |
| `tests/` | Reproduction, number-conserving, and paired test suites |
| `docs/agent-workflows.md` | Concrete task recipes for research and development |
| `docs/model-development.md` | New-model derivations, source map, and independent verification |
| `docs/observables.md` | Physical operator definitions, output names, and normalizations |

The active solvers implement a periodic triangular lattice with one site per
unit cell. The hopping is `RT=1` in each `src/calc_basic.f90`; a `t` field in a
JSON file is not a general runtime hopping control. New lattices or hopping
models require coordinated changes to the solver, ED, and observables.

Implementing those new models is a supported agent workflow. Route changes to
the graph, hopping/pairing matrices, flavor structure, or interaction operators
through [the model-development guide](docs/model-development.md) and the
`bafqmc-new-model` skill. Use `bafqmc-new-calculation` for a scan once the
requested Hamiltonian is implemented. A user who asks for a new physical model
has authorized the corresponding implementation work; continue through a
runnable example and relevant verification, resolving only missing physical
choices that materially determine the calculation.

## Physical contracts

Read [docs/algorithm.md](docs/algorithm.md) and the relevant solver physics
guide before changing a kernel or estimator.

- `U1` multiplies `(n_b+n_c)^2`; `U2` multiplies `(n_b-n_c)^2`. The main-text
  coupling is `U=U2` with `U1=0`. The supplemental scan retains both channels.
  The existing continuous HS implementation uses `U1<=0` and `U2>=0`.
- Positive `t=1` is the frustrated triangular hopping convention. The trace
  uses `H - mu*N`; `energy_density` contains physical energy per site,
  excluding `-mu*N`. The paper plots `-E = -Lx*Ly*energy_density`.
- The paired solver uses `+Delta*(b^+ c^+ + b c)`. The paper uses the negative
  pair term. They are related by `c_code=-c_paper` at the same positive
  `Delta`. The four benchmark observables are invariant; the paper's anomalous
  pair amplitude is minus the code's `pair_equal`.
- The number-conserving equal-time Green function is `<b_i b_j^+>`; its
  reversed order includes the bosonic identity term. The pairing solver uses
  the full `(b,c,b^+,c^+)` Nambu basis. Preserve its block ordering and
  determinant-square-root convention.
- For the current one-site-per-cell model, `Ns=Lx*Ly`; total density is `N/Ns`, and the paper's structure factors are
  normalized by `Ns^2`. The triangular K point is `(4*pi/3,0)`; the existing
  K estimators require commensurate sizes, with both lengths multiples of 3.
  A multisublattice model has `Ns=nsub*Lx*Ly` and requires physical intracell
  positions, bond indices, Fourier form factors, and corresponding normalization.
- For the `U1=0` main model, `mu < -3*t - abs(Delta)` is the sufficient
  convergence condition used by the benchmarks. The attractive `U1<0` scan
  uses a separate finite-occupation reference. Preserve this distinction.
- Published error bars are **standard errors of the mean (SEM)** from ten
  blocks of 10000 measurement bins, not the standard deviation of individual
  measurements. Keep blocking, warmup, seeds, Trotter step, and ED cutoffs
  explicit when defining a new campaign.
- Numerical failures, missing files, nonfinite measurements, and inconsistent
  parameter mappings must be fixed. Statistical residuals are retained as
  diagnostics, without a mandatory sigma threshold for paper reproduction.

For a new interaction, derive its HS channels, normal-ordering shifts and
scalar weights. Establish TRS/RP or conjugate-sector conditions for each
decoupled field configuration, and analyze the physical thermal-trace domain
separately. Preserve the NC conjugate-sector shortcut only when the new
factors obey it. The local updates currently assume a diagonal rank-one site
change (NC) or diagonal four-sector site change (paired); different HS support
needs a corresponding update derivation. New paired models must establish
their Nambu scalar and determinant-square-root branch. Check the new physics
with independent finite-Fock references, dense fixed-field products and
nonzero proposal ratios where relevant, alongside the existing examples.

## Computation and checks

All commands below run from the repository root in the installed environment.
For initial inspection or routine orchestration changes:

```bash
python3 reproduce.py --plan
make check
python3 reproduce.py --mode check
git diff --check
```

The compact test suite skips optional live MPI tests unless enabled. Check
the actual test summary; do not report skipped tests as executed. For changes
to execution, inputs, or dependencies, run a fresh installation check:

```bash
python3 reproduce.py --mode smoke --output /tmp/bafqmc-smoke-unique
```

Choose a new output directory. For continuation changes, also run:

```bash
BAFQMC_RUN_MPI_TESTS=1 python3 -m unittest discover -s tests/reproduction -v
```

Set `BAFQMC_PYTHON_ED=/path/to/python` for a separate QuSpin interpreter.
For kernel and estimator changes, run `make physics` for independent finite-Fock
ED references, thermodynamic identities, and live Fortran analytic checks. Follow
[docs/testing.md](docs/testing.md) for coverage and choose a small interacting
comparison relevant to the changed physics.

The default command is a **full** production calculation:

```bash
python3 reproduce.py --output benchmarks/paper/output/paper-run
```

It runs all 22 points (main 15 + supplement 7), approximately 12–24 hours on
the reference desktop, with 16 GiB RAM and 8 GiB free disk recommended. Read
[benchmarks/paper/RESOURCES.md](benchmarks/paper/RESOURCES.md). Use the full run
when the user requests reproduction; routine documentation checks do not
require production sampling. `--scope main`, `--scope supplement`, and
`--model` select subsets. Continue with the same selection, output directory,
and environment plus `--resume`; completed stages are verified before reuse.
Continuation restarts an interrupted stage from its initial inputs.

## Data and changes

- Develop new models with named examples and campaign outputs so the original
  and new calculations remain easy to run. Update inputs and processed data
  intentionally when the task calls for it; record changed physical definitions,
  calculation settings, and the reason for revised results in normal documentation.
- Keep large raw chains, executables, logs, scratch, and figure previews out
  of Git. Use ignored `runs/`, `benchmarks/paper/output/`, or solver `build/`
  directories. Small curated inputs, means, SEM, block means, and ED scalar
  references are appropriate to track for new reproducible examples.
- The executables append to fixed output filenames in their working
  directories. Every independent chain needs a fresh directory. Preserve the
  archived pairing case order when selecting cases because it determines seeds.
- Keep code and docs self-contained. Use repository-relative paths or explicit
  user-supplied paths; no dependencies on a manuscript checkout or local machine.
- Preserve unrelated work, inspect the final diff, and run checks proportionate
  to the changes. Commit only task-owned work when the user's workflow calls
  for a commit. Publication or remote communication requires user authorization.
- Contributions are under MIT. The paper is not yet assigned a public
  identifier; do not invent an arXiv link, DOI, or citation entry.
