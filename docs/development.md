# Development and verification

Read [AGENTS.md](../AGENTS.md) and [the algorithm guide](algorithm.md) before
changing physical kernels. Commands in this guide run from the repository root
in the installed Linux environment.

## Locate the change

| Task | Starting points |
| --- | --- |
| Add an equal-time observable | `src/*/src/obser_equal.f90`, corresponding ED driver, analysis parser |
| Change HS fields or local proposals | `src/*/src/fields.f90`, `operator_Hubbard.f90`, `localU.f90`, `local_sweep.f90` |
| Change stabilized propagation | `src/*/src/stabilization.f90`, `process_matrix.f90`, `multiply.f90` |
| Change the lattice or hopping | `src/*/src/lattice.f90`, `calc_basic.f90`, `non_interact.f90`, and ED geometry |
| Add or schedule a campaign | Solver `run_paper.py` and `benchmarks/campaign*.py` |
| Change paper orchestration | `benchmarks/paper/production.py`, `reproduce.py` |
| Change statistics or figure output | `benchmarks/paper/analysis.py`, `plot_manuscript.py` |
| Change linear algebra or RNG | `src/common/` and its numerical invariant checks |

Each solver's Makefile lists the active Fortran sources. In particular,
number-conserving `globalK.f90` and `global_update.f90` are retained source
files outside the active executable. Check the build list before extending a
code path.

## Checks proportionate to the change

For documentation, inspect links and commands and run `git diff --check`.
For Python or orchestration changes, start with the fast checks:

```bash
make check
```

For kernels, estimators, numerical dependencies, or execution changes, run
the scientific checks against independent references:

```bash
make physics
```

The [testing guide](testing.md) describes the Gaussian Fortran calculations,
independent finite-Fock ED, thermodynamic identities, and numerical-library
checks. Use `PYTHON_ED=/path/to/python` with Make for a separate QuSpin environment.
Changes to resume/provenance logic also need the live pipeline tests:

```bash
BAFQMC_RUN_MPI_TESTS=1 python3 -m unittest discover -s tests/reproduction -p test_production.py -v
```

If ED is installed separately, set `BAFQMC_PYTHON_ED=/path/to/python` for this
test command. For an installation check through the public reproduction entry
point, run `python3 reproduce.py --mode smoke --output /tmp/bafqmc-development-smoke`
with a fresh output directory and, when needed, `--python-ed /path/to/python`.

## Physical estimators and algorithms

Write the operator definition first: flavor sums, equal-time order,
normalization, connected or disconnected pieces, and pairing phase where
applicable. Derive the estimator using the actual Green-function convention.
Add its counterpart to ED and the analysis schema, so its numerical output
retains a physical definition from solver to plot.

Add a direct small-system comparison appropriate to the changed physics.
The existing finite-Fock tests exercise interacting ED in a fixed Hilbert space;
the Gaussian tests compare actual solver measurements with an untruncated
analytic solution, including pairing phase and observable normalization.

For changes to HS updates or interacting propagation, also run the relevant
`make -C src/<solver> benchmark-dqmc` campaign and extend it to the changed regime
when needed. The Gaussian `Nwrap` comparison covers quadratic propagation;
interacting stability and severe conditioning require a test at those parameters.
Record sampling uncertainty and convergence settings with stochastic comparisons.
Paper reproduction retains every valid result without a universal sigma gate.

## Keep results reusable

Record parameters, initial seeds, compiler and numerical libraries, sample
counts, blocking, and ED cutoffs with new comparisons. Track small input
fixtures and processed reference values needed by tests. Store raw chains,
logs, compiled files, and local previews in ignored output directories.

When updating benchmark data, keep the inputs, reference parameters, block means,
and plotted statistics consistent. Document the calculation and reason for the
change so that readers can follow it through Git history. Use `--mode check`
to reconstruct the means and SEM and compare the ED reference parameters.

See the solver-specific development guides for file formats:
[number-conserving](solvers/number_conserving/development.md) and
[pairing](solvers/pairing/development.md).

## Bilingual documentation

English pages are built from the root README, `docs/`, and solver guides;
Chinese pages live in `docs/zh/`. Keep corresponding physical definitions,
commands, and resource estimates aligned when editing either language. Both
configurations are built and deployed together, with Chinese at `/BAFQMC/zh/`.
Use repository-relative Markdown links and images; the builder resolves them
for the website while keeping them usable on GitHub.

```bash
python3 -m pip install -r docs/requirements.txt
python3 scripts/build_docs.py
mkdir -p .build/preview
ln -sfn ../site .build/preview/BAFQMC
python3 -m http.server 8000 --directory .build/preview
```

Open `http://localhost:8000/BAFQMC/` or `http://localhost:8000/BAFQMC/zh/`.
Check formulas, language links, Chinese search, and code-copy buttons in a
browser. A successful strict build also checks local Markdown links.
