# Development Guide

Commands and component paths below are relative to `src/number_conserving/`
from the repository root. Shared documentation lives in `docs/solvers/number_conserving/`.

## Layout

- `src/` contains Fortran source files.
- `build/` contains generated objects, modules, and `bosonDQMC.out`.
- `../../examples/number_conserving/examples/` contains committed run input directories.
- `scripts/` contains local run helpers.
- `benchmarks/` contains comparison scripts, fixture outputs, ED references, and optional ED tooling.
- `../../docs/solvers/number_conserving/` contains physics and development documentation.

Publication figure conventions are documented in [the plotting guide](plotting.md), with
Matplotlib helpers in `benchmarks/plot_style.py`.

## Build Targets

```bash
make print-config
make build
make run-example
make benchmark
make benchmark-fast
make benchmark-dqmc
make check-fixtures
make benchmark-ed
make clean
```

`make build` compiles the active object list in `Makefile`. `globalK.f90` and `global_update.f90` are kept in `src/` but are not part of the active executable.

## Run Directory Contract

The executable reads these files from the current working directory:

- `paramC_sets.txt`
- `confin.txt`
- `seeds.txt`

The executable writes scalar observables and logs into the current working directory. Run from a dedicated directory to avoid mixing outputs from unrelated parameter sets.

Output filenames are fixed by the Fortran executable. Keep generated outputs inside run directories unless the runtime I/O contract is being deliberately changed.

## Adding A Run Case

Create a directory under `../../examples/number_conserving/examples/` or an untracked working run directory with:

- `paramC_sets.txt`
- `confin.txt`
- `seeds.txt`

Run it locally:

```bash
bash scripts/run_local.sh path/to/run_dir 1
```

The default example target copies the distributed inputs to `build/example` and writes scalar observables, logs, and continuous pole diagnostics there:

```bash
make run-example
```

Validate the pole diagnostic file structure and cross-file consistency for any completed run with:

```bash
python3 benchmarks/check_pole_diagnostics.py <run_dir>
```

## Adding A Benchmark Case

Create:

- `benchmarks/references/<case>.json`
- fixture outputs under `benchmarks/fixtures/<case>_mc_outputs/`
- optional ED params under `benchmarks/ed/params_<case>.txt`

The reference JSON must include `dqmc_fixture` and document whether each expected value is total, per flavor, per site, or divided by `Lq`. After adding a case, run:

```bash
python3 -m unittest discover -s ../../tests/number_conserving -p 'test_compare.py' -v
make check-fixtures
```

For changes that substantively alter the DQMC algorithm, run the live benchmark:

```bash
make benchmark
```

`make benchmark` is an alias for `make benchmark-dqmc`. It uses fresh temporary run directories, so generated output from older runs cannot contaminate the comparison. The default suite is listed in `benchmarks/dqmc_suite.json` and includes one free analytic case plus all four ED reference cases preserved in `benchmarks/dqmc_references/`.

All live suite cases use `dtau = beta / Ltrot = 0.01`. Interacting live cases use `Nbin = 100000` and are compared to ED with block-estimated standard errors of the Monte Carlo mean. On the validation workstation under Linux/WSL with `MPI_NP=1`, the full suite was observed at `real 604.50` seconds, about 10 minutes 5 seconds; budget at least 15 minutes and do not shorten `Ltrot` or `Nbin` unless the benchmark definition is being deliberately changed.

## Generated Files

Generated build files are under `build/` and removed by:

```bash
make clean
```

Generated research outputs under the repository-root `runs/` directory are ignored by git. Fixture outputs under `benchmarks/fixtures/` are committed when they are intentionally used by `make check-fixtures`.
