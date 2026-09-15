# Development Guide

Commands and component paths below are relative to `src/pairing/`
from the repository root. Shared documentation lives in `docs/solvers/pairing/`.

## Layout

- `src/` contains the active Fortran source files.
- `build/` contains generated objects, modules, and `build/bosonDQMC.out`.
- `scripts/` contains local and benchmark run helpers.
- `../../examples/pairing/examples/triangle_pairing_3x2/` contains the committed example input files.
- `benchmarks/` contains comparison scripts, checked-in fixtures, reference JSON, and QuSpin ED scripts.
- `../../docs/solvers/pairing/` contains physics and development documentation. Keep generated research outputs under the repository-root `runs/` directory.

## Runtime Input Contract

The executable reads fixed filenames from the current working directory:

- `paramC_sets.txt`
- `confin.txt`
- `seeds.txt`

`paramC_sets.txt` currently reads these rows:

```text
RU1 RU2 mu RDelta
Nlx Nly Ltrot Beta
NlxTherm NlyTherm LtrotTherm
Nwrap Nbin Nsweep shiftLoc
is_tau Nthermal
is_warm Nwarm shiftWarm1 shiftWarm2
iniType iniAmpl iniBias1 iniBias2
```

Run from a dedicated directory. The program writes scalar observables, correlation files, `info.txt`, and `confout.txt` into that same directory.

## Generated Files

Generated build outputs are ignored:

```text
build/
*.o
*.mod
*.out
*.lst
*.opt-report
```

Generated research outputs under the repository-root `runs/` directory are ignored, including:

```text
density density_up density_do num_up num_do kinetic
doubleOcc squareOcc local_numsquare numsquare_up numsquare_do
density_total density_site_total onsite_n2_up onsite_n2_do
pair_equal interaction_energy_density pairing_energy_density
chemical_energy_density grand_energy_density energy_density
sf_K dw_K psf_Gamma den_upup_sub* den_dodo_sub* den_updo
confout.txt info.txt calc.log
```

Fixture outputs under `benchmarks/fixtures/` are committed only when used by `make check-fixtures`.

## Makefile Targets

Current targets:

```bash
make print-config
make build
make run-example
make check-fixtures
make benchmark-ed
make benchmark-fast
make benchmark-dqmc
make benchmark
make clean
```

`make benchmark` aliases `make benchmark-dqmc`. The live benchmark manifest currently uses:

```text
../../examples/pairing/examples/triangle_pairing_ed_smoke
benchmarks/references/triangle_pairing_ed_smoke.json
../../examples/pairing/benchmarks/triangle_pairing_delta0_nopairing_3x2
benchmarks/references/triangle_pairing_delta0_nopairing.json
```

These are development regression inputs. For new research, keep the complete campaign parameters, seeds, analysis settings, and results together under the ignored repository-root `runs/` directory.

## Local Runs

Local run:

```bash
make run-example RUN_DIR=build/local-example MPI_NP=1
```

The ED driver estimates the constrained basis and dense memory requirement
before building the QuSpin object. Choose the cutoff to match allocated memory;
for a 3x3 system:

```text
nmax ncut  basis.Ns  dense H.eigh RSS estimate
3    4      7297      3.2 GiB
3    5     33307     66.1 GiB
5    5     33649     67.5 GiB
```

A full-trace calculation at `ncut=5` requires a node with sufficient memory.
Sparse low-energy calculations are a separate trace approximation and should
record their included states and temperature range.

## Staged Verification

For documentation-only changes:

```bash
git diff --check
```

For Python or benchmark-comparison changes:

```bash
python3 -m py_compile benchmarks/compare.py benchmarks/run_dqmc_suite.py
python3 -m py_compile benchmarks/pairing_delta_analysis.py benchmarks/campaign_pairing_delta.py
python3 -m unittest discover -s ../../tests/pairing -p 'test_pairing*.py' -v
make check-fixtures
```

For ED changes:

```bash
python benchmarks/ed/compare_ed_implementations.py --case tiny
python benchmarks/ed/test_observable_identities.py --case tiny
```

For Fortran or runtime changes:

```bash
make print-config
make build
make run-example
```

For physics-level changes to observable normalization, Nambu Wick formulas, or `Delta=0` regressions, update the benchmark references and run the applicable live benchmark suite with fresh measurements.

For Delta-sweep research campaigns, initialize ignored run directories with:

```bash
python3 benchmarks/campaign_pairing_delta.py write-manifest \
  --output ../../runs/<dated-campaign>/manifest.json
python3 benchmarks/campaign_pairing_delta.py init-local \
  --manifest ../../runs/<dated-campaign>/manifest.json \
  --output-dir ../../runs/<dated-campaign>/inputs \
  --seed 24681357
```

Run ED in the repository QuSpin environment and keep parallelism within the
allocated memory. The portable `run_paper.py` runner executes ED cases serially.
