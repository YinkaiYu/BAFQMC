# Pairing BAFQMC and exact diagonalization

This solver samples the two-flavor triangular-lattice Bose–Hubbard model with
onsite pairing. Its full four-sector Nambu Green matrix and determinant square
root implement the paired construction used in the manuscript. The Fortran physics kernels are maintained directly in this repository under the manuscript notation.

Commands below run from `src/pairing/` in Linux or WSL. See the repository
[installation guide](../../docs/getting-started.md) for environment setup and the
root README for the combined paper reproduction command.

This solver produces panels (e–h) of the combined main-text benchmark: seven
Delta points at `U1=1`, `U2=0`, `mu=-5`, and `beta=4`. Inputs use the manuscript
convention: `U1` multiplies the repulsive relative-density square and `U2`
multiplies the attractive total-density square. Thus the main scan uses
`U1=U` and `U2=0`; the supplemental two-channel scan keeps `U1=1` and
varies non-positive `U2`.
The solver's positive real pair coefficient is related to the manuscript's
negative pair term by `c_code=-c_paper`. All four plotted observables are
unchanged by this phase convention; the anomalous amplitude `pair_equal`
has the opposite sign in manuscript operators. See [the Hamiltonian and
phase convention](../../docs/solvers/pairing/physics.md).

## Build and first run

```bash
make build
make run-example
```

The build uses the shared `../common/` LAPACK adapters and random-number
generator. It requires an MPI Fortran compiler and BLAS/LAPACK. Intel MPI/MKL
is also supported after its environment has been activated. Override `FC`,
`FFLAGS`, `LDFLAGS`, or `LDLIBS` through Make when needed.

`make run-example` copies the committed inputs to `build/example/`, then runs
one MPI rank. Choose another empty destination on subsequent runs:

```bash
make run-example RUN_DIR=build/my-example MPI_NP=1
```

Modify the copied `paramC_sets.txt` to explore other model or Monte Carlo
parameters; its seven numeric rows are documented in
[the development guide](../../docs/solvers/pairing/development.md). The hopping coefficient is
`RT=1` in `src/calc_basic.f90`. Run output uses fixed filenames in the working
directory. The executable appends observables, so each independent chain needs
a fresh directory. `seeds.txt` contains one scalar seed per MPI rank.

## Paper campaign and analysis

`run_paper.py` accepts the supplied paper manifest or a custom manifest using
the same schema. The stages use a shared output directory:

```bash
python run_paper.py --manifest ../../benchmarks/paper/data/pairing/manifest.json --output build/paper --mode init
python run_paper.py --manifest ../../benchmarks/paper/data/pairing/manifest.json --output build/paper --mode dqmc
python run_paper.py --manifest ../../benchmarks/paper/data/pairing/manifest.json --output build/paper --mode ed
python run_paper.py --manifest ../../benchmarks/paper/data/pairing/manifest.json --output build/paper --mode analyze
```

The production campaign has seven pairing values and 100000 measurement bins
per value. These are production calculations; use the repository smoke command
for a quick installation check. From the repository root,
`python3 reproduce.py --scope main --model pairing` runs this scan and produces
`benchmark_combined_pairing.pdf` with the manuscript's panel labels (e–h).
The default seed policy initializes each case from `base_seed + case_index` and generates
its per-rank seed list. Use `--seed` to start independent chains and `--np` to
set MPI ranks. Results vary with compiler, numerical libraries, and rank count.
Use `--dry-run` to inspect parameters and per-case seeds without writing files.
Repeat `--case CASE_NAME` to execute selected cases; initialization always uses
the full manifest, preserving seed positions for separately scheduled cases. On
WSL, choose an output directory on the Linux filesystem
instead of `/mnt/c`: the solver appends many small observable records.

A calibration on an AMD Ryzen 5 9600X under WSL2, using Intel Fortran/MKL,
one MPI rank and one numerical-library thread, took 46.38 seconds and 65 MiB
peak RSS for 1000 bins at the paper's physical parameters, including the
500 warmup iterations. Linear scaling gives approximately nine hours
for all seven 100000-bin runs. The measured production chains took 12.2 hours
of aggregate case elapsed time. Allow roughly 9–13 hours for the pairing
BAFQMC campaign on a similar machine, with variation from system load and
storage. Its full observable output occupies about 683 MiB.

The stages write `inputs/<case>/`, `ed_results/<case>.json`, and
`summary/comparison_observables.csv`. The JSON comparison additionally retains
block errors, imaginary-part diagnostics and cutoff-sensitive observables.
The reference file `comparison_dqmc_ed_nmax3_ncut4.json` matches the paper
plot input; consult the embedded manifest for the cutoffs of a custom run.
Analysis exports all valid means, standard errors, and reference differences.
`--require-agreement` optionally enables the inherited three-standard-error check.

Energy files store physical energy per site, excluding the chemical-potential
term. The paper plots total physical energy with a minus sign:
`-E = -Lx * Ly * energy_density`.

## ED and verification

Activate the repository QuSpin environment, or pass its interpreter with
`--python /path/to/python` to `run_paper.py`. Standalone ED usage is:

```bash
python benchmarks/ed/ed_pairing_triangle_general.py --params benchmarks/ed/params_triangle_pairing_smoke.json --output build/ed-smoke.json
make benchmark-ed
make check-fixtures
make benchmark-fast
```

`make benchmark-ed` checks two independent tiny-basis ED representations and
finite-difference identities for measured observables. `make benchmark-fast`
runs a small live BAFQMC comparison. `make benchmark-dqmc` also checks the
zero-pairing limit against the included number-conserving reference.

Pairing ED performs the full thermal trace within the stated occupation
cutoffs. The paper uses `nmax=3,ncut=4` (7297 basis states), with an estimated
3.2 GiB dense-diagonalization peak. Cases run serially; the default memory cap
is 12 GiB. Increase `--dense-memory-cap-gib` only to match available memory.
Larger cutoffs grow rapidly and are separate scientific calculations.
For the same machine and one thread, an actual paper-size ED case at
`Delta=0.2` took 40.40 seconds and 3.01 GiB peak RSS; its four paper observables
agreed with the stored ED values to within `5.3e-18` absolute. The seven ED
cases therefore require approximately five minutes under those conditions,
in addition to the BAFQMC time. These timings were measured on 2026-09-15.

See [the physics guide](../../docs/solvers/pairing/physics.md) for the Hamiltonian, Nambu conventions,
HS scalar factor, Green-function blocks and Wick estimators;
[the observable contract](../../docs/solvers/pairing/pairing_observable_contract.md)
defines the observable correspondence with the number-conserving solver.
