# Number-conserving BAFQMC and exact diagonalization

This directory contains the Fortran finite-temperature BAFQMC solver for the
two-flavor triangular-lattice Bose-Hubbard Hamiltonian used in the manuscript,
the QuSpin ED implementations, regression cases, and campaign
analysis tools. See [the physics guide](../../docs/solvers/number_conserving/physics.md) for the Hamiltonian,
Hubbard–Stratonovich fields, Green functions, estimators, and cutoff conditions.

The solver produces the relative-density scan with `U1=U` and `U2=0`, and the
total-density scan with varying `U2` at `U1=1`. Both have `Delta=0`. The input
names `U1` and `U2` identify the
relative-density and total-density channels throughout the implementation.

## Build and run

On Linux/WSL install an MPI Fortran compiler and BLAS/LAPACK; the default build
supports GNU Fortran with `mpifort`. Python analysis requires NumPy and plotting
requires Matplotlib. ED additionally requires QuSpin, SciPy, and Numba.
The [installation guide](../../docs/getting-started.md) provides dependency setup.

```bash
make -C src/number_conserving build
make -C src/number_conserving run-example
make -C src/number_conserving benchmark-fast
make -C src/number_conserving check-fixtures
```

`run-example` copies the eight-bin 3x3 free-boson example into `build/example`, then
runs one MPI rank. It refuses an existing nonempty output directory. Set
`RUN_DIR=/path/to/new/run` for another output location. To develop a new model
point, copy the three input files from `examples/number_conserving/examples/triangle_3x3_free` into a new
run directory, edit `paramC_sets.txt`, and run:

```bash
bash src/number_conserving/scripts/run_local.sh /path/to/new/run 1
```

The executable reads `paramC_sets.txt`, `confin.txt`, and `seeds.txt` in its
working directory. It appends measurements there and writes restart state.
Do not rerun into an old output directory when starting a fresh chain.

The [shared numerical adapter](../common/README.md) supplies the BLAS/LAPACK
interfaces used by the portable build. Processed benchmark data and their
initial inputs are included under `benchmarks/paper/data/`.

## Paper production entry point

The root reproduction command prepares the paper manifest. This mode-specific
entry point can also run any selected case independently:

```bash
python src/number_conserving/run_paper.py \
  --manifest /path/to/manifest.json --output /path/to/new/results --mode dqmc
python src/number_conserving/run_paper.py \
  --manifest /path/to/manifest.json --output /path/to/new/results --mode ed
python src/number_conserving/run_paper.py \
  --manifest /path/to/manifest.json --output /path/to/new/results --mode analyze
```

From the repository root, `python3 reproduce.py --scope combined --model number_conserving`
reproduces the eight relative-density U1 points, and `python3 reproduce.py --scope total_density`
reproduces the seven total-density U2 points. Their physical and sampling
parameters are listed in [the reproduction guide](../../benchmarks/paper/README.md).

Use a QuSpin-enabled interpreter for ED, or pass `--python /path/to/python`.
Use `--case ID` (repeatable) to select cases and `--dry-run` to inspect commands.
Manifest paths are relative to the manifest file:

```json
{
  "analysis": {"block_size": 1000, "skip_samples": 0},
  "cases": [{
    "id": "example",
    "input_dir": "inputs/example",
    "ed_params": "inputs/example/params.json",
    "parameters": {
      "Lx": 3, "Ly": 3, "t": 1.0,
      "U1": 1.5, "U2": 0.0, "beta": 4.0, "mu": -3.5
    }
  }]
}
```

Outputs are `OUTPUT/ID/dqmc/`, `OUTPUT/ID/ed/results.json`, and
`OUTPUT/analysis.json`. ED parameters and convergence policy are read from the
manifest. Stable free-boson points generate a fresh exact analytic reference in the
`ed` stage, matching the manuscript prescription. Analysis exports the means,
standard errors, and ED differences for all valid samples. Missing or malformed
output is an error.

The full 3x2 live regression suite is `make benchmark-dqmc`; it uses the
100000-bin interacting cases and can take approximately 10–15 minutes.
This suite is separate from the longer 3x3 paper production campaign.

## Source map

- `src/`: model, HS fields, local updates, stabilized propagation, observables.
- `benchmarks/ed/EDtriangle_quspin_3x3.py`: paper ED using particle-number and
  translation blocks, shell convergence diagnostics, and incremental results.
- `benchmarks/ed/EDtriangle_symm_NEblock.py`: 3x2 ED reference driver.
- `benchmarks/campaign_analysis.py`: parsers, blocking, normalization,
  exact free reference, and reliability diagnostics.
- `benchmarks/campaign_3x3.py`: general campaign input/report helper;
  its default campaign is broader than the selected manuscript figure.
- `../../examples/number_conserving/`: compact example, smoke, and live-regression inputs.
- `../../tests/number_conserving/`: Python unit and optional runtime tests.
- `../../docs/solvers/number_conserving/`: physics and development guides.
- `benchmarks/fixtures/`, `references/`, `dqmc_references/`: regression data.
- `../common/`: shared portable numerical adapter.

The plotting workflow and paper inputs/data are documented at the
repository's reproduction entry point.

For a fast complete orchestration check, pass
`--manifest src/number_conserving/benchmarks/campaigns/pipeline_smoke.json`
to the three modes above. It runs eight free-boson bins and computes the exact free reference.
The optional QuSpin runtime tests separately exercise a tiny interacting ED
calculation. The interacting 3x2 example is
retained separately; it is not the default introductory run.

Validation: the five-case live regression suite passed
with the portable Intel/MKL build (all four interacting cases retain 100000
bins). The Python test suite and the optional two MPI runtime tests, tiny QuSpin
runtime test, and memory-guard checkpoint test passed.

The ED driver and `run_paper.py` accept `--dense-memory-cap-gib` (default 12).
Before dense diagonalization, a conservative four-complex-matrix estimate is
checked; exceeding the budget stops computation while preserving the last
completed particle shell. This guard does not change physical cutoffs. The
largest stored interacting paper block has dimension 9075 (shell seven),
with a 4.91 GiB dense-workspace estimate; shell eight would reach dimension
27225 (44.18 GiB). A serial full reproduction should use a machine with at
least 16 GiB physical RAM and approximately 8 GiB available to the calculation.
