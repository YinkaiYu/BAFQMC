# Get started

Run the computational workflow in Linux or WSL. All commands below start in
the repository root. On WSL, a checkout on the Linux filesystem gives better
performance for the many small measurement writes.

## Install the environment

Install an MPI Fortran compiler, Make, and BLAS/LAPACK. For example, on
Ubuntu or Debian:

```bash
sudo apt-get update
sudo apt-get install -y git make gfortran openmpi-bin libopenmpi-dev libblas-dev liblapack-dev
```

Create the supplied Python environment with Conda:

```bash
conda env create -f benchmarks/paper/environment.yml
conda activate bafqmc
python3 scripts/doctor.py
```

The environment contains NumPy, Matplotlib, and QuSpin for ED. An existing
Python 3.11 environment can instead install the two requirements files:

```bash
python3 -m pip install -r benchmarks/paper/requirements.txt -r benchmarks/paper/requirements-ed.txt
```

The Makefiles use `mpifort` with GNU Fortran by default. An activated Intel
MPI/MKL environment is also supported. Build configuration and optional
overrides are described in [the shared numerical guide](../src/common/README.md).
For a separate QuSpin environment, pass `--python-ed /path/to/python` to the
root reproduction command and environment doctor.

## Check the installation

```bash
make check
python3 reproduce.py --mode smoke
```

The smoke workflow builds both solvers, executes small Monte Carlo runs,
computes exact references, and checks analysis. It writes results under
`benchmarks/paper/output/smoke/`. For another run, select a fresh destination:

```bash
python3 reproduce.py --mode smoke --output /tmp/bafqmc-smoke-second
```

The smoke cases use reduced sampling and ED basis sizes. They verify the
installation before production.

## Reproduce all benchmark data

```bash
python3 reproduce.py --plan
python3 reproduce.py --output benchmarks/paper/output/paper-run
```

The second command builds and runs BAFQMC, computes the ED references, then
creates tables and figures from the new measurements. It includes all 22
points: the 15 main benchmark points and seven supplemental points. Each
production point uses 100000 measurement bins and the input settings listed in the paper manifest.
For the reference desktop, budget **12–24 hours**, **16 GiB RAM** with about
8 GiB available to the calculation, and **8 GiB free disk**. See
[measured resource estimates](../benchmarks/paper/RESOURCES.md) for details.

Outputs in the selected directory include:

```text
runs/                                raw BAFQMC data and ED results
progress.json                        completed stages and timing
observables.csv                      means, SEM, and ED values
block_means.csv                       blocked measurements
records.json                         results with parameters and comparison values
figures/benchmark_combined.pdf        main benchmark
figures/benchmark_attractive.pdf      supplemental benchmark
environment.json                     environment and selected-case metadata
```

Continue after an interruption with the same output directory and options:

```bash
python3 reproduce.py --output benchmarks/paper/output/paper-run --resume
```

Completed stages are checked and reused. The interrupted stage is rerun from
the manifest inputs. The root runner manages scratch automatically; use
`--work-dir /path/to/empty/linux-scratch` to choose its location.

To compute a selected part of the paper:

```bash
python3 reproduce.py --scope main
python3 reproduce.py --scope supplement
python3 reproduce.py --scope main --model pairing
```

The standalone `python3 reproduce.py` command always computes the full
campaign by default. `--mode plot` is available to inspect the stored
processed data. A list of all stages and data conventions is in
[the reproduction guide](../benchmarks/paper/README.md).

## Start your own calculation

Tell your agent the Hamiltonian, lattice size, temperature, observables,
desired parameter scan, and available resources. Have it read
[AGENTS.md](../AGENTS.md) and follow
[the custom-campaign recipes](agent-workflows.md). The solver-specific runners
accept separate manifests, so research inputs can be developed independently
of the published benchmark package.
