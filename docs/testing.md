# Scientific tests

The tests connect numerical results to independently calculated physics.
Run these commands from the repository root after [environment setup](getting-started.md):

```bash
make check
make physics
```

`make check` checks campaign inputs, data provenance, observable normalization,
blocking, and analysis. `make physics` builds the numerical code and checks
physical results against independent references. With separate Python environments:

```bash
make check PYTHON=python3 PYTHON_ED=/path/to/quspin/python
make physics PYTHON=python3 PYTHON_ED=/path/to/quspin/python
```

## What the physical checks establish

| Check | Independent reference | Physical quantities exercised |
| --- | --- | --- |
| [Finite-Fock ED](../tests/physics/test_ed_physics.py) | NumPy ladder matrices and a direct density-matrix trace in an 81-state, two-site Hilbert space | Both paired ED implementations: spectra, free energy, density, kinetic and interaction energies, pair amplitude, pair structure factor, and number moments |
| [Gaussian Fortran solvers](../tests/physics/test_gaussian_solvers.py) | Momentum-space Bogoliubov solution of the untruncated quadratic model | Both BAFQMC solvers: all four paper observables, anomalous amplitude, zero-pairing limit, pairing-phase transformation, and stabilization interval |
| [Fixed interacting HS field](../tests/physics/test_fixed_field_solver.py) | Direct NumPy time-slice products, inverses, and determinants for a prescribed nonuniform field | Number-conserving weights, Green-function diagnostics, site density, kinetic and interaction energy, and both imaginary-time sweep directions |
| [ED observable identities](../src/pairing/benchmarks/ed/test_observable_identities.py) | Thermodynamic derivatives and operator identities | Pair amplitude from the pairing derivative of free energy, kinetic energy from its hopping derivative, energy decomposition, and bosonic number normalization |
| [Numerical interfaces](../src/common/README.md) | Matrix reconstruction and known random-number recurrence | BLAS/LAPACK operations and the shared random generator used by both solvers |

The finite-Fock reference constructs its own Hamiltonian, operators, and thermal
trace. It compares the same occupation cutoff in both production ED representations.
The tests also verify the explicit flavor phase rotation taking
`Delta` to `-Delta`, and obtain particle number and grand energy from derivatives
with respect to chemical potential and inverse temperature.

The Gaussian checks run six fresh calculations on a periodic 3x3 triangular
lattice at `beta=1`, `mu=-5`, and zero interactions. Each calculation has four
measurement bins. They cover `Delta=0,+0.2,-0.2`, including the zero-pairing
comparison between the two solvers. Every bin is compared with the analytic
expectation for density, physical energy, and both K-point structure factors;
the paired solver also measures the anomalous amplitude. These tests resolve
the bosonic identity term and anomalous contractions in the density structure
factor. No stored benchmark measurements or ED occupation cutoff enter this comparison.

The fixed-field check turns on both interaction channels (`U1=0.1`, `U2=-0.7`)
and reads explicit nonuniform fields on six imaginary-time slices. Independent
dense propagation gives the determinant weight and Wick-contracted observables.
Zero proposal displacement holds the field fixed while the executable traverses
both time directions, at `Nwrap=1,2,6`. It also checks that the final fields equal
the prescribed configuration. Three fresh calculations take about 2 seconds.

Changing `Nwrap` checks agreement between stabilization intervals in both the
quadratic and prescribed interacting-field cases. These checks concern fixed-field
propagation and weights; acceptance decisions and strongly ill-conditioned
low-temperature products need targeted tests in those regimes.

The six Gaussian calculations took approximately 11 seconds with already-built
Intel MPI/MKL executables on the reference workstation. To run just this check,
without QuSpin:

```bash
BAFQMC_RUN_MPI_TESTS=1 python3 -m unittest discover -s tests/physics -p test_gaussian_solvers.py -v
```

The direct command skips live calculations unless the environment flag is set;
`make physics` enables them. Read the test summary when reporting which checks ran.

## Sampling and algorithm changes

Gaussian and finite-Fock comparisons are deterministic. Their tolerances describe
floating-point arithmetic and finite-difference accuracy. Interacting Monte Carlo
calculations have sampling uncertainty as well as Trotter and reference-cutoff
effects, which must be assessed at the chosen physical parameters.

For changes to HS fields, local updates, or the interacting propagation, run a
small interacting calculation with the relevant observable and an independent
reference. Both solvers provide live regression campaigns:

```bash
make -C src/number_conserving benchmark-dqmc
make -C src/pairing benchmark-dqmc
```

The number-conserving campaign includes four interacting 100000-bin cases and
takes approximately 10–15 minutes on the reference workstation. The paired
campaign exercises finite pairing and the zero-pairing limit. Choose further
cases to cover a changed interaction, geometry, or low-temperature regime.

Report means, block-based SEM, reference values, and residuals together with
sampling and cutoff settings. The paper reproduction retains all valid data
and uses no universal sigma threshold. Use deterministic identities for exact
contracts and the actual uncertainties when interpreting sampled comparisons.

## Workflow changes and new tests

Changes to production execution or continuation also need the small live pipeline
tests, which run both BAFQMC and ED:

```bash
BAFQMC_RUN_MPI_TESTS=1 python3 -m unittest discover -s tests/reproduction -p test_production.py -v
```

Set `BAFQMC_PYTHON_ED=/path/to/quspin/python` if ED uses a separate environment.
New physical tests should identify the operator or limiting case, compute an
independent expectation, and exercise the numerical output that can be wrong.
Keep seeds, sample counts, normalization, and any cutoff explicit. Generated
measurements belong in temporary or ignored output directories; compact input
fixtures and reference expectations can be committed.

Full paper reproduction is a separate production calculation, invoked by
`python3 reproduce.py`. Its [resource guide](../benchmarks/paper/RESOURCES.md)
gives the complete computing budget.
