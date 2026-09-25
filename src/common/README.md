# Portable numerical support

Both Fortran implementations compile `mymats.f90` and `random.f90` from this
folder. They require MPI and a BLAS/LAPACK implementation. No NAG, EISPACK,
LINPACK, or other third-party library source is bundled here.

The portable build provides the small matrix-operation and random-number
interfaces required by both solvers. The bundled implementations use standard
BLAS/LAPACK kernels and are released under the MIT license.

| Interface | Bundled implementation |
| --- | --- |
| `mmult(C,A,B)` | BLAS `ZGEMM` |
| `diag(A,U,W)` | LAPACK `ZHEEV`, complex Hermitian |
| `inv(A,Ainv,det)` | LAPACK `ZGETRF`, `ZGETRI` |
| `udv(A,U,D,V,ncon)` | LAPACK QR, `A = U diag(D) V` |
| `ranf(seed)` | Explicit 64-bit modular arithmetic |

The uniform random recurrence remains

```text
seed_next = (48828125 * seed) mod 2147483648
u = seed_next / 2147483648
```

Do not use a zero seed for new runs. Changes in eigensolver roundoff and compiler or
BLAS behavior can still change an interacting Monte Carlo trajectory.
Recomputed error bars and means should be assessed statistically.

The active stabilized propagation and local-update Fortran sources use the same
interfaces described above. The production path uses pivoted LAPACK QR code.

## Build and check

With GNU Fortran, an MPI development package, and BLAS/LAPACK installed:

```bash
make -C src/common check
make -C src/number_conserving build
make -C src/pairing build
```

`compiler.mk` selects `mpifort` when available. Compiler and linker settings
can be overridden, for example `FC=mpifort LDLIBS='-lopenblas'`. The GNU build
uses `-fallow-argument-mismatch` for the inherited `mpif.h` calling convention.
For Intel, first load the Intel oneAPI compiler/MPI environment, then use
`FC=mpiifx`; the default library switch is `-qmkl`.

`make check` verifies Hermitian eigendecomposition, nonsymmetric complex
inversion and determinant, QR reconstruction and orthogonality, and the
uniform random sequence. Each mode also retains live physical benchmark
checks and observable tests.

## Numerical validation

The numerical interface checks passed with Intel `mpiifx` 2025.2.1 and MKL.
Both the number-conserving and finite-pairing live smoke benchmarks passed.
The complete number-conserving regression suite also passed: one analytic
free-boson case and four interacting cases with 100000 bins each.
The number-conserving eight-bin 3x3 free-boson pipeline reproduced the density,
energy density, and two plotted structure factors to within `3e-15` of their
analytic values. GNU/BLAS build settings are supplied, but a GNU compiler was
not available in this validation environment.
