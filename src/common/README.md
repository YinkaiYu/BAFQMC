# Portable numerical support

Both Fortran implementations compile `mymats.f90` and `random.f90` from this
folder. They require MPI and a BLAS/LAPACK implementation. No NAG, EISPACK,
LINPACK, or other third-party library source is bundled here.

The original solver used `MyMats` from `Lib_90_new/Modules/mat_mod.f90` for a
small set of matrix operations and the `ranf` interface from
`Lib_90_new/Ran/ran_imada.f`. These two files are newly written implementations
of the required interfaces, under the software distribution's MIT license.

| Interface | Bundled implementation | Historical backend |
| --- | --- | --- |
| `mmult(C,A,B)` | BLAS `ZGEMM` | BLAS `ZGEMM` |
| `diag(A,U,W)` | LAPACK `ZHEEV`, complex Hermitian | EISPACK `CH` |
| `inv(A,Ainv,det)` | LAPACK `ZGETRF`, `ZGETRI` | LINPACK `ZGEFA`, `ZGEDI` |
| `udv(A,U,D,V,ncon)` | LAPACK QR, `A = U diag(D) V` | NAG QR |
| `ranf(seed)` | Explicit 64-bit modular arithmetic | Implicit 32-bit overflow |

The uniform random recurrence remains

```text
seed_next = (48828125 * seed) mod 2147483648
u = seed_next / 2147483648
```

The archived initial seeds and recurrence are retained verbatim. Do not use
a zero seed for new runs. Changes in eigensolver roundoff and compiler or
BLAS behavior can still change an interacting Monte Carlo trajectory.
Recomputed error bars and means should be assessed statistically.

The active stabilized propagation and local-update Fortran sources are
unchanged. The `udv` adapter supports the retained alternative stabilization
routine; the production path uses its own pivoted LAPACK QR code.

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
