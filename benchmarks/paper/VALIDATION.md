# Numerical validation record

The computational package was validated in Linux/WSL on
2026-09-15. The Fortran tests used Intel `mpiifx` 2025.2.1 and MKL. Python ED
used QuSpin 1.0.1, NumPy 2.4.4, SciPy 1.17.1, and Numba 0.65.1.

- All Fortran model, propagation, and measurement sources compile with the
  shared numerical interfaces described in `src/common/README.md`.
- Matrix reconstruction, Hermitian eigensystems, inversion, and the random
  sequence passed independent invariant checks. The modular random recurrence
  also matched 150000 transitions under the integer-wrap convention.
- The number-conserving live regression completed its analytic free case and
  four interacting 100000-bin cases. Its Python, parameter-routing, and ED
  memory checks passed.
- Paired live checks covered finite Delta and Delta=0. Independent small ED
  formulations agreed on 25 observables and energy derivatives.
- Small complete DQMC -> ED -> analysis runs passed for both solvers. Production
  orchestration checks cover per-case copying, hashes, continuation, and
  restarting an interrupted chain from the manifest input.
- A full paper-size paired ED calculation (7297 basis states at Delta=0.2)
  reproduced the four stored reference values within 5.3e-18 absolute.
- All 22 stored cases were reblocked from their 100000-sample measurement
  sequences before retaining the compact processed tables. The stored means
  and standard errors agree with that independent reconstruction. The compact
  package also passed without any raw chains present.
These checks cover the solvers, reference calculations, stored data, and
production orchestration. The full 22-point campaign is the default reproduction
command; its time and memory budget is based on the stage measurements in
[RESOURCES.md](RESOURCES.md).

GNU Fortran/BLAS build settings and a GitHub Actions workflow are supplied.
That compiler was not installed on the validation workstation, so local runtime
results above refer to the Intel/MKL build. CI exercises the GNU setup when run.

## Reproduction checks

The reproduction interface was checked against the benchmark data layout on
2026-09-15. The default production comprises eight relative-density points,
seven pairing points, and seven total-density points.

- All 13 reproduction tests passed, including real eight-bin BAFQMC/ED runs,
  continuation and interrupted-chain recovery, and selecting a single case
  for fresh reference calculation. Changing the selected cases is detected
  when resuming a campaign.
- Scope selection preserves the supplied inputs, case ordering, and seeds.
  The 88 stored observable rows retain every numeric value; added columns
  identify the current manuscript figures and notation.
- A 172-state paired ED calculation at Delta=+0.2 and -0.2 gave identical
  density, physical energy, structure factors, pairing energy, and free energy.
  The anomalous `pair_equal` amplitude reversed sign, confirming the documented
  phase conversion between the solver and manuscript conventions.
