# Changelog

## 0.1.1 — 2026-09-15

- State the paper Hamiltonian, triangular lattice, supported parameters, and
  code-to-paper conventions in both READMEs and the algorithm guide.
- Document all active physical outputs, operator definitions, normalization,
  pair phase, and measurement availability with GitHub-rendered LaTeX.
- Add a model-development guide and agent skill for coordinated lattice,
  hopping, interaction, solver, ED, and observable extensions.
- Make README agent requests individually copyable.
- Remove archived-data checksum manifests and byte-level validation gates;
  retain numerical data checks and calculation continuation consistency.

## 0.1.0 — 2026-09-15

Initial public release of BAFQMC under the MIT license.

- Number-conserving and full Nambu bosonic auxiliary-field Monte Carlo
  solvers, with portable MPI Fortran and BLAS/LAPACK support.
- Exact-diagonalization references and complete production workflows for
  the main and supplemental benchmark figures: 22 parameter points.
- Original inputs, seeds, compact means and standard errors, block means,
  reference values, and reusable plotting code. Raw measurements are generated
  by the reproduction command.
- Custom-campaign interfaces, physics and development guides, numerical
  tests, and continuous integration.
- Independent finite-Fock ED, analytic Gaussian, fixed interacting-HS-field,
  and thermodynamic-identity tests, with real MPI execution in CI.
- Agent instructions, research task skills, an environment doctor, and a
  development container for GNU Fortran, MPI, and Python ED.
