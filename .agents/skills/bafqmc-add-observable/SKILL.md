---
name: bafqmc-add-observable
description: Implement or revise a physical observable across BAFQMC measurements, exact-diagonalization references, and statistical analysis, with a small-system verification.
---

# Add a physical observable

Read [AGENTS.md](../../../AGENTS.md) and the relevant solver's
`docs/solvers/<solver>/physics.md`. For paired observables also read
[the shared observable contract](../../../docs/solvers/pairing/pairing_observable_contract.md).
Use [the development guide](../../../docs/development.md) for verification.

First state the operator, flavor sums, equal-time ordering, momentum,
normalization, and any connected subtraction. Derive its Wick estimator using
the implemented Green function: the number-conserving matrix is
`<b_i b_j^+>`, so reversing its order includes the bosonic identity term.
The paired solver uses the full `(b,c,b^+,c^+)` Nambu ordering. Its anomalous
pair amplitude changes sign when expressed in the paper's operator phase.

Trace the entire measurement path before editing:

| Layer | Starting files |
| --- | --- |
| Accumulation, initialization, and reset | `src/*/src/obser_equal.f90` |
| MPI reduction, normalization, and output | `src/*/src/fourier_trans.f90` |
| Number-conserving ED | `src/number_conserving/benchmarks/ed/EDtriangle_quspin_3x3.py` |
| Paired ED | `src/pairing/benchmarks/ed/ed_pairing_triangle_general.py` |
| Parsing, blocking, and comparisons | `src/number_conserving/benchmarks/campaign_analysis.py`, `src/pairing/benchmarks/pairing_delta_analysis.py` |

Keep the operator definition beside the estimator and add its ED counterpart
with the same normalization. Extend the applicable output schema and parser.
The paper's `benchmarks/paper/analysis.py` exports four fixed benchmark
observables; extend it only if the task includes a change to those exports.
A new research observable does not require rewriting published data or figures.

For structure factors, retain reciprocal-lattice coordinates and check the
allowed momentum for the chosen size. For energy, distinguish physical
energy per site from total energy and from the chemical-potential contribution.
Preserve complex measurements where needed and report the physical channel
with block-based SEM.

Verify a meaningful small-system identity, analytic limit, or independent ED
representation. Add a focused test that would catch the wrong operator order,
normalization, or phase. Existing useful checks include:

```bash
make physics
```

Run the checks applicable to the changed solver, plus a live comparison that
actually measures the new observable. For a separate ED environment use
`PYTHON_ED=/path/to/quspin/python` with `make physics`.
Retain all valid means and statistical uncertainties; do not turn a stochastic
residual into an automatic three-SEM rejection of the paper reproduction.

Deliver the definition, derivation or identity, implementation, repeatable
small input, and observed BAFQMC/ED comparison. Track compact fixtures where
they support verification; keep generated raw chains and build output ignored.
