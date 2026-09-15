# Algorithm and physical conventions

BAFQMC samples continuous Hubbard–Stratonovich (HS) fields in a
finite-temperature imaginary-time path integral. For each field configuration,
the bosonic problem becomes quadratic and its trace is evaluated through
single-particle matrices. Symmetry of the decoupled problem establishes
nonnegative weights, including for frustrated hopping. The implementations
here provide both number-conserving propagation and a full Nambu formulation
for onsite pairing.

## Implemented model

The computational model is a two-flavor Bose–Hubbard model on a periodic
triangular lattice, with positive nearest-neighbor hopping `t=1`, chemical
potential `mu`, and interactions

```text
U1 * sum_i (n_b,i + n_c,i)^2 + U2 * sum_i (n_b,i - n_c,i)^2.
```

The main benchmark sets `U1=0` and calls `U2` simply `U`. Its interaction
scan has zero pairing. The second main scan adds onsite pairing at `U=1`.
The supplemental scan varies the attractive total-density channel `U1` at
`U2=1`. The published parameters and reference definitions are listed in
[the benchmark guide](../benchmarks/paper/README.md).

For `U1<=0`, the total-density HS field has a real coefficient. For `U2>=0`,
the relative-density field has an imaginary coefficient. In the
number-conserving solver, the two flavor propagators are complex conjugates.
The paired solver evaluates the bosonic Gaussian trace in the full Nambu
basis `(b,c,b^+,c^+)`, retaining anomalous as well as normal Green functions.
The Fortran code performs local field updates, stabilized matrix propagation,
and Wick evaluation of observables; Python handles campaigns, reference
calculations, statistics, and plotting.

The main benchmark model has a finite trace throughout the auxiliary-field
domain under `mu < -3*t - abs(Delta)`, which every main benchmark point
satisfies. The supplemental attractive-density comparison uses its stated
finite-occupation reference. New physical regimes should carry their own
trace/convergence and sampling analysis.

## Mapping between paper and implementation

| Quantity | Implementation convention |
| --- | --- |
| Main-text interaction `U` | `U2=U`, `U1=0` |
| Pair term in the paper | `-Delta*(b^+ c^+ + b c)` |
| Pair term in code | `+Delta*(b^+ c^+ + b c)`, with `c_code=-c_paper` |
| Number of sites | `Ns=Lx*Ly` |
| Density | `density_total = <N>/Ns` |
| Physical energy | `energy_density = <H>/Ns`, excluding `-mu*N` |
| Plotted energy | `-E = -Ns*energy_density` |
| Structure factors | `S_SF(K)` and `S_DW(K)`, with `Ns^-2` normalization |
| Ordering momentum | `K=(4*pi/3,0)` |
| Plotted uncertainty | SEM from consecutive block means |

The pair phase leaves density, energy, and both plotted structure factors
unchanged. The anomalous amplitude changes sign:
`P_paper = -pair_equal`. Its physical energy contribution per site is
`-Delta*P_paper = Delta*pair_equal`.

The active lattice has one site per unit cell. Its forward bonds are
`(x,y)->(x+1,y)`, `(x,y)->(x,y+1)`, and `(x,y)->(x-1,y+1)`, with periodic
boundaries and the reverse hopping included. The existing K-point estimators
use lattice lengths divisible by three. Hopping is set by `RT` in each
solver's `src/calc_basic.f90`; the present command-line workflow is not an
arbitrary-lattice model builder.

## Exact references and errors

The number-conserving ED code uses particle-number and translation blocks;
the noninteracting benchmark uses a freshly evaluated analytic free-boson
reference. Pairing ED evaluates the thermal trace in an occupation-constrained
basis. Its paper setting, `nmax=3,ncut=4`, contains 7297 basis states. Cutoffs
and completed particle shells are recorded with the reference observables.

The workflow reports BAFQMC means, SEM, ED values, and their differences. It
retains valid data independently of a chosen sigma threshold. New calculations
can examine block-size dependence, independent seeds, Trotter steps, and ED
cutoffs as part of the scientific comparison.

For derivations and exact source mappings, see the
[number-conserving physics guide](solvers/number_conserving/physics.md),
[pairing physics guide](solvers/pairing/physics.md), and
[pairing observable contract](solvers/pairing/pairing_observable_contract.md).
