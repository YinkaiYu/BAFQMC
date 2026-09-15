# Models, lattices, and physical conventions

BAFQMC samples continuous Hubbard–Stratonovich (HS) fields in a
finite-temperature imaginary-time path integral. At each field configuration,
the bosons are traced out through quadratic single-particle propagation.
Time-reversal symmetry (TRS) or reflection positivity (RP) of the decoupled
problem establishes nonnegative weights, including for frustrated hopping.
The supplied implementations treat number-conserving bosons and onsite pairing.

## Hamiltonian implemented by the solvers

We use the paper's operators and signs: two boson flavors $`\hat b_i`$ and
$`\hat c_i`$, with $`\hat n_{b,i}=\hat b_i^+\hat b_i`$ and
$`\hat n_{c,i}=\hat c_i^+\hat c_i`$. The main benchmark Hamiltonian is

```math
\begin{aligned}
\hat H={}&\hat H_t+\hat H_U+\hat H_\Delta,\\
\hat H_t={}&t\sum_{\langle ij\rangle}
\left(\hat b_i^+\hat b_j+\hat b_j^+\hat b_i
+\hat c_i^+\hat c_j+\hat c_j^+\hat c_i\right),\\
\hat H_U={}&U\sum_i(\hat n_{b,i}-\hat n_{c,i})^2,\\
\hat H_\Delta={}&-\sum_i\left(\Delta\,\hat b_i^+\hat c_i^+
+\Delta^*\,\hat b_i\hat c_i\right).
\end{aligned}
```

The ensemble is grand canonical, at inverse temperature $`\beta`$:

```math
\hat N=\sum_i(\hat n_{b,i}+\hat n_{c,i}),\qquad
\hat H_\mu=\hat H-\mu\hat N,\qquad
Z=\mathrm{Tr}e^{-\beta\hat H_\mu}.
```

The number-conserving solver uses $`\Delta=0`$. The paired solver accepts real
$`\Delta`$ and retains normal and anomalous Green functions. Its internal
operator phase is $`\hat b_{\mathrm{code}}=\hat b_{\mathrm{paper}}`$,
$`\hat c_{\mathrm{code}}=-\hat c_{\mathrm{paper}}`$, so the code's pair term
has a plus sign at the same input $`\Delta`$. This is the same physical model.
The implementation's Nambu order is $`(b,c,b^+,c^+)`$ in that code basis.

Both solvers also implement the two density channels used in the Supplemental
Material, replacing $`\hat H_U`$ by

```math
\begin{aligned}
\hat H_U={}&\sum_i\left[
U_1(\hat n_{b,i}+\hat n_{c,i})^2
+U_2(\hat n_{b,i}-\hat n_{c,i})^2\right]\\
={}&\sum_i\left[(U_1+U_2)(\hat n_{b,i}^2+\hat n_{c,i}^2)
+2(U_1-U_2)\hat n_{b,i}\hat n_{c,i}\right].
\end{aligned}
```

The main model sets $`U_1=0`$, $`U_2=U`$. The continuous-HS implementation uses
$`U_1\leq0`$ and $`U_2\geq0`$.

The squared densities include their linear number terms. When translating
from a model written with $`n(n-1)`$, carry the resulting chemical-potential
shift into both BAFQMC and ED.

## Triangular lattice

The active geometry has one site per primitive cell, with periodic boundaries:

```math
\begin{gathered}
\mathbf a_1=(1,0),\qquad
\mathbf a_2=\left(\frac12,\frac{\sqrt3}{2}\right),\\
\mathbf r_{x,y}=x\mathbf a_1+y\mathbf a_2,\qquad
N_s=L_xL_y,\\
(x,y)\equiv(x+L_x,y)\equiv(x,y+L_y).
\end{gathered}
```

Each cell supplies the three forward bonds

```math
(x,y)\longrightarrow(x+1,y),\quad
(x,y)\longrightarrow(x,y+1),\quad
(x,y)\longrightarrow(x-1,y+1).
```

Their Hermitian conjugates supply the reverse directions. On very small
periodic cells, distinct bonds can connect the same pair of sites or return
to the original site; their contributions are added. The independent ED and
fixed-field tests use this bond multiplicity explicitly.

With the implemented positive hopping $`t=1`$, the triangular loops are
frustrated. The dispersion is

```math
\varepsilon(\mathbf k)=2t\left[
\cos(\mathbf k\cdot\mathbf a_1)
+\cos(\mathbf k\cdot\mathbf a_2)
+\cos\bigl(\mathbf k\cdot(\mathbf a_2-\mathbf a_1)\bigr)\right].
```

The reciprocal basis and one band minimum are

```math
\begin{gathered}
\mathbf b_1=2\pi\left(1,-\frac1{\sqrt3}\right),\qquad
\mathbf b_2=2\pi\left(0,\frac2{\sqrt3}\right),\\
\mathbf b_i\cdot\mathbf a_j=2\pi\delta_{ij},\qquad
\mathbf K=\frac23\mathbf b_1+\frac13\mathbf b_2
=\left(\frac{4\pi}{3},0\right).
\end{gathered}
```

Allowed momenta are $`\mathbf k=(m/L_x)\mathbf b_1+(n/L_y)\mathbf b_2`$.
The existing $`K`$-point estimators require both lengths to be multiples of
three. On other sizes their Fortran output remains zero for these channels;
choose a compatible momentum and implement its estimator for a new study.
Density and energy measurements remain available on other lattice sizes.

## Parameters and current interfaces

| Physical quantity | Existing control |
| --- | --- |
| Lattice dimensions $`L_x,L_y`$ | Runtime inputs `Nlx,Nly`; BAFQMC accepts different positive lattice lengths |
| Temperature $`\beta`$ and time step $`\Delta\tau=\beta/L_\tau`$ | Runtime inputs `Beta,Ltrot` |
| Density interactions $`U_1,U_2`$ and chemical potential $`\mu`$ | Runtime inputs `RU1,RU2,mu` |
| Real onsite pairing $`\Delta`$ | Runtime input `RDelta` in the paired solver |
| Hopping $`t`$ | Set to `RT=1.d0` by `Params_set` in each solver's `src/calc_basic.f90`; change the implementation and rebuild to vary it |
| Crystal geometry, bond amplitudes, extra interactions or pairing patterns | Model development; use the [extension guide](model-development.md) |

A JSON `t` field alone does not change the current Fortran hopping. The
provided number-conserving interacting ED workflow targets the $`3\times3`$
benchmark lattice; the paired ED implementations accept small $`L_x\times L_y`$
cells subject to their basis cutoffs. Extending a BAFQMC calculation to a new
geometry includes choosing and implementing a matching reference calculation.

## From interactions to auxiliary fields

Write $`n_+=n_b+n_c`$ and $`n_-=n_b-n_c`$. The local HS identities are

```math
\begin{aligned}
e^{-\Delta\tau U_1 n_+^2}
&=\int\frac{d\phi_1}{\sqrt{2\pi}}e^{-\phi_1^2/2}
  e^{\sqrt{-2U_1\Delta\tau}\,\phi_1n_+},\\
e^{-\Delta\tau U_2 n_-^2}
&=\int\frac{d\phi_2}{\sqrt{2\pi}}e^{-\phi_2^2/2}
  e^{\mathrm i\sqrt{2U_2\Delta\tau}\,\phi_2n_-}.
\end{aligned}
```

For number-conserving propagation the two flavors see conjugate onsite
potentials. Consequently $`B_{c,\ell}=\overline{B_{b,\ell}}`$, and the code
can obtain the $`c`$-sector Green function from the $`b`$ sector by conjugation.
The paired implementation evaluates the full Nambu Gaussian trace, including
normal-ordering scalar factors and the determinant-square-root
weight. Local field updates, stabilized propagation, and Wick estimators
are implemented in Fortran; Python handles ED, campaigns, and analysis.

For the main model $`U_1=0`$, the sufficient condition
$`\mu<-3t-|\Delta|`$ gives a finite trace throughout the auxiliary-field domain;
all main benchmark points satisfy it. The attractive $`U_1<0`$ supplemental
benchmark uses its separately specified finite-occupation comparison.
For a new Hamiltonian, establish the trace domain and the appropriate HS
symmetry together with its implementation, as explained in the
[model-development guide](model-development.md).

## Benchmark models and observables

| Scan | Hamiltonian parameters | Geometry and ensemble |
| --- | --- | --- |
| Main interaction scan | $`U_1=0`$, $`U_2=U`$, $`\Delta=0`$ | $`3\times3`$, $`t=1`$, $`\beta=4`$, $`\mu=-3.5`$ |
| Main pairing scan | $`U_1=0`$, $`U_2=U=1`$, variable $`\Delta`$ | $`3\times3`$, $`t=1`$, $`\beta=4`$, $`\mu=-5`$ |
| Supplemental density-channel scan | Variable $`U_1\leq0`$, $`U_2=1`$, $`\Delta=0`$ | $`3\times3`$, $`t=1`$, $`\beta=1`$, $`\mu=-7`$ |

The paper writes the pair term with a minus sign. Its operators and the
implementation's operators are related by
$`b_{\mathrm{code}}=b_{\mathrm{paper}}`$ and
$`c_{\mathrm{code}}=-c_{\mathrm{paper}}`$. Density, energy, and both benchmark
structure factors are invariant under this phase change; the anomalous pair
amplitude changes sign.

For $`\hat n_i=\hat n_{b,i}+\hat n_{c,i}`$, the main benchmark observables are

```math
\begin{aligned}
\rho&=\frac{\langle\hat N\rangle}{N_s},\qquad E=\langle\hat H\rangle,\\
S_{\mathrm{SF}}(\mathbf k)&=\frac1{N_s^2}\sum_{ij}
e^{\mathrm i\mathbf k\cdot(\mathbf r_i-\mathbf r_j)}
\left\langle\hat b_i^+\hat b_j+\hat c_i^+\hat c_j\right\rangle,\\
S_{\mathrm{DW}}(\mathbf k)&=\frac1{N_s^2}\sum_{ij}
e^{\mathrm i\mathbf k\cdot(\mathbf r_i-\mathbf r_j)}
\left\langle\hat n_i\hat n_j\right\rangle.
\end{aligned}
```

The [output-observable guide](observables.md) defines every current physical
output and its normalization. `energy_density` is $`E/N_s`$, excluding the chemical-potential term; the paper
plots $`-E=-N_s\,\texttt{energy\_density}`$. The code's anomalous amplitude is
$`\texttt{pair\_equal}=N_s^{-1}\sum_i\langle b_i^+c_i^++b_ic_i\rangle_{\mathrm{code}}`$,
so its pairing energy per site is $`\Delta\,\texttt{pair\_equal}`$.
The paper's anomalous amplitude is the negative of this value.

The workflow reports means, block-based standard errors of the mean (SEM),
reference values, and differences. The main paired ED setting
`nmax=3,ncut=4` contains 7297 basis states; number-conserving ED accumulates
particle-number and translation sectors. Noninteracting references are
also calculated directly from the free-boson dispersion. Each comparison
records its sampling parameters and ED cutoffs.

See [benchmark inputs and settings](../benchmarks/paper/README.md),
[scientific tests](testing.md), and the detailed
[number-conserving](solvers/number_conserving/physics.md) and
[paired](solvers/pairing/physics.md) implementation guides.
