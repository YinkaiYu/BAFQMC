# Physics Notes

Implementation paths and commands in this guide are relative to
`src/pairing/` from the repository root.

This solver implements finite-temperature grand-canonical DQMC for a two-flavor Bose-Hubbard model with onsite pairing on a triangular lattice. The active triangular implementation has one orbital per cell, so $`N_s=L_xL_y`$.

## Hamiltonian

The paper's main-text model has the physical Hamiltonian

```math
\begin{aligned}
\hat H={}&t\sum_{\langle ij\rangle}
\left(\hat b_i^+\hat b_j+\hat b_j^+\hat b_i
+\hat c_i^+\hat c_j+\hat c_j^+\hat c_i\right)\\
&-\sum_i\left(\Delta\,\hat b_i^+\hat c_i^+
+\Delta^*\,\hat b_i\hat c_i\right)
+U\sum_i(\hat n_{b,i}-\hat n_{c,i})^2,
\end{aligned}
```

where $`\hat n_{b,i}=\hat b_i^+\hat b_i`$, $`\hat n_{c,i}=\hat c_i^+\hat c_i`$,
$`t>0`$, and $`U\ge0`$. The trace is grand canonical:

```math
Z=\operatorname{Tr}e^{-\beta\hat H_\mu},\qquad
\hat H_\mu=\hat H-\mu(\hat N_b+\hat N_c),\qquad
\hat N_b=\sum_i\hat n_{b,i},\quad \hat N_c=\sum_i\hat n_{c,i}.
```

The code retains two density-interaction channels,

```math
\hat H_U=\sum_i\left[
U_1(\hat n_{b,i}+\hat n_{c,i})^2
+U_2(\hat n_{b,i}-\hat n_{c,i})^2\right].
```

The main-text benchmark sets $`U_1=0`$ and $`U_2=U`$. The paired scan has $`U=1`$,
$`\mu=-5`$, $`\beta=4`$, and $`\Delta=0,0.05,\ldots,0.30`$; its four observables
form panels (e–h) of the combined benchmark figure. The code reads $`U_1`$,
$`U_2`$, $`\mu`$, and a real $`\Delta`$ from `paramC_sets.txt`; $`t`$ is `RT=1` in
`src/calc_basic.f90`. Both BAFQMC and ED use these same input coefficients.
The interaction is equivalently

```math
(U_1+U_2)(\hat n_{b,i}^2+\hat n_{c,i}^2)
+2(U_1-U_2)\hat n_{b,i}\hat n_{c,i}
```

on each site, as used in the ED construction and energy estimator.

### Manuscript notation and pair phase

For the real pairing coefficient used in the calculations, the implementation
uses the flavor phase convention

```math
\hat b_{i,\rm code}=\hat b_{i,\rm paper},\qquad
\hat c_{i,\rm code}=-\hat c_{i,\rm paper}.
```

In these code operators, the same physical pair term is

```math
\hat H_{\Delta,\rm code}
=+\Delta\sum_i\left(
\hat b_{i,\rm code}^+\hat c_{i,\rm code}^+
+\hat b_{i,\rm code}\hat c_{i,\rm code}\right).
```

**All implementation derivations below use the code operators.** Their hats
and `code` subscripts are omitted for readability, including in the Nambu
basis, Green-function block table, HS scalar factor, and ED identities.
The positive pair coefficient in those formulas is the code representation
of the paper's negative pair term.

Hopping and density interactions are unchanged by this phase rotation.
The four plotted observables $`\rho`$, $`-E`$, $`S_{\rm SF}(K)`$, and
$`S_{\rm DW}(K)`$ therefore use the stored values directly. The anomalous
amplitude changes sign:

```math
P_{\rm paper}=\frac1{N_s}\sum_i
\left\langle\hat b_i\hat c_i+\hat b_i^+\hat c_i^+\right\rangle_{\rm paper}
=-\texttt{pair\_equal},\qquad
e_\Delta=-\Delta P_{\rm paper}=\Delta\,\texttt{pair\_equal}.
```

The total pair energy is $`E_\Delta=N_s e_\Delta`$. The physical energy excludes
$`-\mu(\hat N_b+\hat N_c)`$ throughout both implementations.

For this $`U_1=0`$ model, the condition $`\mu<-3t-|\Delta|`$ guarantees a finite
trace for every real auxiliary-field configuration, as proved in the SM
subsection "Convergence throughout the auxiliary-field domain". All seven
paired benchmark points satisfy it. Only the relative-density HS field is
coupled to the bosons, so the total-density normal-ordering factor is inactive.

## Triangular Lattice Vectors And Bonds

The Hamiltonian's triangular lattice is defined by integer cell topology. The spatial cell index is $`(x,y)`$, with periodic boundaries, and each site stores three forward bonds:

```math
(x,y)\rightarrow(x+1,y),\qquad
(x,y)\rightarrow(x,y+1),\qquad
(x,y)\rightarrow(x-1,y+1).
```

This is a standard triangular nearest-neighbor graph. The pairing and number-conserving solvers
use the same real-space and reciprocal-space convention:

```math
\mathbf a_1=(1,0),\qquad
\mathbf a_2=\left(\frac12,\frac{\sqrt3}{2}\right),
```

```math
\mathbf b_1=2\pi\left(1,-\frac{1}{\sqrt3}\right),\qquad
\mathbf b_2=2\pi\left(0,\frac{2}{\sqrt3}\right),
\qquad \mathbf b_i\cdot\mathbf a_j=2\pi\delta_{ij}.
```

The three forward displacements are $`\mathbf a_1`$, $`\mathbf a_2`$, and
$`-\mathbf a_1+\mathbf a_2`$, all with length $`1`$. The kinetic estimator
uses the integer bond table and adds both directions explicitly. In
`src/lattice.f90`, the zero-based displacement labels are used to populate

```math
\texttt{aimj\_v}(x,y)=(x-1)\mathbf a_1+(y-1)\mathbf a_2,
```

```math
\texttt{xk\_v}(m,n)=
\frac{m}{L_x}\mathbf b_1+
\frac{n}{L_y}\mathbf b_2,
```

with zero-based momentum indices $`m=x-1,n=y-1`$, and `k_dot_r` is the
ordinary dot product. For $`3\times3`$, the code-facing $`K`$ index
`Latt%inv_cell_list(2*Nlx/3+1, Nly/3+1)` represents

```math
K=\frac{2}{3}\mathbf b_1+\frac{1}{3}\mathbf b_2
=\left(\frac{4\pi}{3},0\right),
```

which is one of the triangular-lattice single-particle band minima for the
current hopping sign convention.

## Continuous Hubbard-Stratonovich Fields

The Trotter step uses two continuous Gaussian fields $`\phi_{1,i,\tau}`$ and $`\phi_{2,i,\tau}`$. For the intended signs $`U_1\le0`$ and $`U_2\ge0`$,

```math
e^{-\Delta\tau U_1(n_b+n_c)^2}
\propto \int d\phi_1\, e^{-\phi_1^2/2}
e^{\alpha_1\phi_1(n_b+n_c)},
\qquad
\alpha_1=\sqrt{-2U_1\Delta\tau},
```

```math
e^{-\Delta\tau U_2(n_b-n_c)^2}
\propto \int d\phi_2\, e^{-\phi_2^2/2}
e^{i\alpha_2\phi_2(n_b-n_c)},
\qquad
\alpha_2=\sqrt{2U_2\Delta\tau}.
```

The auxiliary fields live on spatial sites, not Nambu sectors:
`Conf%phi_list(ns, i_site, ntau)`, where `ns=1` is the $`U_1`$ channel and `ns=2` is the $`U_2`$ channel.
The paper defines $`\alpha_2`$ as real; `OperatorHubbard%alpha` stores the complex coefficient $`i\alpha_2`$ for this channel. At zero coupling the corresponding coefficient vanishes.

## Nambu Basis Convention

The active pairing implementation uses the four-component bosonic Nambu basis

```math
\chi=
\begin{pmatrix}
b\\ c\\ b^+\\ c^+
\end{pmatrix},
\qquad
H=-\frac12\chi^T\Omega A\chi+\mathrm{const.},
```

with

```math
\Omega=
\begin{pmatrix}
0&0&I&0\\
0&0&0&I\\
-I&0&0&0\\
0&-I&0&0
\end{pmatrix}.
```

For the one-body triangular hopping, chemical potential, and onsite pairing
piece, the commutator matrix has the block form

```math
A=
\begin{pmatrix}
T&0&0&\Delta I\\
0&T&\Delta I&0\\
0&-\Delta I&-T&0\\
-\Delta I&0&0&-T
\end{pmatrix}.
```

Here $`T`$ is the normal one-body matrix, including the hopping and chemical
potential terms in the code convention. This is the matrix convention used by
`src/non_interact.f90` and by the HS update matrices. It is not the reduced
two-component convention

```math
H=\frac12\tilde\chi^+ \tilde A\tilde\chi,\qquad
\tilde\chi=(b,c^+)^T.
```

That reduced convention may be a future optimization, but it is not the current
implementation. Do not switch formulas between the two conventions without
rederiving the determinant weight, Green matrix, and observable block table.

Although the physical Hamiltonian is Hermitian, the commutator matrix $`A`$ in
this bosonic Nambu representation is generally non-Hermitian when
$`\Delta\neq 0`$. The one-body propagator must therefore be built as a general
matrix exponential of $`A`$, not by a Hermitian diagonalization of $`A`$. For a
single normal mode with normal coefficient $`h`$, the paired Nambu block has
eigenvalues

```math
\sqrt{h^2-\Delta^2},\qquad -\sqrt{h^2-\Delta^2},
```

in the stable regime. A Hermitian diagonalization would instead produce the
wrong $`\sqrt{h^2+\Delta^2}`$ scale and gives a visibly biased finite-$`\Delta`$
benchmark even in the noninteracting limit. `src/non_interact.f90` therefore
uses LAPACK `zgeev` to diagonalize the general complex matrix before forming
$`\exp(\mp\Delta\tau A)`$.

## HS Constant In The Nambu Convention

Let $`\lambda_1`$ and $`\lambda_2`$ denote the actual complex coefficients
that multiply the normal-ordered density channels on one space-time site:

```math
\lambda_1=\alpha_1\phi_1,\qquad
\lambda_2=i\alpha_2\phi_2,\qquad
\alpha_2=\sqrt{2U_2\Delta\tau}.
```

In the Fortran code, `OperatorHubbard%alpha` stores $`i\alpha_2`$, so its
product with the real field is $`\lambda_2`$. The HS density term satisfies

```math
\lambda_1(b^+ b+c^+ c)
+\lambda_2(b^+ b-c^+ c)
=
-\frac12\chi^T\Omega
\begin{pmatrix}
\lambda_1+\lambda_2&0&0&0\\
0&\lambda_1-\lambda_2&0&0\\
0&0&-(\lambda_1+\lambda_2)&0\\
0&0&0&-(\lambda_1-\lambda_2)
\end{pmatrix}
\chi
-\lambda_1 .
```

The final $`-\lambda_1`$ is the bosonic commutator c-number. It does not appear
as a matrix element in the Nambu propagator, but it must appear in the sampling
weight. Therefore a local $`U_1`$ update contributes

```math
\frac{e^{-\alpha_1\phi_1'}}{e^{-\alpha_1\phi_1}}
=
\exp[-\alpha_1(\phi_1'-\phi_1)]
```

to `ratio_constant`. In `src/operator_Hubbard.f90` this is implemented as
`expalpha_old / expalpha_new` for `IUflag == 1`. The $`U_2`$ channel has no
scalar c-number because its two flavor coefficients sum to zero.

## Green Matrix And Local Determinant Factor

For one full imaginary-time product, define

```math
M=I-e^{A_1}e^{A_2}\cdots e^{A_N},
\qquad
G=M^{-1}\Omega.
```

The physical Green matrix is $`G`$. The quantity propagated and updated in the
Fortran code is instead

```math
\tilde G\equiv G\Omega^{-1}=M^{-1}.
```

Despite the historical variable name, `Prop%Gr` stores $`\tilde G`$, not the
physical $`G`$. The local Metropolis probability is

```math
p(\phi\to\phi')=
\min\left[
1,\left|
\texttt{ratio\_constant}\,
\texttt{ratio\_gaussian}\,
\texttt{ratio\_det}
\right|\right],
```

where

```math
\texttt{ratio\_constant}
=
\frac{e^{-\alpha_1\phi_1'}}{e^{-\alpha_1\phi_1}},
\qquad
\texttt{ratio\_gaussian}
=
\frac{e^{-\phi'^2/2}}{e^{-\phi^2/2}}.
```

The determinant factor is the current Nambu determinant contribution. In
`src/localU.f90`, the local update builds a doubled-sector $`4\times4`$ block
$`P`$ and evaluates

```math
\texttt{ratio\_det}=\det(P)^{-1/2}.
```

The exponent $`-1/2`$ is part of the present doubled Nambu convention. Do not
change it to $`-1`$ or $`-1/4`$ without a full rederivation and benchmark
update.

## Program Green-Function Block Table

The raw `Prop%Gr` matrix is $`\tilde G`$. In the sector order
$`(b,c,b^+,c^+)`$, its blocks represent

```math
\tilde G=
\left\langle
\begin{pmatrix}
bb^+ & bc^+ & -bb & -bc\\
cb^+ & cc^+ & -cb & -cc\\
b^+ b^+ & b^+ c^+ & -b^+ b & -b^+ c\\
c^+ b^+ & c^+ c^+ & -c^+ b & -c^+ c
\end{pmatrix}
\right\rangle .
```

`src/obser_equal.f90` converts these raw blocks into physical contractions by
undoing the explicit minus signs in the table. For example,
$`G_{b,c}=-\texttt{Prop\%Gr(b,cdag)}`$ and
$`G_{b^+,b}=-\texttt{Prop\%Gr(bdag,bdag)}`$.

Define

```math
n_b(i,j)=G_{b^+,b}(i,j),\qquad
n_c(i,j)=G_{c^+,c}(i,j).
```

These are the contractions used directly in physical normal-ordered
observables. In particular, cross-flavor products use the two sectors that are
present in the explicit Nambu matrix, $`n_b n_c`$. They should not be replaced
by an ad hoc conjugated product when evaluating `doubleOcc`.

## Observable Estimator Formulas

The scalar equal-time estimators are accumulated per configuration and then averaged over observations, MPI ranks, and bins. The density and number estimators are

```math
N_b=\sum_i \mathrm{Re}\,n_b(i,i),\qquad
N_c=\sum_i \mathrm{Re}\,n_c(i,i),
```

```math
\rho_b=\frac{N_b}{N_s},\qquad
\rho_c=\frac{N_c}{N_s},\qquad
\rho=\rho_b+\rho_c.
```

The kinetic output is per spatial site and includes the hopping coefficient:

```math
K_{\mathrm{DQMC}} =
\frac{t}{N_s}\sum_{i}\sum_{\delta}
\mathrm{Re}\left[
n_b(i,j_\delta)+n_b(j_\delta,i)
+n_c(i,j_\delta)+n_c(j_\delta,i)
\right].
```

The Wick contractions used by `D_bb`, `D_cc`, and `D_bc` are

```math
D_{bb}(i,j)=
n_b(i,i)n_b(j,j)
+G_{b^+,b}(i,j)G_{b,b^+}(i,j)
+G_{b^+,b^+}(i,j)G_{b,b}(i,j),
```

```math
D_{cc}(i,j)=
n_c(i,i)n_c(j,j)
+G_{c^+,c}(i,j)G_{c,c^+}(i,j)
+G_{c^+,c^+}(i,j)G_{c,c}(i,j),
```

```math
D_{bc}(i,j)=
n_b(i,i)n_c(j,j)
+G_{b^+,c^+}(i,j)G_{b,c}(i,j)
+G_{b^+,c}(i,j)G_{b,c^+}(i,j).
```

The occupancy and pairing outputs are

```math
\texttt{doubleOcc}=
\frac{1}{N_s}\sum_i \mathrm{Re}\,D_{bc}(i,i),
```

```math
\texttt{squareOcc}=
\frac{1}{N_s}\sum_i \mathrm{Re}\,[n_b(i,i)^2+n_c(i,i)^2],
```

```math
\texttt{local\_numsquare}=
\frac{1}{N_s}\sum_i \mathrm{Re}\,[D_{bb}(i,i)+D_{cc}(i,i)],
```

```math
\texttt{numsquare\_up}=
\sum_{ij}\mathrm{Re}\,D_{bb}(i,j),\qquad
\texttt{numsquare\_do}=
\sum_{ij}\mathrm{Re}\,D_{cc}(i,j),
```

```math
\texttt{pair\_equal}=
\frac{1}{N_s}\sum_i\mathrm{Re}\,
\left[G_{b,c}(i,i)+G_{b^+,c^+}(i,i)\right].
```

`squareOcc` preserves the density-product estimator above. In the current
model, conservation of $`N_b-N_c`$ removes same-flavor anomalous contractions,
so its physical meaning is half the per-site normal-ordered same-flavor
onsite pair density. `local_numsquare` includes the full Wick contractions
for $`N_s^{-1}\sum_i\langle n_{b,i}^2+n_{c,i}^2\rangle`$ and is compared to ED.
See the [observable reference](../../observables.md) for the current output
inventory and the corresponding distinction for more general paired models.

### `Delta=0` And Corrected No-Pairing `doubleOcc`

The no-pairing reference code samples only the `b` determinant and reconstructs
the `c` sector by complex conjugation. In
`src/number_conserving/src/obser_equal.f90`, the relevant definitions are

```math
\texttt{Grupc}=G_b^T-I,\qquad
\texttt{Grdo}=G_b^\ast,\qquad
\texttt{Grdoc}=\texttt{Grdo}^T-I.
```

Therefore

```math
\texttt{Grdoc}
=(G_b^\ast)^T-I
=(G_b^T-I)^\ast
=\texttt{Grupc}^\ast,
```

and the corrected no-pairing output named `doubleOcc` is, configuration by
configuration, the physical cross-flavor estimator

```math
\texttt{doubleOcc}_{\mathrm{no\ pairing}}
=\frac{1}{N_s}\sum_i\mathrm{Re}\,
[\texttt{Grupc}(i,i)\texttt{Grupc}(i,i)^\ast]
=\frac{1}{N_s}\sum_i|\texttt{Grupc}(i,i)|^2.
```

This is the same object as the explicit-Nambu `doubleOcc` at $`\Delta=0`$,
where the anomalous blocks vanish and the `c` normal sector is the complex
conjugate of the `b` sector. Therefore Benchmark A compares `doubleOcc`
directly between the two repositories:

```math
\texttt{doubleOcc}_{\mathrm{pairing},\,\Delta=0}
=\frac{1}{N_s}\sum_i\mathrm{Re}\,D_{bc}(i,i),
```

against the corrected no-pairing `doubleOcc`. Older no-pairing output produced
with `Grdoc = dconjg(transpose(Grdo)) - ZKRON` should not be used as a
reference, because that bug instead gave $`\texttt{Grdoc}=\texttt{Grupc}`$ and
turned `doubleOcc` into $`\mathrm{Re}[\texttt{Grupc}^2]`$.

For density correlations, `den_upup_sub11`, `den_dodo_sub11`, and `den_updo` are the $`k=0`$ Fourier components of $`D_{bb}`$, $`D_{cc}`$, and $`D_{bc}`$ after the real-space accumulation by separation $`r_i-r_j`$.

## ED Finite-Temperature Trace

Because $`H_\Delta`$ changes total particle number, finite-$`\Delta`$ ED uses the grand-canonical Hamiltonian directly rather than fixed-$`N_b,N_c`$ blocks. The general ED script builds a no-symmetry `boson_basis_general` basis over $`2N_s`$ boson modes with local cutoff `nmax` and optional total cutoff `ncut`.

For the full eigensystem $`H_\mu|\alpha\rangle=E_\alpha|\alpha\rangle`$,

```math
Z=\sum_\alpha e^{-\beta E_\alpha},\qquad
\langle O\rangle=
\frac{1}{Z}\sum_\alpha e^{-\beta E_\alpha}
\langle\alpha|O|\alpha\rangle.
```

The implementation shifts energies by $`E_0=\min_\alpha E_\alpha`$ for numerical stability and reports

```math
F=E_0-\frac{1}{\beta}\log\sum_\alpha e^{-\beta(E_\alpha-E_0)}.
```

The ED identity tests check

```math
\frac{\partial F}{\partial \Delta}
=N_s\,\texttt{pair\_equal},
\qquad
\frac{\partial F}{\partial t}
=\frac{\texttt{kinetic\_total}}{t}
\quad (t\ne0).
```

Both solvers use the same benchmark energy convention. The
paper-facing plotted energy is the total physical/internal energy

```math
E=\langle H_{\rm phys}\rangle,
\qquad
\text{paper panel}=-E .
```

It excludes the chemical-potential term used in the grand-canonical sampling
Hamiltonian. The DQMC and ED JSON files store the per-site value

```math
e=\frac{E}{N_s}=\texttt{energy\_density},
```

so report and plotting scripts must multiply `energy_density` by $`N_s`$
before drawing $`-E`$. The chemical-potential contribution and the diagnostic
grand-canonical energy density remain per-site diagnostics:

```math
e_\mu=-\mu\rho,\qquad e_{\rm grand}=e+e_\mu=e-\mu\rho.
```

## DQMC/ED Normalization Conversions

| output file | DQMC estimator | DQMC normalization | ED observable | comparison operation | Benchmark A | Benchmark B |
| --- | --- | --- | --- | --- | --- | --- |
| `num_up` | $`\sum_i \mathrm{Re}\,n_b(i,i)`$ | total | $`N_b`$ | `mean` | yes | yes |
| `num_do` | $`\sum_i \mathrm{Re}\,n_c(i,i)`$ | total | $`N_c`$ | `mean` | yes | yes |
| `density_up` | `num_up / Nsite` | per site | $`N_b/N_s`$ | `mean` | yes | yes |
| `density_do` | `num_do / Nsite` | per site | $`N_c/N_s`$ | `mean` | yes | yes |
| `density` | `density_up + density_do` | per site | $`(N_b+N_c)/N_s`$ | `mean` | yes | yes |
| `density_total` | alias of `density` | per site | `density_total` | `mean` | yes | yes |
| `density_site_total` | $`\rho_i=\mathrm{Re}[G_{b^+ b}(i,i)+G_{c^+ c}(i,i)]`$ | site vector | IPR postprocessing | fixed-width row | no | yes |
| `kinetic` | hopping expectation including `RT` | per site | ED JSON `kinetic_total` = total $`H_t`$ | `mean_times_nsite` | yes | yes |
| `doubleOcc` | $`N_s^{-1}\sum_i \mathrm{Re}\,D_{bc}(i,i)`$ | per site | $`N_s^{-1}\sum_i n_{b,i}n_{c,i}`$ | `mean` | yes | yes |
| `squareOcc` | $`N_s^{-1}\sum_i \mathrm{Re}[n_b(i,i)^2+n_c(i,i)^2]`$ | per site | half normal-ordered same-flavor onsite pair density in the current model | `mean` | yes | no |
| `local_numsquare` | $`N_s^{-1}\sum_i \mathrm{Re}[D_{bb}(i,i)+D_{cc}(i,i)]`$ | per site | $`N_s^{-1}\sum_i(n_{b,i}^2+n_{c,i}^2)`$ | `mean` | yes when reference exists | yes |
| `onsite_n2_up` | $`\sum_i\mathrm{Re}\,D_{bb}(i,i)`$ | total | `onsite_n2_up` | `mean` | yes when reference exists | yes |
| `onsite_n2_do` | $`\sum_i\mathrm{Re}\,D_{cc}(i,i)`$ | total | `onsite_n2_do` | `mean` | yes when reference exists | yes |
| `numsquare_up` | $`\sum_{i,j}\mathrm{Re}\,D_{bb}(i,j)`$ | total | $`N_b^2`$ | `mean` | yes | yes |
| `numsquare_do` | $`\sum_{i,j}\mathrm{Re}\,D_{cc}(i,j)`$ | total | $`N_c^2`$ | `mean` | yes | yes |
| `pair_equal` | $`N_s^{-1}\sum_i\mathrm{Re}\langle b_i c_i+b_i^+ c_i^+\rangle`$ | per site | $`N_s^{-1}\partial F/\partial\Delta`$ | `mean` | zero check | yes |
| `interaction_energy_density` | $`(U_1+U_2)(M_b^{(2)}+M_c^{(2)})/N_s+2(U_1-U_2)\texttt{doubleOcc}`$ | per site | `interaction_energy_density` | `mean` | no | yes |
| `pairing_energy_density` | $`\Delta\,\texttt{pair\_equal}`$ | per site | `pairing_energy_density` | `mean` | zero check | yes |
| `energy_density` | $`e=E/N_s=e_t+e_U+e_\Delta`$, excluding $`-\mu N`$; paper plots use $`-E=-N_s e`$ | per site storage | `energy_density` | `mean` | no | yes |
| `chemical_energy_density` | $`-\mu\rho`$ | per site | `chemical_energy_density` | `mean` | diagnostic | diagnostic |
| `grand_energy_density` | `energy_density + chemical_energy_density` | per site | `grand_energy_density` | `mean` | diagnostic | diagnostic |
| `sf_K` | $`N_s^{-2}\sum_{ij}e^{iK\cdot(r_i-r_j)}\langle b_i^+ b_j+c_i^+ c_j\rangle`$ | complex scalar | `S_SF_K` | real/imag block stats | no | yes |
| `dw_K` | $`N_s^{-2}\sum_{ij}e^{iK\cdot(r_i-r_j)}[D_{bb}(i,j)+D_{cc}(i,j)+D_{bc}(i,j)+D_{cb}(i,j)]`$ | complex scalar | `S_DW_K` | real/imag block stats | no | yes |
| `psf_Gamma` | $`N_s^{-2}\sum_{ij}\langle b_i^+ c_i^+ c_j b_j\rangle`$ | complex scalar | `S_PSF_Gamma` | real/imag block stats | no | yes |

Here

```math
M_b^{(2)}=\sum_i\langle n_{b,i}^2\rangle,\qquad
M_c^{(2)}=\sum_i\langle n_{c,i}^2\rangle.
```

For $`K`$-point observables the pairing code uses the same $`K`$ index and
phase convention as the number-conserving solver when both lattice dimensions are
multiples of three. ED entrypoints reject incompatible $`K`$ labels, and
campaign analysis must mark `sf_K` and `dw_K` unavailable rather than treating
placeholders or zeros as physics data.
At compatible sizes `sf_K` and `dw_K` are standard Hermitian structure-factor
channels:

```math
S_{\rm SF}(K)=\frac{1}{N_s}\langle b_K^+b_K+c_K^+c_K\rangle,
\qquad
S_{\rm DW}(K)=\frac{1}{N_s^2}\langle n_K n_{-K}\rangle ,
\quad n_K=\sum_i e^{-iK\cdot r_i}(n_{b,i}+n_{c,i}).
```

Their exact thermal averages are real and non-negative. DQMC still writes two
columns because the estimator is accumulated as a complex Fourier sum; the
real part is the comparison channel and a statistically significant imaginary
part marks the observable unreliable. The same imaginary-part reliability rule
applies to the $`\Gamma`$-point `psf_Gamma`.
The IPR computed from `density_site_total` is a DQMC spatial-uniformity
diagnostic. It is useful for spotting broken site indexing or severe sampling
inhomogeneity, but it should not be used as a strict ED agreement gate in
low-density pilots because the ratio
$`\sum_i \rho_i^2/(\sum_i\rho_i)^2`$ has a finite-sample upward bias.

The kinetic naming contract is:

```text
DQMC file `kinetic`: per spatial site.
ED JSON `kinetic_total`: total hopping expectation including the coefficient `t`.
DQMC-vs-ED and DQMC-vs-DQMC comparison operation: `mean_times_nsite`.
```
