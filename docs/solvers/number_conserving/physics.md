# Physics And Code Map

Implementation paths and commands in this guide are relative to
`src/number_conserving/` from the repository root.

## Model

The code simulates the two-flavor Bose-Hubbard model on a periodic triangular
lattice. In the paper's creation-operator notation, its physical Hamiltonian is

```math
\begin{aligned}
\hat H={}&t\sum_{\langle ij\rangle}
\left(\hat b_i^+\hat b_j+\hat b_j^+\hat b_i
+\hat c_i^+\hat c_j+\hat c_j^+\hat c_i\right)\\
&+U_1\sum_i(\hat n_{b,i}+\hat n_{c,i})^2
+U_2\sum_i(\hat n_{b,i}-\hat n_{c,i})^2,
\end{aligned}
```

where $`\hat n_{b,i}=\hat b_i^+\hat b_i`$, $`\hat n_{c,i}=\hat c_i^+\hat c_i`$, and $`N_s=L_xL_y`$ for the
single orbital per cell. The code names this site count `Lq`.

The current convention is:

- $`t=1>0`$
- $`U_1\le0`$
- $`U_2\ge0`$
- flavors `b` and `c` are also called `up` and `down`

The main-text benchmark sets $`U_1=0`$ and calls the relative-density coupling
$`U=U_2`$. Its eight interaction-scan points form panels (a–d) of the combined
benchmark figure. The Supplemental Material retains the two-channel
Hamiltonian above and scans seven $`U_1`$ values at $`U_2=1`$. The Monte Carlo trace is

```math
Z=\operatorname{Tr}e^{-\beta\hat H_\mu},\qquad
\hat H_\mu=\hat H-\mu(\hat N_b+\hat N_c),\qquad
\hat N_b=\sum_i\hat n_{b,i},\quad \hat N_c=\sum_i\hat n_{c,i}.
```

The reported physical energy $`E=\langle\hat H\rangle`$ excludes the chemical-potential term.
Operator hats are omitted in the implementation derivations below for
readability. The [observable reference](../../observables.md) lists every
active physical output and its normalization.

For the main-text model, $`\mu<-3t`$ ensures a finite trace for every
auxiliary-field configuration. The relative-density HS factors below are
unitary, and the quadratic propagator damps high occupations. The SM proves
the resulting uniform trace bound, together with its paired extension
$`\mu<-3t-|\Delta|`$. This proof concerns $`U_1=0`$; the supplemental
attractive-density scan uses its stated finite-occupation references.

## Hubbard-Stratonovich Fields

The Trotter decomposition uses two real Gaussian fields on each site and time
slice. With the paper's real coefficients $`\alpha_1,\alpha_2\ge0`$,

```math
e^{-\Delta\tau U_1(n_b+n_c)^2}
=\int\frac{d\phi_1}{\sqrt{2\pi}}\,
e^{-\phi_1^2/2+\alpha_1\phi_1(n_b+n_c)},
\qquad \alpha_1=\sqrt{-2U_1\Delta\tau},
```

```math
e^{-\Delta\tau U_2(n_b-n_c)^2}
=\int\frac{d\phi_2}{\sqrt{2\pi}}\,
e^{-\phi_2^2/2+i\alpha_2\phi_2(n_b-n_c)},
\qquad \alpha_2=\sqrt{2U_2\Delta\tau}.
```

The code stores the total-density coefficient as the real number $`\alpha_1`$
and the relative-density coefficient as the imaginary number $`i\alpha_2`$;
`Dtau` is $`\Delta\tau`$.

At zero interaction strength the corresponding coupling vanishes. Thus the
main-text scan has only the relative-density phase field coupled to the bosons;
at $`U=0`$ both couplings vanish.

After decoupling, the two flavor Hamiltonians are complex conjugates. The code samples the `b` flavor explicitly. The `c` flavor Green matrix is reconstructed with complex conjugation.

## Green Function Convention

For a fixed auxiliary-field configuration, the equal-time Green matrix is

```math
(G_b)_{ij}(\tau)=\langle b_i b_j^+\rangle_\phi,\qquad
G_b(\tau)=[I-\mathcal B_b(\tau)]^{-1},\qquad
G_c(\tau)=G_b(\tau)^*.
```

Here $`\mathcal B_b(\tau)`$ is the cyclic imaginary-time propagator product at
the measurement slice. Bosonic commutation gives the normal-ordered matrix

```math
\langle b_i^+b_j\rangle_\phi=(G_b)_{ji}-\delta_{ij}.
```

In code, `Prop%Gr` stores the `b` flavor matrix. `ObserEqual_mod` constructs:

```text
Grup  = Prop%Gr
Grupc = transpose(Grup) - ZKRON
Grdo  = dconjg(Prop%Gr)
Grdoc = transpose(Grdo) - ZKRON
```

## Continuous Pole Diagnostics

The continuous-HS run writes configuration-level pole diagnostics alongside the scalar observables:

| DQMC output file | Meaning |
| --- | --- |
| `pole_z` | Re/Im pairs for $`z_a=1/\lambda_a(G)`$ |
| `pole_distance` | $`\min_a\lvert z_a\rvert`$ |
| `pole_x` | `-log10(pole_distance)` |
| `green_spectral_radius` | $`\max_a\lvert\lambda_a(G)\rvert`$ |
| `green_smax` | largest singular value of $`G`$ |
| `log_weight` | $`\log P_{\rm HS}+2\sum_a\log s_a(G)`$ |

Here $`\lambda_a(G)`$ are eigenvalues of the $`b`$-flavor Green matrix and
$`s_a(G)`$ are its singular values. For the continuous fields,

```math
\log P_{\rm HS}=-\frac12\sum_{a=1}^{2}\sum_{i,\ell}\phi_{a,i,\ell}^{\,2},
```

with normalization constants omitted from `log_weight`.

These diagnostics are sampled once per bin per MPI rank. For `MPI_NP > 1`, each rank contributes one independent configuration sample per bin, and rank 0 gathers and writes `ISIZE` rows per bin in rank order. The rows are not rank averaged. Use `MPI_NP=1` when following a single Markov-chain time series or pole-spike trace.

Pole diagnostics are configuration diagnostics. They do not change the density, number, kinetic, or occupancy observable normalizations below.

## Physical Symbols To Code Variables

| Physical symbol | Code variable or location |
| --- | --- |
| `Lx`, `Ly` | `Nlx`, `Nly` in `src/calc_basic.f90`, read from `paramC_sets.txt` |
| `Lq = Lx * Ly` | `Lq` in `CalcBasic` |
| `beta` | `Beta` in `CalcBasic` |
| `Delta tau` | `Dtau = Beta / Ltrot` in `Params_set` |
| `Ltrot` | `Ltrot` in `CalcBasic`, read from `paramC_sets.txt` |
| `t` | `RT` in `CalcBasic`, currently set to `1.d0` in `Params_set` |
| `U1`, `U2` | `RU1`, `RU2` in `CalcBasic`, read from `paramC_sets.txt` |
| main-text `U` | `RU2` with `RU1=0` |
| `mu` | `mu` in `CalcBasic`, read from `paramC_sets.txt` |
| auxiliary field flavor index | `ns = 1` for `U1`, `ns = 2` for `U2` |
| auxiliary fields | `Conf%phi_list(ns, ii, nt)` in `src/fields.f90` |
| triangular lattice nearest-neighbor bonds | `Latt%L_bonds(ii, nb)` in `src/lattice.f90` |
| space-time bonds | `Latt%LT_bonds(iit, nb)` in `src/lattice.f90` |
| `b` flavor Green matrix | `Prop%Gr` |
| `c` flavor Green matrix | `dconjg(Prop%Gr)` |
| local update shift | `shiftLoc` in `CalcBasic`, read from `paramC_sets.txt` |
| warm-up shift | `shiftWarm(1:2)` in `CalcBasic`, read from `paramC_sets.txt` |

## Observable Map

| DQMC output file | Fortran field | Meaning |
| --- | --- | --- |
| `density_up` | `Obs%density_up` | per-site density for the `b` flavor |
| `density_do` | `Obs%density_do` | per-site density for the `c` flavor |
| `num_up` | `Obs%num_up` | total `b`-flavor particle number |
| `num_do` | `Obs%num_do` | total `c`-flavor particle number |
| `kinetic` | `Obs%kinetic` | kinetic observable accumulated with division by `Lq` in `Obs_equal_calc` |
| `doubleOcc` | `Obs%doubleOcc` | per-site cross-flavor onsite density product, `<sum_i n_b,i n_c,i> / Lq` |
| `squareOcc` | `Obs%squareOcc` | per-site half normal-ordered same-flavor onsite pair, `0.5 * <sum_i,s n_s,i (n_s,i - 1)> / Lq` |
| `numsquare_up` | `Obs%numsquare_up` | total `b`-flavor number-square estimator |
| `numsquare_do` | `Obs%numsquare_do` | total `c`-flavor number-square estimator |
| `den_upup_sub11` | `Obs%den_corr_up` after Fourier transform | `b-b` density correlation for the single orbital case |
| `den_dodo_sub11` | `Obs%den_corr_do` after Fourier transform | `c-c` density correlation for the single orbital case |
| `den_updo` | `Obs%den_corr_updo` after Fourier transform | cross-flavor density correlation |

## 3x3 Campaign Observables

The 3x3 BAFQMC-vs-ED campaign is a separate research workflow from the
default 3x2 regression benchmark. It adds equal-time DQMC files and a 3x3
QuSpin ED reference for observables needed in the paper-scale comparison.
The K point is available only when both `Lx` and `Ly` are multiples of 3.
For the 3x3 triangular lattice, the code-facing K index is
`Latt%inv_cell_list(2*Nlx/3 + 1, Nly/3 + 1)`.

The lattice uses the triangular primitive and reciprocal vectors

```math
\mathbf a_1=(1,0),\qquad
\mathbf a_2=\left(\frac12,\frac{\sqrt3}{2}\right),
```

```math
\mathbf b_1=2\pi\left(1,-\frac1{\sqrt3}\right),\qquad
\mathbf b_2=2\pi\left(0,\frac2{\sqrt3}\right),\qquad
\mathbf b_i\cdot\mathbf a_j=2\pi\delta_{ij}.
```

Each site has three forward bonds, with periodic cell indices:

```math
(x,y)\longrightarrow(x+1,y),\quad(x,y+1),\quad(x-1,y+1).
```

The Hamiltonian includes their Hermitian conjugates explicitly. For $`3\times3`$,
the code-facing K index above has integer momentum coordinate $`(2,1)`$ and represents

```math
\mathbf K=\frac23\mathbf b_1+\frac13\mathbf b_2
=\left(\frac{4\pi}{3},0\right).
```
Older campaign data generated with the equivalent primitive-cell relabeling
$`\mathbf a_2'=\mathbf a_2-\mathbf a_1`$ describe the same triangular-lattice K channel, so those
benchmark results do not require rerunning solely for this coordinate-label
change.

The campaign uses $`N_s=L_xL_y`$ (`Lq`) and structure-factor normalization
$`N_s^{-2}`$. For $`3\times3`$, $`N_s=9`$ and $`N_s^2=81`$.

| DQMC output file | ED/schema name | Definition and normalization |
| --- | --- | --- |
| `density_total` | `density_total` | $`\rho=(\langle N_b\rangle+\langle N_c\rangle)/N_s`$ |
| `density_site_total` | postprocessed `IPR` input | site-resolved total density profile $`\rho_i=\langle n_{b,i}+n_{c,i}\rangle`$ |
| `energy_density` | `energy_density` | physical energy density without the grand-canonical `-mu N` term: `kinetic + interaction_energy_density` |
| `doubleOcc` | `doubleOcc` | $`N_s^{-1}\sum_i\langle n_{b,i}n_{c,i}\rangle`$ |
| `onsite_n2_up` | `onsite_n2_up` | $`\sum_i\langle n_{b,i}^2\rangle`$; DQMC estimator `sum_i (2*Grupc_ii^2 + Grupc_ii)` |
| `onsite_n2_do` | `onsite_n2_do` | $`\sum_i\langle n_{c,i}^2\rangle`$; DQMC estimator `sum_i (2*Grdoc_ii^2 + Grdoc_ii)` |
| `numsquare_up` | `numsquare_up` | $`\langle N_b^2\rangle`$, where $`N_b=\sum_i n_{b,i}`$ |
| `numsquare_do` | `numsquare_do` | $`\langle N_c^2\rangle`$, where $`N_c=\sum_i n_{c,i}`$ |
| `sf_K` | `S_SF_K` | $`N_s^{-2}\sum_{ij}e^{i\mathbf K\cdot(\mathbf r_i-\mathbf r_j)}\langle b_i^+b_j+c_i^+c_j\rangle`$ |
| `psf_Gamma` | `S_PSF_Gamma` | $`N_s^{-2}\sum_{ij}\langle b_i^+c_i^+c_jb_j\rangle`$ |
| `dw_K` | `S_DW_K` | $`N_s^{-2}\sum_{ij}e^{i\mathbf K\cdot(\mathbf r_i-\mathbf r_j)}\langle(n_{b,i}+n_{c,i})(n_{b,j}+n_{c,j})\rangle`$ |

Equivalently, the plotted structure factors use the standard convention
$`S_O(\mathbf q)=\langle O_{\mathbf q}^+O_{\mathbf q}\rangle/N_s^2`$ with the phase
$`e^{i\mathbf q\cdot(\mathbf r_i-\mathbf r_j)}`$. With this convention they are real and
non-negative in the exact average. The DQMC files are complex because they
store finite-statistics estimators; the imaginary part is a diagnostic and
must be statistically consistent with zero.

`sf_K`, `psf_Gamma`, and `dw_K` are complex-valued DQMC files written as
two columns, `real imag`. For non-K lattices such as the committed 3x2
benchmark inputs, `sf_K` and `dw_K` remain zero rather than stopping the
program; `psf_Gamma` is still computed.
Campaign postprocessing compares the real part with ED and requires the
imaginary part to be statistically consistent with zero before a point is
marked as trusted.

The interaction-energy density is reconstructed consistently in DQMC and ED:

```text
interaction_energy_density =
  (U1 + U2) * (onsite_n2_up + onsite_n2_do) / Lq
  + 2 * (U1 - U2) * doubleOcc

energy_density = kinetic + interaction_energy_density
```

Here DQMC `kinetic` is already per site. ED reports `total_kinetic`, so the
ED conversion is `energy_density = total_kinetic / Lq + interaction_energy_density`.

The inverse participation ratio is computed from the site-resolved total
density distribution, not from same-flavor second moments:

```math
\operatorname{IPR}_\rho=\frac{\sum_i\rho_i^2}{(\sum_i\rho_i)^2},\qquad
\rho_i=\langle n_{b,i}+n_{c,i}\rangle.
```

For a non-negative density profile, $`\operatorname{IPR}_\rho`$ lies between $`1/N_s`$ and $`1`$.
On a uniform $`3\times3`$ profile it is $`1/9`$; if the density is localized on one
site it is $`1`$. If the blocked mean total density is too small, the IPR is
marked unreliable rather than plotted as a trusted comparison.

For the $`3\times3`$ hopping convention, the single-particle minimum is $`\epsilon_{\min}=-3t`$.
Grand-canonical $`U_1=0`$ campaign points use $`\mu<-3t`$; otherwise
the free/equal-flavor growth channel makes the bosonic partition function
diverge. For stable $`U_1=U_2=0`$ points, the campaign can use the exact
grand-canonical free-boson reference instead of a many-body cutoff ED run:

```math
\epsilon_{\mathbf k}=2t\left[
\cos(\mathbf k\cdot\mathbf a_1)+\cos(\mathbf k\cdot\mathbf a_2)
+\cos\bigl(\mathbf k\cdot(\mathbf a_2-\mathbf a_1)\bigr)\right],\qquad
n_{\mathbf k}=\frac1{e^{\beta(\epsilon_{\mathbf k}-\mu)}-1},
```

```math
\begin{aligned}
\rho&=\frac2{N_s}\sum_{\mathbf k}n_{\mathbf k},&
S_{\rm SF}(\mathbf K)&=\frac2{N_s}n_{\mathbf K},\\
S_{\rm DW}(\mathbf K)&=\frac2{N_s^2}\sum_{\mathbf k}n_{\mathbf k+\mathbf K}(1+n_{\mathbf k}),&
S_{\rm PSF}(\Gamma)&=\frac1{N_s^2}\sum_{\mathbf k}n_{\mathbf k}n_{-\mathbf k}.
\end{aligned}
```

The factor of two in $`\rho`$, $`S_{\rm SF}`$, and $`S_{\rm DW}`$ accounts for the two flavors.
Reliable real ED is still preferred when present, but missing or unreliable
free ED is replaced in campaign postprocessing by `ed_status=exact_free_boson`.
Negative `U1` points are finite-window/cutoff comparisons, not fully
converged grand-canonical ED results unless the cutoff diagnostics are
explicitly satisfied. In the campaign postprocessing, a negative-`U1`
finite-window ED point is allowed into trusted curves only when the total
density is at most `5e-2` and at least six nonzero total-particle shells have
been accumulated; earlier `low_density_cutoff_accepted` checkpoints remain
diagnostic rather than publication references.

## Benchmark Normalization

The ED reference script loops over two-species fixed-`NE` blocks and reports:

- total particle number `NE = NE_b + NE_c`
- total kinetic expectation value from both flavor hopping layers
- `doubleOcc`, `squareOcc`, `numsquare_up`, and `numsquare_do`

The DQMC code writes:

- `num_up` for the `b` flavor
- `num_do` for the `c` flavor
- `kinetic` divided by `Lq`
- `doubleOcc` and `squareOcc` with the per-site normalizations in the observable map
- `numsquare_up` and `numsquare_do` as total flavor number-square estimators

Therefore the benchmark comparison should use:

```text
total_NE_DQMC = last(num_up) + last(num_do)
total_kinetic_DQMC = last(kinetic) * Lq
doubleOcc_DQMC = last(doubleOcc)
squareOcc_DQMC = last(squareOcc)
numsquare_up_DQMC = last(numsquare_up)
numsquare_do_DQMC = last(numsquare_do)
```

For real Monte Carlo runs, use sample means rather than the last bin:

```text
total_NE_DQMC = mean(num_up) + mean(num_do)
total_kinetic_DQMC = mean(kinetic) * Lq
doubleOcc_DQMC = mean(doubleOcc)
squareOcc_DQMC = mean(squareOcc)
numsquare_up_DQMC = mean(numsquare_up)
numsquare_do_DQMC = mean(numsquare_do)
```

The live DQMC benchmark estimates the uncertainty of these means by blocking the time series. The reported `stderr` is the standard error of the mean computed from block means, not the standard deviation of raw per-bin samples.

Changing this conversion is a physics-level change and must be reflected in `benchmarks/references/*.json`, `benchmarks/dqmc_references/*.json`, and `benchmarks/README.md`.
