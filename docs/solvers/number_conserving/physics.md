# Physics And Code Map

Implementation paths and commands in this guide are relative to
`src/number_conserving/` from the repository root.

## Model

The code simulates the two-flavor Bose-Hubbard model on a triangular lattice:

```text
H = t sum_<ij> (b_i^dagger b_j + c_i^dagger c_j + h.c.)
  + U1 sum_i (n_{b,i} + n_{c,i})^2
  + U2 sum_i (n_{b,i} - n_{c,i})^2
```

The current convention is:

- `t = 1 > 0`
- `U1 <= 0`
- `U2 >= 0`
- flavors `b` and `c` are also called `up` and `down`

The main-text benchmark sets `U1=0` and calls the relative-density coupling
`U=U2`. Its eight interaction-scan points form panels (a–d) of the combined
benchmark figure. The Supplemental Material retains the two-channel
Hamiltonian above and scans seven `U1` values at `U2=1`. The Monte Carlo trace
uses `H - mu*(N_b+N_c)`; the reported physical energy excludes the chemical
potential term.

For the main-text model, `mu < -3*t` ensures a finite trace for every
auxiliary-field configuration. The relative-density HS factors below are
unitary, and the quadratic propagator damps high occupations. The SM proves
the resulting uniform trace bound, together with its paired extension
`mu < -3*t - abs(Delta)`. This proof concerns `U1=0`; the supplemental
attractive-density scan uses its stated finite-occupation references.

## Hubbard-Stratonovich Fields

The Trotter decomposition uses two continuous auxiliary fields:

- `phi_1` couples to `n_b + n_c` with real coefficient `sqrt(-2 * Dtau * U1)`.
- `phi_2` couples to `n_b - n_c` with imaginary coefficient `i * sqrt(2 * Dtau * U2)`.

Equivalently, the continuous-HS coupling constants are:

```text
alpha_cont(U1) = sqrt(-2 * U1 * Dtau)
alpha_cont(U2) = i * sqrt(2 * U2 * Dtau)
```

At zero interaction strength the corresponding coupling vanishes. Thus the
main-text scan has only the relative-density phase field coupled to the bosons;
at `U=0` both couplings vanish.

After decoupling, the two flavor Hamiltonians are complex conjugates. The code samples the `b` flavor explicitly. The `c` flavor Green matrix is reconstructed with complex conjugation.

## Green Function Convention

For a fixed auxiliary-field configuration, the equal-time Green matrix is:

```text
G_ij(tau) = < b_i b_j^dagger >
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
| `pole_z` | Re/Im pairs for `z_a = 1 / mu_a(G)` |
| `pole_distance` | `min_a |z_a|` |
| `pole_x` | `-log10(pole_distance)` |
| `green_spectral_radius` | `max_a |mu_a(G)|` |
| `green_smax` | largest singular value of `G` |
| `log_weight` | `log_P_HS + 2 * sum log(s_a(G))` |

Here `mu_a(G)` are eigenvalues of the `b`-flavor Green matrix and `s_a(G)` are its singular values. For the continuous fields,

```text
log_P_HS = -0.5 * sum phi^2
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

The lattice uses the common triangular primitive vectors `a1=(1,0)` and
`a2=(1/2,sqrt(3)/2)`. The reciprocal basis uses the standard `2*pi`
convention, with `b_i dot a_j = 2*pi delta_ij`:
`b1=2*pi*(1,-1/sqrt(3))` and `b2=2*pi*(0,2/sqrt(3))`. For 3x3, the
code-facing K index above has integer momentum coordinate `(2,1)` and
therefore represents `K = (2/3)*b1 + (1/3)*b2 = (4*pi/3,0)`.
Older campaign data generated with the equivalent primitive-cell relabeling
`a2' = a2 - a1` describe the same triangular-lattice K channel, so those
benchmark results do not require rerunning solely for this coordinate-label
change.

The campaign uses `Lq = Lx * Ly` and
structure-factor normalization `1 / Lq^2`. For 3x3, `Lq = 9` and
`Lq^2 = 81`.

| DQMC output file | ED/schema name | Definition and normalization |
| --- | --- | --- |
| `density_total` | `density_total` | `( <N_b> + <N_c> ) / Lq` |
| `density_site_total` | postprocessed `IPR` input | site-resolved total density profile `rho_i = <n_{b,i} + n_{c,i}>` |
| `energy_density` | `energy_density` | physical energy density without the grand-canonical `-mu N` term: `kinetic + interaction_energy_density` |
| `doubleOcc` | `doubleOcc` | `(1/Lq) * sum_i <n_{b,i} n_{c,i}>` |
| `onsite_n2_up` | `onsite_n2_up` | `sum_i <n_{b,i}^2>`; DQMC estimator `sum_i (2*Grupc_ii^2 + Grupc_ii)` |
| `onsite_n2_do` | `onsite_n2_do` | `sum_i <n_{c,i}^2>`; DQMC estimator `sum_i (2*Grdoc_ii^2 + Grdoc_ii)` |
| `numsquare_up` | `numsquare_up` | `<N_b^2>`, where `N_b=sum_i n_{b,i}` |
| `numsquare_do` | `numsquare_do` | `<N_c^2>`, where `N_c=sum_i n_{c,i}` |
| `sf_K` | `S_SF_K` | `(1/Lq^2) * sum_{ij} exp(i K dot (r_i-r_j)) * <b_i^dagger b_j + c_i^dagger c_j>` |
| `psf_Gamma` | `S_PSF_Gamma` | `(1/Lq^2) * sum_{ij} <b_i^dagger c_i^dagger c_j b_j>` |
| `dw_K` | `S_DW_K` | `(1/Lq^2) * sum_{ij} exp(i K dot (r_i-r_j)) * <(n_{b,i}+n_{c,i})(n_{b,j}+n_{c,j})>` |

Equivalently, the plotted structure factors use the standard convention
`S_O(q) = <O_q^\dagger O_q> / Lq^2` with the phase
`exp(i q dot (r_i-r_j))`. With this convention they are real and
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

```text
IPR_rho = sum_i rho_i^2 / (sum_i rho_i)^2,
rho_i = <n_{b,i} + n_{c,i}>
```

For a non-negative density profile, `IPR_rho` lies between `1/Lq` and `1`.
On a uniform 3x3 profile it is `1/9`; if the density is localized on one
site it is `1`. If the blocked mean total density is too small, the IPR is
marked unreliable rather than plotted as a trusted comparison.

For the 3x3 hopping convention, the single-particle minimum is `epsilon_min=-3`.
Grand-canonical `U1=0` campaign points must therefore use `mu < -3`; otherwise
the free/equal-flavor growth channel makes the bosonic partition function
diverge. For stable `U1=U2=0` points, the campaign can use the exact
grand-canonical free-boson reference instead of a many-body cutoff ED run:

```text
n_k = 1 / (exp(beta * (epsilon_k - mu)) - 1)
rho = (2 / Lq) * sum_k n_k
S_SF(K) = (2 / Lq) * n_K
S_DW(K) = (2 / Lq^2) * sum_k n_{k+K} * (1 + n_k)
S_PSF(Gamma) = (1 / Lq^2) * sum_k n_k * n_{-k}
```

The factor of two in `rho`, `S_SF`, and `S_DW` is the two flavor layers.
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
