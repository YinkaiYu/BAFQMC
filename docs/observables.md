# Physical observables and output files

Both solvers measure the quantities below on the periodic, single-orbital
triangular lattice. The definitions use the paper's operators and
$`N_s=L_xL_y`$ (`Lq` in the code):

```math
\hat n_{b,i}=\hat b_i^+\hat b_i,\qquad
\hat n_{c,i}=\hat c_i^+\hat c_i,\qquad
\hat N_b=\sum_i\hat n_{b,i},\qquad
\hat N_c=\sum_i\hat n_{c,i}.
```

Expectations are grand-canonical thermal averages. `up` and `do` in filenames
mean $`b`$ and $`c`$. All files listed here are written by both solvers unless
marked **paired only**. The paired anomalous amplitude has a phase convention
explained below; density, hopping, and the three structure factors agree
directly with the paper's convention.

## Densities and particle-number moments

Define the mean density and the onsite second-moment sums by

```math
\rho_i=\langle\hat n_{b,i}+\hat n_{c,i}\rangle,\qquad
\rho=\frac{\langle\hat N_b+\hat N_c\rangle}{N_s},\qquad
M_b=\sum_i\langle\hat n_{b,i}^2\rangle,\quad
M_c=\sum_i\langle\hat n_{c,i}^2\rangle.
```

| Output file | Definition |
| --- | --- |
| `density_up`, `density_do` | $`\langle\hat N_b\rangle/N_s`$, $`\langle\hat N_c\rangle/N_s`$ |
| `density_total` | $`\rho`$ |
| `density` (**paired only**) | Alias for $`\rho`$ |
| `density_site_total` | The $`N_s`$ entries $`\rho_i`$ |
| `num_up`, `num_do` | $`\langle\hat N_b\rangle`$, $`\langle\hat N_c\rangle`$ |
| `numsquare_up`, `numsquare_do` | $`\langle\hat N_b^2\rangle`$, $`\langle\hat N_c^2\rangle`$ |
| `onsite_n2_up`, `onsite_n2_do` | $`M_b`$, $`M_c`$ |
| `local_numsquare` (**paired only**) | $`(M_b+M_c)/N_s`$ |
| `doubleOcc` | $`D=N_s^{-1}\sum_i\langle\hat n_{b,i}\hat n_{c,i}\rangle`$ |
| `squareOcc` | $`Q=(2N_s)^{-1}\sum_i\langle\hat n_{b,i}(\hat n_{b,i}-1)+\hat n_{c,i}(\hat n_{c,i}-1)\rangle`$ for the current models |

`numsquare_*` includes correlations between different sites; `onsite_n2_*`
contains only onsite terms. These files store totals, whereas `doubleOcc`,
`squareOcc`, and `local_numsquare` are per site.

The exact auxiliary-field estimator retained under the name `squareOcc` is

```math
Q_\phi=\frac1{N_s}\sum_i\mathrm{Re}\left[
\langle\hat n_{b,i}\rangle_\phi^2+
\langle\hat n_{c,i}\rangle_\phi^2\right].
```

For the current Hamiltonians, conservation of $`\hat N_b-\hat N_c`$ makes
same-flavor anomalous contractions vanish, so Wick's theorem gives the
normal-ordered definition in the table. A model with same-flavor pairing
requires the additional anomalous contractions. The paired `local_numsquare`
and `onsite_n2_*` estimators already include those contractions explicitly.

The density-profile IPR is a postprocessed quantity, with no separate solver
output file:

```math
\mathrm{IPR}_\rho=
\frac{\sum_i\rho_i^2}{\left(\sum_i\rho_i\right)^2}.
```

It uses the averaged density profile and equals $`1/N_s`$ for a uniform
nonzero profile. It is distinct from a particle-number second moment.

## Energy and the anomalous pair amplitude

The physical Hamiltonian is $`\hat H=\hat H_t+\hat H_U+\hat H_\Delta`$, where

```math
\hat H_t=t\sum_{\langle ij\rangle}\left(
\hat b_i^+\hat b_j+\hat b_j^+\hat b_i+
\hat c_i^+\hat c_j+\hat c_j^+\hat c_i\right),
```

```math
\hat H_U=U_1\sum_i(\hat n_{b,i}-\hat n_{c,i})^2+
U_2\sum_i(\hat n_{b,i}+\hat n_{c,i})^2,
```

```math
\hat H_\Delta=-\Delta\sum_i
\left(\hat b_i^+\hat c_i^++\hat b_i\hat c_i\right).
```

The paired executable accepts real $`\Delta`$; the number-conserving solver
has $`\Delta=0`$. The implemented model uses the two coefficients $`U_1`$ and
$`U_2`$ directly, with values specified by each benchmark input file.
Each nearest-neighbor bond in $`\hat H_t`$ is counted once, with both hopping
directions shown explicitly. Define $`e_t=\langle\hat H_t\rangle/N_s`$,
$`e_U=\langle\hat H_U\rangle/N_s`$, and $`e_\Delta=\langle\hat H_\Delta\rangle/N_s`$.

The paired implementation uses

```math
\hat b_{i,\mathrm{code}}=\hat b_i,\qquad
\hat c_{i,\mathrm{code}}=-\hat c_i,
```

so its quadratic matrix contains a positive $`\Delta`$ pair term. Consequently,

```math
\begin{aligned}
P_{\mathrm{paper}}&=\frac1{N_s}\sum_i
\langle\hat b_i\hat c_i+\hat b_i^+\hat c_i^+\rangle,\\
P_{\mathrm{code}}&=\frac1{N_s}\sum_i
\langle\hat b_{i,\mathrm{code}}\hat c_{i,\mathrm{code}}+
\hat b_{i,\mathrm{code}}^+\hat c_{i,\mathrm{code}}^+\rangle
=-P_{\mathrm{paper}},\\
e_\Delta&=\Delta P_{\mathrm{code}}=-\Delta P_{\mathrm{paper}}.
\end{aligned}
```

| Output file | Definition |
| --- | --- |
| `kinetic` | $`e_t`$, including $`t`$ and both flavors |
| `interaction_energy_density` | $`e_U=(U_1+U_2)(M_b+M_c)/N_s+2(U_2-U_1)D`$ |
| `pair_equal` (**paired only**) | $`P_{\mathrm{code}}`$; negate it for the paper's pair amplitude |
| `pairing_energy_density` (**paired only**) | $`e_\Delta=\Delta P_{\mathrm{code}}`$ |
| `energy_density` | $`e=\langle\hat H\rangle/N_s=e_t+e_U+e_\Delta`$ |
| `chemical_energy_density` (**paired only**) | $`e_\mu=-\mu\rho`$ |
| `grand_energy_density` (**paired only**) | $`\langle\hat H-\mu(\hat N_b+\hat N_c)\rangle/N_s=e-\mu\rho`$ |

`energy_density` excludes the chemical-potential term. The total physical
energy is $`E=N_s e`$; the paper plots $`-E`$. Multiplying by $`N_s`$ also multiplies
the standard error by $`N_s`$. The number-conserving solver does not write
separate chemical or grand-energy files; they follow from the same identities.

## Structure factors at K and Gamma

The triangular-lattice primitive vectors and the measured momenta are

```math
\mathbf a_1=(1,0),\qquad
\mathbf a_2=\left(\frac12,\frac{\sqrt3}{2}\right),\qquad
\mathbf K=\left(\frac{4\pi}{3},0\right),\qquad
\boldsymbol\Gamma=(0,0).
```

All three structure factors use $`N_s^{-2}`$ normalization:

```math
S_{\mathrm{SF}}(\mathbf K)=\frac1{N_s^2}\sum_{ij}
e^{i\mathbf K\cdot(\mathbf r_i-\mathbf r_j)}
\langle\hat b_i^+\hat b_j+\hat c_i^+\hat c_j\rangle,
```

```math
S_{\mathrm{DW}}(\mathbf K)=\frac1{N_s^2}\sum_{ij}
e^{i\mathbf K\cdot(\mathbf r_i-\mathbf r_j)}
\langle(\hat n_{b,i}+\hat n_{c,i})(\hat n_{b,j}+\hat n_{c,j})\rangle,
```

```math
S_{\mathrm{PSF}}(\boldsymbol\Gamma)=\frac1{N_s^2}\sum_{ij}
\langle\hat b_i^+\hat c_i^+\hat c_j\hat b_j\rangle.
```

| Output file | Observable | Availability |
| --- | --- | --- |
| `sf_K` | $`S_{\mathrm{SF}}(\mathbf K)`$ | Both $`L_x`$ and $`L_y`$ divisible by 3 |
| `dw_K` | $`S_{\mathrm{DW}}(\mathbf K)`$ | Both $`L_x`$ and $`L_y`$ divisible by 3 |
| `psf_Gamma` | $`S_{\mathrm{PSF}}(\boldsymbol\Gamma)`$ | Every supported lattice size |

`dw_K` measures the total density. These are full correlation functions,
including disconnected contributions, without subtraction of products of
mean values. The pair structure factor is unchanged by the code-to-paper
phase rotation.

When the lattice cannot represent K, the `sf_K` and `dw_K` channels are
unavailable on that lattice; the executable writes zero entries to preserve the
output schema. This does not affect `psf_Gamma`.

## Density-correlation file families

The writer Fourier-transforms translation-averaged density correlations with
the positive phase convention

```math
C_{bb}(\mathbf q)=\frac1{N_s^2}\sum_{ij}
e^{i\mathbf q\cdot(\mathbf r_i-\mathbf r_j)}
\langle\hat n_{b,i}\hat n_{b,j}\rangle,
```

and analogously defines $`C_{cc}`$ and $`C_{bc}`$ by replacing the two density
operators. The translation average contributes $`1/N_s`$ and the Fourier
transform contributes another $`1/N_s`$. The current writer selects only
$`\mathbf q=\boldsymbol\Gamma`$:

| Output file | Quantity written |
| --- | --- |
| `den_upup_sub11` | $`C_{bb}(\boldsymbol\Gamma)=\langle\hat N_b^2\rangle/N_s^2`$ |
| `den_dodo_sub11` | $`C_{cc}(\boldsymbol\Gamma)=\langle\hat N_c^2\rangle/N_s^2`$ |
| `den_updo` | $`C_{bc}(\boldsymbol\Gamma)=\langle\hat N_b\hat N_c\rangle/N_s^2`$ |

The `sub11` suffix denotes the single orbital's two correlation indices.
These files contain one complex value per bin, without a momentum coordinate
or a full momentum grid. They do not subtract disconnected contributions.
In particular, `den_updo` includes different-site correlations, whereas
`doubleOcc` contains only onsite correlations.

## Rows, columns, and uncertainty

Files are appended in the run's working directory. Start each independent
run in a fresh directory. Each equal-time row is a bin mean: the local sweep
averages measurements over imaginary-time slices and sweeps, and rank zero
then writes the average over MPI ranks. For `Nsweep` bidirectional sweeps,
each rank contributes `2 * Nsweep * Ltrot` equal-time measurements per bin.

| File type | Columns in each row |
| --- | --- |
| Scalar densities, moments, energies, and `pair_equal` | One real value |
| `density_site_total` | $`N_s`$ real values, one for each site |
| `sf_K`, `dw_K`, `psf_Gamma`, and `den_*` | Two values: real part, imaginary part |

For `density_site_total`, zero-based site index $`i=x+L_x y`$ orders the
cells, with $`x=0,\ldots,L_x-1`$ running fastest. No column labels are stored
in these raw files. The exact structure factors and density correlations
above are real; the imaginary columns retain the finite-sample estimator
and provide a consistency diagnostic.

Postprocessing blocks the retained bin series to estimate the uncertainty
of its mean. For $`B`$ block means $`\bar x_j`$,

```math
\bar x=\frac1B\sum_{j=1}^B\bar x_j,\qquad
\mathrm{SEM}(\bar x)=
\sqrt{\frac{\sum_{j=1}^B(\bar x_j-\bar x)^2}{B(B-1)}}.
```

The reported `stderr` is this standard error, rather than the standard
deviation of raw bins. Complex columns are analyzed separately. Raw samples
are generated locally; the small processed tables retain means and their
standard errors. See [reproduction](../benchmarks/paper/README.md) for the paper workflow.

## Active measurements and diagnostics

Equal-time measurements begin with the first output bin. `Nthermal` controls
when the time-displaced path is entered; it does not discard equal-time rows.
Auxiliary-field warmup and any analysis-time bin discard are separate steps.

Both current `obser_tau.f90` estimators and `m_write_obs_tau` writers are empty.
Enabling `is_tau` therefore produces no time-displaced physical observables,
susceptibilities, or Matsubara-frequency data. The paired solver accumulates
an internal single-particle correlation, but its `green` write call is
commented out. The generic Fourier-output helper methods do not by themselves
enable additional output files.

The number-conserving files `pole_z`, `pole_distance`, `pole_x`,
`green_spectral_radius`, `green_smax`, and `log_weight` describe configuration
weights and numerical conditioning. Unlike equal-time observable rows, these
contain one configuration sample per bin **per MPI rank**, gathered in rank
order. Their definitions are in the [number-conserving physics guide](solvers/number_conserving/physics.md#continuous-pole-diagnostics).
Acceptance, phase, and matrix-stability statistics in `info.txt` are also
diagnostics. Logs, seeds, and field configurations record runtime state.

## Implementation map

| Responsibility | Number-conserving solver | Paired solver |
| --- | --- | --- |
| Wick estimators and normalizations | [obser_equal.f90](../src/number_conserving/src/obser_equal.f90) | [obser_equal.f90](../src/pairing/src/obser_equal.f90) |
| Active files, column layouts, MPI averaging | [fourier_trans.f90](../src/number_conserving/src/fourier_trans.f90) | [fourier_trans.f90](../src/pairing/src/fourier_trans.f90) |
| Lattice ordering and Fourier normalization | [lattice.f90](../src/number_conserving/src/lattice.f90) | [lattice.f90](../src/pairing/src/lattice.f90) |
| Bin loop and measurement flags | [main.f90](../src/number_conserving/src/main.f90) | [main.f90](../src/pairing/src/main.f90) |
| Time-displaced estimator status | [obser_tau.f90](../src/number_conserving/src/obser_tau.f90) | [obser_tau.f90](../src/pairing/src/obser_tau.f90) |

For extensions, update each estimator, its normalization, and the active
writer together. The [model-development guide](model-development.md) maps
the corresponding Hamiltonian, basis, and update changes.
