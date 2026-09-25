# Pairing Benchmark Observable Contract

This solver shares observable definitions with the number-conserving solver
at `src/number_conserving/`.
Definitions that exist in the no-pairing benchmark must keep the same
normalization, symbols, and \(K\)-point convention here.

## Shared Triangular Geometry

Both solvers use the triangular primitive vectors

\[
\mathbf a_1=(1,0),\qquad
\mathbf a_2=\left(\frac12,\frac{\sqrt3}{2}\right),
\]

and reciprocal vectors

\[
\mathbf b_1=2\pi\left(1,-\frac{1}{\sqrt3}\right),\qquad
\mathbf b_2=2\pi\left(0,\frac{2}{\sqrt3}\right).
\]

Thus \(\mathbf b_i\cdot\mathbf a_j=2\pi\delta_{ij}\). On a \(3\times3\)
lattice the benchmark \(K\) point is

\[
K=\frac{2}{3}\mathbf b_1+\frac{1}{3}\mathbf b_2
=\left(\frac{4\pi}{3},0\right),
\]

implemented as `Latt%inv_cell_list(2*Nlx/3+1, Nly/3+1)`.

## Observables Inherited From The No-Pairing Benchmark

| symbol | DQMC/ED key | definition |
| --- | --- | --- |
| \(\rho\) | `density_total` | \(\langle N_b+N_c\rangle/N_s\) |
| \(E\) | `energy_density` | paper-facing total physical energy \(E=N_s\,\texttt{energy\_density}=\langle H_{\rm phys}\rangle\), excluding \(-\mu N\); plots use \(-E\) |
| \(D_{bc}\) | `doubleOcc` | \(N_s^{-1}\sum_i\langle n_{b,i}n_{c,i}\rangle\) |
| \(M_b^{(2)}\) | `onsite_n2_up` | \(\sum_i\langle n_{b,i}^2\rangle\) |
| \(M_c^{(2)}\) | `onsite_n2_do` | \(\sum_i\langle n_{c,i}^2\rangle\) |
| \(N_b^2\) | `numsquare_up` | \(\langle(\sum_i n_{b,i})^2\rangle\) |
| \(N_c^2\) | `numsquare_do` | \(\langle(\sum_i n_{c,i})^2\rangle\) |
| \(S_{\rm SF}(K)\) | `S_SF_K` / `sf_K` | \(N_s^{-2}\sum_{ij,\alpha}e^{iK\cdot(r_i-r_j)}\langle a^\dagger_{\alpha i}a_{\alpha j}\rangle\) |
| \(S_{\rm DW}(K)\) | `S_DW_K` / `dw_K` | \(N_s^{-2}\langle n_K n_{-K}\rangle\), \(n_i=n_{b,i}+n_{c,i}\) |
| \(S_{\rm PSF}(\Gamma)\) | `S_PSF_Gamma` / `psf_Gamma` | \(N_s^{-2}\sum_{ij}\langle b_i^\dagger c_i^\dagger c_j b_j\rangle\) |
| \(\mathrm{IPR}_\rho\) | `IPR` | \(\sum_i\rho_i^2/(\sum_i\rho_i)^2\), \(\rho_i=\langle n_{b,i}+n_{c,i}\rangle\) |

The \(K\)-point structure factors are Hermitian positive structure-factor
channels with exact real non-negative thermal averages. DQMC files remain
two-column complex files because the Fourier sums are accumulated as complex
numbers. The real part is the comparison channel; the imaginary part is checked
with the analysis gate
\(|{\rm Im}|\le \max(10^{-14},3\,{\rm SEM}_{\rm Im})\), where the absolute
floor prevents roundoff-scale values from being amplified by a numerically tiny
SEM. A violation of this imaginary-part gate, or a real part that is negative
beyond the same roundoff/block-error scale, marks the observable unreliable.

## Pairing-Specific Extensions

The pairing benchmark adds observables that do not exist in the no-pairing
model. The operators in this table follow the solver convention
`b_code=b_paper`, `c_code=-c_paper`, in which the real pair coefficient is
positive. Consequently the manuscript's anomalous amplitude per site is
\(P_{\rm paper}=N_s^{-1}\sum_i\langle b_i c_i+b_i^+c_i^+\rangle_{\rm paper}
=-\texttt{pair_equal}\), and its energy contribution per site is
\(e_\Delta=-\Delta P_{\rm paper}=\Delta\,\texttt{pair_equal}\).
The four plotted benchmark observables are invariant under this phase change.
See [the Hamiltonian conventions](physics.md#manuscript-notation-and-pair-phase).

| symbol | key | definition |
| --- | --- | --- |
| \(P_\Delta\) | `pair_equal` | \(N_s^{-1}\sum_i\langle b_i c_i+b_i^\dagger c_i^\dagger\rangle\) |
| \(e_\Delta\) | `pairing_energy_density` | \(\Delta P_\Delta\), per site; total \(E_\Delta=N_s e_\Delta\) |
| \(e_U\) | `interaction_energy_density` | \((U_1+U_2)(M_b^{(2)}+M_c^{(2)})/N_s+2(U_2-U_1)D_{bc}\), per site |
| \(e_\mu\) | `chemical_energy_density` | \(-\mu\rho\), per-site diagnostic only |
| \(E_{\rm grand}\) | `grand_energy_density` | \(N_s\,\texttt{grand\_energy\_density}=E+N_s e_\mu\), diagnostic only |

The paper-facing energy panel remains \(-E\), not
\(-e\) and not \(-E_{\rm grand}\). The internal DQMC/ED key
`energy_density` is kept only as a storage and comparison convention.

## Extending observables

When changing a shared observable, compare the implementations in
`src/pairing/src/obser_equal.f90` and
`src/number_conserving/src/obser_equal.f90` as well as their ED operators.
For K-point quantities, check reciprocal coordinates, Hermiticity, and
normalization. Add the corresponding small ED identity or limiting-case test.
