# Implement a new model

BAFQMC can be developed beyond the two triangular-lattice models distributed
here. This guide maps a new Hamiltonian onto the actual solver, ED, and
measurement code. It covers changes to the lattice, hopping, pairing, flavors,
and interactions. For a parameter scan of an already implemented model, use
[the calculation recipes](agent-workflows.md). The
[observable reference](observables.md) defines the current measured operators
and output normalizations that extensions can build on.

Start from the Hamiltonian the researcher requested and carry the extension
through a working example and independent numerical checks. The worked
anisotropic-lattice example below is an **implementation blueprint**; the
current executable has no generic lattice or hopping-matrix input flag.

## Specify sites, operators, and the Hamiltonian

Distinguish unit cells, physical sites, flavor modes, and Nambu components.
For a crystal with $`n_{\mathrm{sub}}`$ sites per unit cell,

```math
N_{\mathrm{cell}}=L_xL_y,\qquad
N_s=n_{\mathrm{sub}}N_{\mathrm{cell}},\qquad
\mathbf r_{\mathbf R a}=R_1\mathbf a_1+R_2\mathbf a_2+\boldsymbol\delta_a.
```

Here $`a`$ is the sublattice index and $`\boldsymbol\delta_a`$ is its position
within the cell. With $`F`$ boson flavors there are $`M=F N_s`$ annihilation
operators. A full Nambu representation has $`2M`$ components. The existing
number-conserving solver propagates one flavor and reconstructs its conjugate
partner; its matrix dimension is consequently $`N_s`$ for the current model.

Write a site/flavor index convention and an explicit bond list, including
periodic image displacements. State whether each bond appears once with its
Hermitian conjugate or as two directed entries. Small periodic clusters can
have repeated image bonds and self-image bonds: summing their amplitudes and
deduplicating their endpoints are different Hamiltonians. The current code
sums the image contributions.

A useful general quadratic specification, using the paper's negative pairing
convention, is

```math
H_0
=\sum_{\alpha\beta} a_\alpha^+ h^{(0)}_{\alpha\beta}a_\beta
-\frac12\sum_{\alpha\beta}
\left(\Delta_{\alpha\beta}a_\alpha^+a_\beta^+
+\Delta_{\alpha\beta}^{*}a_\beta a_\alpha\right)+C,
\qquad h^{(0)}=\bigl(h^{(0)}\bigr)^\dagger,\quad \Delta=\Delta^T.
```

The physical Hamiltonian and its grand-canonical counterpart are

```math
H=H_0+H_{\mathrm{int}},\qquad
H_\mu=H-\mu N,\qquad N=\sum_\alpha a_\alpha^+a_\alpha,
\qquad h=h^{(0)}-\mu I.
```

The combined indices $`\alpha,\beta`$ include physical site and flavor, and
$`h`$ is the normal matrix used for the thermal trace. Hermitian hopping
requires the reverse amplitude to be the complex conjugate. Bosonic pairing is symmetric under
exchange of the two combined indices; the factor $`1/2`$ avoids double counting.
State the physical pairing phase explicitly. The present solver uses a
positive real interflavor pair coefficient, related to the paper's negative
coefficient by $`c_{\mathrm{code}}=-c_{\mathrm{paper}}`$.

Specify the interaction in operators, including its linear and constant
terms. For example, $`n^2=n(n-1)+n`$, so replacing one by the other also changes
the chemical-potential contribution. Keep physical energy and grand-canonical
energy separately defined throughout the implementation. For the paper's convention, use $`b,c`$ for the two flavors and write
$`H_{\mathrm{int}}=U_1(n_b-n_c)^2+U_2(n_b+n_c)^2`$, with
$`U_1\geq0`$ and $`U_2\leq0`$. The distributed examples use $`U_2=0`$
for the relative-density scan and $`U_1=1`$ while varying $`U_2<0`$ for the
total-density scan; both channels are part of the implemented model.

## Derive the decoupling and its sign protection

For an interaction written as a sum of Hermitian quadratic channels
$`H_{\mathrm{int}}=\sum_\ell g_\ell Q_\ell^2`$, derive the HS representation
used in each Trotter factor with time step $`\Delta\tau`$. A continuous
Gaussian channel has the identity

```math
e^{-\Delta\tau g Q^2}
=\int_{-\infty}^{\infty}\frac{d\phi}{\sqrt{2\pi}}
e^{-\phi^2/2}
e^{\sqrt{-2\Delta\tau g}\,\phi Q}.
```

For $`g<0`$ the field coefficient is real; for $`g>0`$ it is imaginary. If the
channels or kinetic factors do not commute, document their order and the
Trotter approximation. Carry every normal-ordering shift and scalar factor
into the weight, rather than placing all information in the matrix
propagator.

For each independent field configuration, establish the applicable
time-reversal symmetry (TRS), reflection positivity (RP), or conjugate-sector
construction **after decoupling**. Write the symmetry operator or reflection,
its action on sites/flavors/Nambu indices, and the conditions on each generated
quadratic factor. For a TRS route, record the antiunitary square and all
additional hypotheses of the sign-positivity criterion being used. For an RP
route, specify the reflected partition and the required form and sign of
cross-partition couplings. Symmetry of the original interaction alone does not
check these configuration-level conditions.

The simplest existing construction gives an explicit example. For equal real
hopping of the two flavors, $`U_1\geq0`$ and $`U_2\leq0`$ generate conjugate
single-flavor factors at the **same** field values. If $`B_c(\phi)=B_b(\phi)^*`$,
then, within the convergent trace domain,

```math
w(\phi)=p_{\mathrm{HS}}(\phi)
\left|\det\left[I-B_b(\phi)\right]\right|^{-2}\geq0.
```

The number-conserving code uses this relation directly. A model with different
flavor hoppings, flavor mixing, or a new HS channel needs its own matrix
representation and weight derivation whenever that relation changes. Taking
an absolute value of a new determinant is not a proof of its sign protection.
For complex hopping, this shortcut requires conjugate normal matrices,
$`h_c=h_b^*`$, together with conjugate HS factors. Equal complex hopping with
the same flux in both flavors generally does not satisfy that condition.
Replacing one flavor's flux by its opposite changes the requested Hamiltonian;
assess the actual same-flux model using its applicable symmetry construction.
If a proposed Hamiltonian lies outside an established class, identify which
condition is affected and develop the requested model accordingly; do not
silently change its physical couplings to recover the old class.

If a requested model instead requires complex-weight sampling, implement that
estimator explicitly. The paired solver's local ratio-phase diagnostics are
not phase reweighting: its current measurements are ordinary averages. A
complex-weight calculation needs the full configuration phase
$`w(\phi)=|w(\phi)|e^{i\theta(\phi)}`$ and the corresponding estimator

```math
\langle O\rangle
=\frac{\langle O(\phi)e^{i\theta(\phi)}\rangle_{|w|}}
{\langle e^{i\theta(\phi)}\rangle_{|w|}}.
```

This entails phase-aware measurement and analysis, including the scalar and
branch contributions to $`\theta`$, or a separately established positive
construction for the same physical Hamiltonian.

## Establish the physical trace and matrix representation

Sign protection and trace convergence are separate parts of the construction.
Specify the parameter region in which the grand-canonical physical model has
a finite thermal trace, and establish the domain required by the decoupled
quadratic traces. A finite ED occupation cutoff defines a finite-dimensional
reference; convergence toward an untruncated model is a further calculation.

For a quadratic model, positivity of the Hermitian bosonic energy matrix is a
useful sufficient condition for a finite thermal trace. One sufficient bound
for the quadratic form above is

```math
\lambda_{\min}(h)>\|\Delta\|_2.
```

Real Bogoliubov frequencies alone do not establish positivity of that energy
matrix. For the current $`U_2=0`$ construction with uniform onsite pairing and
real hopping, the corresponding sufficient bound is

```math
\mu<\varepsilon_{\min}-|\Delta|,
```

where $`\varepsilon_{\min}`$ is the hopping-band minimum. The relative-density
HS factors are unitary. The triangular benchmark has
$`\varepsilon_{\min}=-3t`$; a different graph needs its own band bound. New
nonunitary HS channels require the corresponding trace-domain analysis.

Keep the quadratic energy matrix distinct from the bosonic commutator matrix
that is exponentiated. In the full Nambu formulation the latter is generally
non-Hermitian even for a Hermitian physical Hamiltonian. The paired solver's
`exp_general_matrix` uses general-complex diagonalization. A new pairing
matrix must retain the correct particle/hole signs and conjugations; using a
Hermitian eigensolver on that commutator matrix changes the propagator.

The paired representation also has a scalar normal-ordering contribution.
For the current total-density HS field $`x_2=\sqrt{-2U_2\Delta\tau}\,\phi_2`$,
the scalar factor is $`e^{-x_2}`$, and a proposal contributes
$`e^{-(x_2'-x_2)}`$ to `ratio_constant`. Re-derive this factor for a different
quadratic generator or Nambu convention.

The local determinant factor currently uses
`det_Pblock**(-0.5d0)`. The code records the phase and samples the magnitude;
it does not maintain an explicit branch-continuation state. Establish the
square-root branch associated with the physical trace for a new paired
model, anchored in a known convergent limit, and test the full ratio including
its scalar. If the model needs branch tracking, implement and verify that
state as part of the extension. Changing the determinant exponent alone does
not convert between a full and a reduced Nambu representation.

## Change the actual implementation

The following map points to active source files. Both components use similar
names, with different matrix dimensions and update formulas.

| Responsibility | Files and symbols | Required coordination |
| --- | --- | --- |
| Input and dimensions | `src/<solver>/src/calc_basic.f90`: `read_input`, `Params_set` | Read and broadcast new parameters; separate cells, sites, flavors, and sectors; record them in run metadata. `RT=1` is currently assigned here. |
| Crystal and graph | `src/<solver>/src/lattice.f90`: `Lattice_make` | Update `L_bonds`, `LT_bonds`, site/sector index maps, real/reciprocal vectors, and Fourier phases. |
| Quadratic propagation | `src/<solver>/src/non_interact.f90`: `def_hamT`, `opT_set`; paired `exp_general_matrix` | Build the new hopping/pairing matrix and its inverse-time factors, keeping the physical and commutator conventions distinct. |
| Fields and channel schedule | `fields.f90`: `AuxConf_make`; `model.f90`: `Model_init`; `local_sweep.f90` | Update field support, channel allocation, initial/restart I/O, and both sweep directions. The active schedule explicitly uses two channels. |
| HS factors | `operator_Hubbard.f90`: `opU_set`, `opU_get_delta`, `opU_mmult_L/R` | Implement the derived operator, field measure, scalar, and local matrix change. In the paired code, coupling sign currently selects the channel type. |
| Ratios and Green updates | `localU.f90`: `LocalU_metro`; `multiply.f90`; `stabilization.f90` | Match update support and rank, derive the determinant ratio, and verify stabilized propagation against direct matrix products. |
| Measurements | `obser_equal.f90`: `Obs_equal_calc`; `fourier_trans.f90` | Change kinetic, interaction, pairing, and momentum estimators with the Hamiltonian; carry definitions into file output and analysis. |
| Campaigns and references | `src/<solver>/run_paper.py`, `src/<solver>/benchmarks/campaign*.py`, ED files below | Generate consistent solver/ED inputs, validate the implemented model, and report every new parameter and reference cutoff. |

Consult the source directly through the
[number-conserving component](../src/number_conserving/src/) and
[paired component](../src/pairing/src/). In particular:

- **A sublattice count is not a complete lattice definition.** The current
  `Norb=1` bond construction explicitly targets destination orbital 1.
  A kagome extension needs three physical sublattices, intracell positions,
  intercell bonds, and consistent orbital indices. Audit allocations and
  normalizations using `Lq`, `Ndim`, `Nsite`, `Norb`, `Nsub`, and `Nbond`.
  Sublattice count, coordination number, and forward-bond storage count are
  distinct. `Nbond=3` and the spatial/time columns of `LT_bonds` are currently
  fixed to the triangular graph. Check tensor extents in Fourier helpers
  such as `m_write_k_3` and `m_write_reciprocal_3` against their callers;
  orbital-correlation dimensions must follow the new orbital indexing. The paired
  `Nsec=4` labels $`(b,c,b^+,c^+)`$; it is not the sublattice count.
- **Quadratic and interaction changes affect different update support.**
  New deterministic hopping or bond pairing may retain the existing onsite
  density-HS support. A bond HS field, flavor-off-diagonal channel, or new
  interaction changes it. The number-conserving `LocalU_metro` applies a
  rank-one diagonal-site update. The paired version selects four sectors of
  one site and uses only diagonal entries of its $`4\times4`$ `Delta`.
  Adding off-diagonal entries to that array alone does not implement the
  required Woodbury update.
- **Channel signs encode operators.** The paired `opU_set` uses the positive $`U_1`$ channel for relative density
  and the nonpositive $`U_2`$ channel for total density. Opposite signs or
  additional operators require explicit channel definitions; renaming `Op_U1`
  or changing a JSON number does not change that dispatch.
- **Measurements contain the old Hamiltonian explicitly.** Kinetic energy
  sums the scalar `RT` over `L_bonds`; paired energy uses onsite
  `RDelta*pair_equal`. Update these when propagation changes. The current
  paired implementation extracts the full four-sector contractions; it does
  not reconstruct its interacting paired $`c`$ sector by conjugating $`b`$.
- **Momentum and normalization follow the new crystal.** Use
  $`e^{i\mathbf q\cdot(\mathbf r_{\mathbf R a}-\mathbf r_{\mathbf R'b})}`$,
  including sublattice positions, and state whether the result is a matrix
  in sublattice indices or a summed physical structure factor. Audit every
  $`N_s^{-1}`$ and $`N_s^{-2}`$ factor. The triangular K index is not a generic
  ordering momentum; current noncommensurate Fortran runs leave its output
  zero while paired ED omits it. That zero represents an unavailable K
  estimator. Define and implement the new model's actual momenta.

The existing paired time-dependent observable routine `Obs_tau_calc` is outside
the supported equal-time workflow. A research task involving imaginary-time
correlations must add an estimator together with its propagation checks.

## Extend ED and analysis with the solver

The [number-conserving ED driver](../src/number_conserving/benchmarks/ed/EDtriangle_quspin_3x3.py)
contains `_hamiltonian_static`, `_hopping_two_species`,
`_translation_permutations`, the two-species basis construction, and K-point
operators. Its public validator currently restricts interacting reference
calculations to $`3\times3`$. Rework those restrictions and symmetry blocks
when introducing other sizes or graphs; do not retain translation blocks
that the new Hamiltonian no longer conserves. The free reference in
`src/number_conserving/benchmarks/campaign_analysis.py` also contains the
triangular dispersion and needs the new spectrum.

For pairing, update
[geometry_triangle.py](../src/pairing/benchmarks/ed/geometry_triangle.py),
the general driver's `build_basis`, `build_hamiltonian`, and observable
operators, and the [tensor ED implementation](../src/pairing/benchmarks/ed/ed_pairing_triangle_tensor.py).
The present basis has two flavors per physical site. Complex hopping/pairing
also requires complex coefficients and a suitable Hamiltonian dtype. Preserve
the distinction between physical site count and QuSpin's `basis.Ns`, which is
the many-body Hilbert-space dimension.

Add a model identifier and explicit parameter fields in the new input schema
when several implementations coexist. Extend the readers and validators that
consume those fields. A new manifest flag has an effect only after those
readers and numerical kernels implement it. Keep generated tables
self-describing: model, geometry, parameters, units, normalization, sampling,
SEM, ED cutoffs, and energy convention. Use separate examples and campaign
outputs so the original triangular cases and new model are both runnable.
Update reference data intentionally when definitions or calculations change,
and explain the change alongside the result.

## Worked blueprint: anisotropic triangular hopping

Consider a requested extension with the same two flavors and onsite channels,
but three independent real hopping amplitudes:

```math
H_{t}=\sum_{\mathbf R,\sigma}\sum_{\nu=1}^{3}t_\nu
\left(a_{\mathbf R\sigma}^+a_{\mathbf R+\mathbf d_\nu,\sigma}
+a_{\mathbf R+\mathbf d_\nu,\sigma}^+a_{\mathbf R\sigma}\right),
\qquad
\mathbf d_1=(1,0),\quad\mathbf d_2=(0,1),\quad\mathbf d_3=(-1,1).
```

The displacement pairs are integer cell coordinates. Keep the physical
triangular primitive vectors from the algorithm guide. The dispersion becomes

```math
\varepsilon(\mathbf k)=2\left[
t_1\cos(\mathbf k\cdot\mathbf a_1)
+t_2\cos(\mathbf k\cdot\mathbf a_2)
+t_3\cos\bigl(\mathbf k\cdot(\mathbf a_2-\mathbf a_1)\bigr)
\right].
```

This change preserves the conjugate-flavor relation for equal real hoppings
and the original onsite HS channels. For $`U_2=0`$ and uniform real onsite
pairing, a simple sufficient thermal bound is

```math
\mu<-2\left(|t_1|+|t_2|+|t_3|\right)-|\Delta|.
```

It follows from a lower bound on the hopping spectrum and need not be the
tight band minimum. As a concrete implementation test, use
$`(t_1,t_2,t_3)=(1,0.8,0.35)`$, $`\mu=-5`$, $`\beta=1`$ and first set
$`U_1=U_2=0`$. Test $`\Delta=0`$ and, for the paired extension, $`\Delta=0.1`$.
Then exercise the existing relative-density interaction at $`U_1=0.7`$.
These are proposed extension tests, not bundled runnable cases.

An agent implementing this request should:

1. Add and broadcast the three amplitudes in the input representation;
   retain the isotropic values $`(1,1,1)`$ as the original-model default.
   Connect each `L_bonds(:,nu)` entry to its amplitude in `def_hamT`.
2. Make the kinetic estimator use the identical amplitudes and Hermitian
   conjugates. Keep onsite HS algebra and its local support unchanged after
   checking the conjugate-sector relation for the new quadratic factors.
3. Update both ED hopping constructions and the analytic free dispersion.
   Choose a small finite-occupation reference or extend the current NC ED
   size interface as needed; record its basis definition.
4. Compute free thermal density, energy, and momentum occupation independently
   from the dispersion. For pairing, use
   $`\omega_{\mathbf k}=\sqrt{(\varepsilon_{\mathbf k}-\mu)^2-\Delta^2}`$
   in the stable Gaussian case. Check kinetic energy against derivatives
   with respect to each $`t_\nu`$.
5. Verify a nonuniform fixed HS field by direct dense multiplication, then
   run an interacting BAFQMC/ED comparison. Recover the original triangular
   answers at $`(1,1,1)`$ with the existing paper cases.
6. Add a small curated example, exact commands, resource measurements, and
   a model definition to the documentation. Keep larger generated scans in
   their own output directory.

Setting $`t_3=0`$ removes the diagonal graph bonds. To describe a physical
square lattice as well, set orthogonal primitive vectors and update the
reciprocal coordinates and momentum observables. This illustrates why graph
connectivity and its real-space embedding are separate pieces of a model.

## Independent verification and deliverables

Build a small dense reference directly from the new model definition. It
should assemble site operators, hopping, pairing, and interaction terms
independently of the production helpers. Compare Hamiltonian spectra and
thermal observables in a common finite Fock basis; also check Hermiticity,
pairing symmetry, energy decomposition, and relevant thermodynamic or gauge
identities. The approach in
[test_ed_physics.py](../tests/physics/test_ed_physics.py) is a starting point.

For explicit nonuniform fields, form the ordered short-time matrices and
their full products independently. Compare Green functions, the complete
weight including scalars, and proposed-update ratios. Test several field
configurations and stabilization intervals. Extend
[test_fixed_field_solver.py](../tests/physics/test_fixed_field_solver.py) to
the new geometry/channel; its current coverage checks interacting HS
propagation but does not test a nonzero Metropolis proposal. A change to the
update formula additionally needs old/new dense ratios and updated Green
matrices for deterministic nonzero proposals.

Use analytic free or Bogoliubov limits as in
[test_gaussian_solvers.py](../tests/physics/test_gaussian_solvers.py), then
perform an interacting run that actually exercises the new terms. Report
means, SEM, and reference differences without a universal statistical sigma
gate. Test the new model explicitly: passing the old triangular benchmarks
only establishes that those cases remain consistent.

Run the existing checks alongside the added model tests:

```bash
make check
make physics
```

Use `PYTHON_ED=/path/to/quspin/python` when ED has a separate environment.
These checks have independent physical references; a new model should add
equally direct coverage. Full 22-point production remains available through
`python3 reproduce.py` when required by the requested validation scope.

Deliver a documented Hamiltonian and lattice, the HS/symmetry and trace
derivation appropriate to its parameters, working solver and reference
implementations, a runnable small example, independent tests, and measured
computational cost. Record intentional changes to existing inputs or results
as part of normal development. The researcher should be able to run both
the original examples and the new physical model from this checkout.
