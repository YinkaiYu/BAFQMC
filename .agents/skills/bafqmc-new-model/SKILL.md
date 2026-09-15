---
name: bafqmc-new-model
description: Implement a different BAFQMC lattice, hopping, pairing, flavor structure, or interaction, with a derived HS construction and coordinated solver, ED, observable, input, and physical-test changes. Use for model extensions requiring source changes; use bafqmc-new-calculation for scans of an already implemented model.
---

# Implement a new physical model

Read [AGENTS.md](../../../AGENTS.md) and
[the model-development guide](../../../docs/model-development.md). Use its
source map and worked anisotropic-lattice blueprint for the parts relevant to
the requested Hamiltonian. This is an implementation workflow, not a
parameter-scan recipe.

Establish the physical operators, crystal/boundaries, site and flavor counts,
hopping/pairing amplitudes, interactions, observables, and desired validation
from the user's request. Continue the authorized research work through code
and a runnable example. Resolve genuinely missing physical choices; make
routine implementation decisions independently. Do not substitute another
Hamiltonian or create an approval stage merely because the work changes a model.

Before encoding a new term, derive its quadratic/HS factors, normal-ordering
shifts and scalar weight. Establish the relevant sign-protection conditions
for every independent auxiliary-field configuration. State the physical
thermal-trace regime separately. Use the existing conjugate-sector shortcut
only where the new factors actually satisfy that relation; a paired extension
must also establish its Nambu representation and determinant-root branch.
In particular, equal complex hopping with the same flux in both flavors
generally does not give the required conjugate factors. Opposite flavor fluxes
would define a different Hamiltonian; implement and assess the requested one.
If that implementation uses complex weights, implement full configuration-phase
reweighting in measurements and analysis. Existing local ratio-phase diagnostics
do not perform it.

Follow the actual code contracts:

- The current crystal has one physical site per cell. For multiple
  sublattices define positions and intercell bonds, and separate
  $N_s=n_{\mathrm{sub}}L_xL_y$ from flavor and Nambu dimensions. Increasing
  `Norb` alone leaves the current destination-orbital and estimator assumptions.
  Audit `Nbond`, space-time bond slots, and Fourier tensor dimensions separately.
- Hopping is currently uniform real `RT=1`; paired inputs supply real onsite
  `RDelta`. Implement readers, MPI broadcasts, quadratic matrices, ED and
  estimators when introducing other amplitudes or a new graph. There is no
  generic lattice flag to set in a manifest.
- Deterministic quadratic changes can retain onsite HS support when its
  derivation remains valid. New HS operators may require new support/rank,
  field storage, channel scheduling, scalars and update formulas. The NC
  update is diagonal rank one; the paired update uses a diagonal four-sector
  block at one site. Coupling sign currently selects the paired channel type.
- Keep Hermitian hopping, symmetric bosonic pairing, full Nambu block
  conventions, physical energy, and new momentum/form-factor normalizations
  consistent. The NC interacting ED entry point currently restricts its
  reference to $3\times3$; adapt that interface and its symmetry blocks when
  the new model requires it.

Implement all connected solver, ED, observable, input and analysis changes.
Use separate named examples/campaigns so both models remain runnable; update
existing data intentionally when the task calls for it and document the reason.
Large raw measurements stay in generated output directories.

Add independent small-system checks built from the requested Hamiltonian:
finite-Fock spectra/thermal traces; analytic free or Bogoliubov limits; dense
fixed-field products and full weight ratios; and deterministic nonzero
proposal checks if an update rule changes. Exercise the new interaction in a
live calculation. Run `make check` and `make physics` as applicable, retaining
the original-model regression coverage. Report physical results with SEM and
reference differences without an arbitrary sigma acceptance rule.

Deliver the implementation, its physical definition and derivation, runnable
inputs/commands, independent verification and a measured resource estimate.
Once the model is implemented, route ordinary scans to
[bafqmc-new-calculation](../bafqmc-new-calculation/SKILL.md).
