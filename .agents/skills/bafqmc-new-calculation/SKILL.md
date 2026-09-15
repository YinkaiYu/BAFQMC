---
name: bafqmc-new-calculation
description: Prepare and run a BAFQMC point or parameter scan for an already implemented model, with consistent solver and ED inputs and a stated computation budget. Route different lattices, hopping, pairing or interaction operators requiring source changes to bafqmc-new-model first.
---

# Start a research calculation

Read [AGENTS.md](../../../AGENTS.md), [the algorithm guide](../../../docs/algorithm.md),
and the applicable section of [the research recipes](../../../docs/agent-workflows.md).
Commands run from the repository root in Linux or WSL.

Establish the Hamiltonian, lattice size, temperature, observables, requested
precision, and available resources from the user's task. Make routine setup
choices independently; resolve missing physical assumptions before computing
a different model. If the request changes the lattice, hopping/pairing
matrices, flavor structure, or interaction operators, use
[bafqmc-new-model](../bafqmc-new-model/SKILL.md) and
[the extension guide](../../../docs/model-development.md) to implement it,
then return here for the requested scan. Continue the authorized model work;
do not reduce an implementation request to a scan of the existing Hamiltonian.
The current solvers use periodic triangular hopping with `RT=1` assigned in
each `src/<solver>/src/calc_basic.f90`, so a JSON `t` edit alone does not
change the executable's hopping.

For a model already added to this checkout, start from its documented example,
input schema, solver and ED entry points. Keep its geometry and operator
definitions when making the requested scan. The templates below describe the
two distributed triangular-lattice models; adapt them to the implemented
model's interface instead of replacing that model with the template.

For the distributed models, use `src/number_conserving/` for zero pairing, or
`src/pairing/` for onsite pairing. The main-text coupling maps to `U1=0, U2=U`.
`U1` and `U2` always denote the total-density and relative-density channels.
The pairing solver's operator phase is documented in its physics guide;
account for it when reporting an anomalous pair amplitude.

Create a separate campaign under an ignored `runs/<project>/` directory:

- **Number conserving:** follow the input-copy recipe in
  [agent-workflows.md](../../../docs/agent-workflows.md). Keep
  `paramC_sets.txt`, ED `params.json`, and manifest `parameters` consistent.
  Relative input paths resolve from the manifest's directory.
- **Pairing:** copy
  `src/pairing/benchmarks/campaigns/triangle_pairing_pipeline_smoke.json`
  into the new workspace, then edit its physical parameters, `deltas`,
  sampling/warmup, analysis block size, and ED cutoffs. This manifest generates
  both sets of inputs. Preserve case order when selecting a subset because
  case position determines its seed.

Both supplied starting examples have only eight measurement bins. Choose
statistics and ED cutoffs appropriate to the requested calculation, and
calibrate a modest case before scheduling a larger scan. Select enough
measurement blocks to estimate the SEM; record `dtau=beta/Ltrot` explicitly.

The solver runners accept the same stage pattern, with separate manifest
schemas. For an example paired workspace prepared at `runs/custom-pairing/`:

```bash
python3 src/pairing/run_paper.py --manifest runs/custom-pairing/manifest.json --output runs/custom-pairing/results --mode dqmc --seed 24681357 --dry-run
python3 src/pairing/run_paper.py --manifest runs/custom-pairing/manifest.json --output runs/custom-pairing/results --mode dqmc --seed 24681357
python3 src/pairing/run_paper.py --manifest runs/custom-pairing/manifest.json --output runs/custom-pairing/results --mode ed --seed 24681357
python3 src/pairing/run_paper.py --manifest runs/custom-pairing/manifest.json --output runs/custom-pairing/results --mode analyze --seed 24681357
```

For number-conserving runs use `src/number_conserving/run_paper.py`,
the number-conserving manifest, and its input `seeds.txt` (no `--seed` flag).
Add `--python /path/to/quspin/python` to the ED stage if necessary. Use one
numerical-library thread for initial calibration. The executables append to
fixed filenames, so independent chains need fresh output directories.

Deliver the prepared inputs, executable commands, output paths, and processed
means/SEM with reference settings and measured cost. Use named examples and
outputs to keep calculations distinct; update existing inputs or processed
data intentionally when requested and document the reason. Curated small
inputs and summaries can become tracked examples;
large raw measurements remain local and regenerable.
