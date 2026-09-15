---
name: bafqmc-new-calculation
description: Prepare and run a new triangular-lattice BAFQMC point or parameter scan, with consistent solver and ED inputs and a stated computation budget.
---

# Start a research calculation

Read [AGENTS.md](../../../AGENTS.md), [the algorithm guide](../../../docs/algorithm.md),
and the applicable section of [the research recipes](../../../docs/agent-workflows.md).
Commands run from the repository root in Linux or WSL.

Establish the Hamiltonian, lattice size, temperature, observables, requested
precision, and available resources from the user's task. Make routine setup
choices independently; resolve missing physical assumptions before computing
a different model. A new geometry or hopping pattern requires an implementation
change: the current solvers use periodic triangular hopping with `RT=1` in
`src/calc_basic.f90`, so changing a JSON `t` value alone is insufficient.

Use `src/number_conserving/` for zero pairing, or `src/pairing/`
for onsite pairing. The main-text coupling maps to `U1=0, U2=U`.
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
means/SEM with reference settings and measured cost. Keep original paper
inputs intact. Curated small inputs and summaries can become tracked examples;
large raw measurements remain local and regenerable.
