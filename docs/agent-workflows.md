# Research workflows for agents

Read [AGENTS.md](../AGENTS.md) first. Commands here run from the repository
root in Linux or WSL. The repository provides working triangular-lattice
solvers and editable campaign inputs; new geometries require implementation
work rather than a generic lattice-selection flag.

The portable task skills in `.agents/skills/` route these workflows:
[reproduce benchmarks](../.agents/skills/bafqmc-reproduce/SKILL.md),
[run a new calculation](../.agents/skills/bafqmc-new-calculation/SKILL.md), and
[add an observable](../.agents/skills/bafqmc-add-observable/SKILL.md).

## Set up and establish a working baseline

1. Follow [environment setup](getting-started.md), then run
   `python3 scripts/doctor.py` and `make check`.
2. Run `python3 reproduce.py --mode smoke` to exercise both BAFQMC and ED.
3. Inspect the run logs and result tables. Report the executed cases and output
   directory to the user, including which optional checks were skipped.

With separate analysis and ED environments:

```bash
python3 scripts/doctor.py --python-ed /path/to/quspin/python
python3 reproduce.py --mode smoke --python-ed /path/to/quspin/python
```

## Reproduce the full paper

Run `python3 reproduce.py --plan`, read
[RESOURCES.md](../benchmarks/paper/RESOURCES.md), and relay the computing budget.
When the user has asked for the full reproduction, carry out the calculation:

```bash
python3 reproduce.py --output benchmarks/paper/output/paper-run
```

The default includes all 22 production points, freshly generated references,
statistical processing, and figures. Preserve the original sampling and ED
settings. If execution is interrupted, continue with:

```bash
python3 reproduce.py --output benchmarks/paper/output/paper-run --resume
```

Use the same scope, model selection, interpreter, and numerical-library
thread setting. Completed stages are verified before reuse; partial stages
restart from the initial inputs. Inspect `progress.json`, `observables.csv`,
and both generated figures before reporting completion.

## Prepare a number-conserving research point

Use a separate manifest and inputs. This preparation snippet copies the tiny
free-boson example into an ignored research workspace:

```bash
python3 - <<'PY'
import json
from pathlib import Path
import shutil

root = Path.cwd()
source = root / 'examples/number_conserving/examples/triangle_3x3_free'
work = root / 'runs/custom-number'
inputs = work / 'inputs/my-point'
inputs.mkdir(parents=True)
for name in ('paramC_sets.txt', 'confin.txt', 'seeds.txt', 'params.json'):
    shutil.copy2(source / name, inputs / name)
parameters = json.loads((inputs / 'params.json').read_text())['parameters']
manifest = {
    'cases': [{
        'id': 'my-point',
        'input_dir': 'inputs/my-point',
        'ed_params': 'inputs/my-point/params.json',
        'parameters': parameters,
        'analysis': {'block_size': 2, 'skip_samples': 0},
    }]
}
(work / 'manifest.json').write_text(json.dumps(manifest, indent=2) + '\n')
PY
```

This is an eight-bin installation-scale starting point. To change the
physics, edit `paramC_sets.txt`, the `parameters` object in `params.json`, and
the case's `parameters` in the manifest consistently. The first numeric row
is `U1 U2 mu`; the second is `Lx Ly Ltrot beta`. The fourth row includes
`Nwrap Nbin Nsweep shiftLoc`; `dtau=beta/Ltrot`. Copy appropriate production
warmup settings when moving beyond the installation example. ED's particle
shells and convergence policy are in `params.json`. Choose a block size that
leaves enough blocks to estimate the SEM from the requested sample count.

Inspect the commands, then execute all three stages against the same output:

```bash
python3 src/number_conserving/run_paper.py --manifest runs/custom-number/manifest.json --output runs/custom-number/results --mode dqmc --dry-run
python3 src/number_conserving/run_paper.py --manifest runs/custom-number/manifest.json --output runs/custom-number/results --mode dqmc
python3 src/number_conserving/run_paper.py --manifest runs/custom-number/manifest.json --output runs/custom-number/results --mode ed
python3 src/number_conserving/run_paper.py --manifest runs/custom-number/manifest.json --output runs/custom-number/results --mode analyze
```

Use `--python /path/to/quspin/python` for ED in a separate environment. A
manifest can contain multiple cases with their own inputs. `--case ID` is
repeatable and selects a subset. Output includes `my-point/dqmc/`,
`my-point/ed/results.json`, and `analysis.json`.

For a Monte Carlo calculation without ED, the direct input-directory runner
is also available; see the
[number-conserving solver guide](../src/number_conserving/README.md).

## Prepare a paired research scan

The pairing manifest is self-contained: its parameters generate both solver
and ED input files. Copy the small pipeline example:

```bash
mkdir -p runs/custom-pairing
cp src/pairing/benchmarks/campaigns/triangle_pairing_pipeline_smoke.json runs/custom-pairing/manifest.json
```

Edit `parameters` for `Lx`, `Ly`, `U1`, `U2`, `mu`, and `beta`; edit `deltas`
for the scan. The current hopping remains `t=1`. `dqmc_defaults` contains
`dtau`, measurement bins, warmup, proposal widths, and blocking; `ed_defaults`
contains `nmax` and `ncut`. Set `analysis_defaults.block_size` consistently
with the desired analysis. The copied example starts at eight bins with a
small ED basis. The published
[pairing manifest](../benchmarks/paper/data/pairing/manifest.json) provides
the complete production settings as a reference.

```bash
python3 src/pairing/run_paper.py --manifest runs/custom-pairing/manifest.json --output runs/custom-pairing/results --mode init --seed 24681357 --dry-run
python3 src/pairing/run_paper.py --manifest runs/custom-pairing/manifest.json --output runs/custom-pairing/results --mode init --seed 24681357
python3 src/pairing/run_paper.py --manifest runs/custom-pairing/manifest.json --output runs/custom-pairing/results --mode dqmc --seed 24681357
python3 src/pairing/run_paper.py --manifest runs/custom-pairing/manifest.json --output runs/custom-pairing/results --mode ed --seed 24681357
python3 src/pairing/run_paper.py --manifest runs/custom-pairing/manifest.json --output runs/custom-pairing/results --mode analyze --seed 24681357
```

Add `--python /path/to/quspin/python` for a separate ED interpreter. The
generated case names appear in the dry run. Use repeatable `--case` arguments
to select cases while retaining the full manifest and its order: the order
sets each case's initial seed. Output includes `inputs/<case>/`,
`ed_results/<case>.json`, and `summary/comparison_observables.csv`.

Run one modest case first and measure time and memory before scaling the
requested research scan. Choose a new output directory for an independent
chain. The low-level executables append measurements to fixed filenames.

## Add an observable

1. Write its operator definition and units. Specify equal-time order, flavor
   sums, normalization, momentum, and connected subtraction if present.
2. Read the solver's `docs/solvers/<solver>/physics.md` and derive the Wick estimator. Preserve
   the bosonic identity contribution in `<b_i b_j^+>` and the full Nambu
   block convention for pairing.
3. Extend `src/obser_equal.f90` and its output, the corresponding ED observable,
   and the analysis parser/schema. Keep the definition beside the implementation.
4. Verify a small case with an analytic limit, ED identity, or independent ED
   representation; then run an appropriate live BAFQMC comparison.
5. Provide a reproducible input and report mean, SEM, and reference value.
   Track compact fixtures where useful, keeping large raw output local.

Use [the development guide](development.md) for exact check commands.
For a new physical model, coordinate the Hamiltonian, HS factorization,
symmetry conditions, geometry, ED construction, and estimators before extending
the production workflow.
