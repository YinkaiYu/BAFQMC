---
name: bafqmc-reproduce
description: Reproduce BAFQMC paper benchmarks, check an installation, or resume an interrupted paper calculation using the repository production workflow.
---

# Reproduce the benchmarks

Read [AGENTS.md](../../../AGENTS.md) and run commands from the repository root
in Linux or WSL. Use [the installation guide](../../../docs/getting-started.md)
for dependencies and [the resource measurements](../../../benchmarks/paper/RESOURCES.md)
for planning.

1. Establish the requested scope. The default is all **22 production points**:
   15 main-text points and 7 supplemental points. `--scope main` or
   `--scope supplement` selects a section; `--model number_conserving` or
   `--model pairing` intersects that selection.
2. Inspect the environment and relay the cost before starting the requested
   calculation. The reference desktop budget is approximately **12–24 hours,
   16 GiB RAM, and 8 GiB free disk** for all points.

   ```bash
   python3 scripts/doctor.py
   python3 reproduce.py --plan
   ```

3. For installation checks, use a fresh output directory:

   ```bash
   python3 reproduce.py --mode smoke --output runs/install-check
   ```

   For a requested full reproduction, execute:

   ```bash
   python3 reproduce.py --output runs/paper-reproduction
   ```

   This generates BAFQMC chains and ED/analytic references, then computes
   statistics and figures. `--mode plot` redraws archived data and does not
   satisfy a request to reproduce the calculations. Documentation or routine
   setup work does not require a full production run.

4. To continue an interrupted production run:

   ```bash
   python3 reproduce.py --output runs/paper-reproduction --resume
   ```

   Retain the original scope, model, thread setting, and ED interpreter.
   Completed stages are verified before reuse; an interrupted stage restarts
   from its initial inputs. Inspect `progress.json` and the failed stage log
   before deciding whether an environment or numerical failure needs a fix.

For separate Python environments, add `--python-ed /path/to/quspin/python`
to the doctor, smoke, and production commands. The root runner sets one
numerical-library thread by default and uses Linux temporary scratch.

Inspect `observables.csv`, `records.json`, `environment.json`, and `figures/`
under the selected output. Report executed points, mean/SEM, references,
resource use, and paths. Preserve seeds, warmup, Trotter settings, blocking,
ED cutoffs, and all valid comparison points. Numerical residuals are
diagnostics; there is no mandatory three-SEM acceptance gate.

Keep raw chains and scratch in ignored output directories. The compact
published package in `benchmarks/paper/data/` remains the reference dataset;
new runs do not replace its recorded means.
