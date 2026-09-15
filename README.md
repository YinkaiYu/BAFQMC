<div align="center">

# BAFQMC

### From local positivity to global symmetry

**Bosonic auxiliary-field quantum Monte Carlo**

**English** · [简体中文](README_zh-CN.md)

[![Benchmarks](https://github.com/YinkaiYu/BAFQMC/actions/workflows/reproduce-benchmarks.yml/badge.svg)](https://github.com/YinkaiYu/BAFQMC/actions/workflows/reproduce-benchmarks.yml)
[![License: MIT](https://img.shields.io/badge/License-MIT-417A69.svg)](LICENSE)
[![arXiv: 2609.XXXXX](https://img.shields.io/badge/arXiv-2609.XXXXX-B31B1B.svg)](#paper)
[![Fortran + Python](https://img.shields.io/badge/Fortran%20%2B%20Python-535B86.svg)](docs/development.md)
[![Agent ready](https://img.shields.io/badge/Agent-ready-927043.svg)](AGENTS.md)

[Documentation](https://www.yykspace.com/BAFQMC/) · [Get started](docs/getting-started.md) · [The algorithm](docs/algorithm.md) · [Reproduce results](benchmarks/paper/README.md) · [Guide for agents](AGENTS.md)

</div>

Worldline and stochastic-series-expansion methods organize sign-free sampling
around local matrix-element positivity. **BAFQMC establishes nonnegative bosonic
weights through global symmetry after Hubbard–Stratonovich decoupling.**
Reflection positivity (RP) and time-reversal symmetry (TRS) open a new
computable region: frustrated bosonic models protected by these symmetries.

<p align="center">
  <img src="docs/assets/territory.svg" alt="Sign-problem-free territory classified by frustration and RP/TRS: BAFQMC reaches the frustrated, symmetry-protected region beyond the unfrustrated region of worldline and SSE methods." width="760">
</p>

This repository puts the construction to work: finite-temperature solvers for
two-flavor bosons on a triangular lattice, with and without onsite pairing;
exact diagonalization (ED); and a complete workflow to reproduce the paper's
benchmarks. **Run new calculations. Extend the method. Build on the code.**
Everything here is available under [MIT](LICENSE).

## Work with your agent

Clone this repository and open it in your preferred coding agent. The
[agent instructions](AGENTS.md) and [task recipes](docs/agent-workflows.md)
provide the implementation map, physical conventions, commands, and checks.
Start with a request such as:

> Read AGENTS.md, set up BAFQMC in Linux or WSL, and run the small BAFQMC + ED
> installation check. Tell me where the results are.

> Reproduce all paper benchmarks. Use the published parameters and seeds,
> save the new results in a dedicated directory, and show me the final figures.

> I want to study the triangular-lattice model at U = 1, beta = 4 and mu = -5
> while varying the pairing strength. Prepare a separate campaign, check one
> small case with ED, and explain the computing budget before production.

> Add an observable for my research. Derive its estimator using the existing
> Green-function convention, implement it in BAFQMC and ED, and verify it on a
> small system.

Describe the physics you want to study and the computing resources you can use.
Your agent can work through the details in the linked guides.

## Reproduce the paper

After [environment setup](docs/getting-started.md), run:

```bash
python3 reproduce.py
```

This builds the solvers, runs fresh BAFQMC and ED calculations for **all 22
benchmark points**, processes the measurements, and produces both figures.
Allow approximately **12–24 hours**, **16 GiB RAM**, and **8 GiB free disk**
on a modern desktop CPU. The default uses one MPI rank and one numerical-library
thread, with cases run sequentially. The [resource guide](benchmarks/paper/RESOURCES.md)
gives measured timings and instructions for resuming a campaign.

```bash
python3 reproduce.py --plan        # Show the scope and resource estimate
python3 reproduce.py --mode smoke  # Small BAFQMC + ED installation check
```

Small processed benchmark data, including means and standard errors, are included.
Large simulation outputs are generated locally and excluded from Git.

## Explore

- [Get started](docs/getting-started.md): install, run, and find your results.
- [Algorithm and conventions](docs/algorithm.md): models, symmetry, and observables.
- [Research tasks for agents](docs/agent-workflows.md): new calculations and extensions.
- [Benchmark details](benchmarks/paper/README.md): parameters, data, and reproduction modes.
- [Contributing](CONTRIBUTING.md): develop and validate changes.

## Paper

arXiv: **2609.XXXXX** (forthcoming).
