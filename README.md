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
  <img src="docs/assets/territory.svg" alt="Schematic of sign-problem-free bosonic quantum Monte Carlo: BAFQMC covers the frustrated region satisfying TRS/RP criteria, while WLQMC/SSE covers the unfrustrated region; the upper-left frustrated symmetric sector is the new sign-free regime." width="760">
</p>

This repository puts the construction to work: finite-temperature solvers for
two-flavor bosons on a triangular lattice, with and without onsite pairing;
exact diagonalization (ED); and a complete workflow to reproduce the paper's
benchmarks. **Run new calculations. Extend the method. Build on the code.**
Everything here is available under [MIT](LICENSE).

## What can I simulate?

The supplied solvers calculate finite-temperature properties of **two boson
flavors on a periodic triangular lattice**. In the paper's notation, the
implemented Hamiltonian is

```math
\begin{aligned}
\hat H ={}& t\sum_{\langle ij\rangle}
\left(\hat b_i^+\hat b_j+\hat b_j^+\hat b_i
+\hat c_i^+\hat c_j+\hat c_j^+\hat c_i\right)\\
&-\sum_i\left(\Delta\,\hat b_i^+\hat c_i^+
+\Delta^*\,\hat b_i\hat c_i\right)\\
&+\sum_i\left[U_1\left(\hat n_{b,i}-\hat n_{c,i}\right)^2
+U_2\left(\hat n_{b,i}+\hat n_{c,i}\right)^2\right].
\end{aligned}
```

Here $`U_1\geq0`$ is the repulsive relative-density coupling and
$`U_2\leq0`$ is the attractive total-density coupling. Here
$`\hat n_{b,i}=\hat b_i^+\hat b_i`$ and
$`\hat n_{c,i}=\hat c_i^+\hat c_i`$. Simulations use the grand-canonical
ensemble $`Z=\mathrm{Tr}e^{-\beta(\hat H-\mu\hat N)}`$, with
$`\hat N=\sum_i(\hat n_{b,i}+\hat n_{c,i})`$.
**Positive $`t=1`$ gives frustrated hopping.** The number-conserving solver
sets $`\Delta=0`$; the Nambu solver supports real onsite pairing.

The lattice has $`N_s=L_xL_y`$ sites and primitive vectors

```math
\mathbf a_1=(1,0),\qquad
\mathbf a_2=\left(\frac12,\frac{\sqrt3}{2}\right).
```

Nearest-neighbor bonds run along $`\mathbf a_1`$, $`\mathbf a_2`$, and
$`\mathbf a_2-\mathbf a_1`$, including their reverse directions and periodic
images. The benchmarks use $`3\times3`$ clusters and measure $`\rho`$, $`-E`$,
$`S_{\mathrm{SF}}(K)`$, and $`S_{\mathrm{DW}}(K)`$.

You can set the two channels independently: $`U_2=0`$ gives the
relative-density-only case, while fixing $`U_1=1`$ and varying $`U_2`$ probes
the total-density channel.

The [model, lattice, and parameter guide](docs/algorithm.md) explains the
ensemble, HS fields, and reference calculations.

**Want to study another model?** The [model-development guide](docs/model-development.md)
and dedicated agent skill explain how to extend the lattice, hopping,
interactions, and pairing, with matching ED and physics tests.

## Work with your agent

If you use a coding agent, give it a request such as one of the following. The
[agent instructions](AGENTS.md) and [task recipes](docs/agent-workflows.md)
provide the implementation map and verification steps.

```text
Get https://github.com/YinkaiYu/BAFQMC, read AGENTS.md, check the Linux or WSL
prerequisites, build both solvers, and run the small installation tests. Report
the setup path and every check result.
```

```text
Reproduce every paper benchmark from the default BAFQMC command. Read
AGENTS.md, keep
large raw chains outside Git, retain the processed benchmark data, and report
the resource budget, output paths, and final comparison tables and figures.
```

```text
Get https://github.com/YinkaiYu/BAFQMC and reproduce all paper benchmarks with
the published parameters and seeds. Save new results in a dedicated directory
and show me the final figures.
```

```text
I want to study the triangular-lattice model at U1 = 1, U2 = 0, beta = 4 and
mu = -5 while varying the pairing strength. Prepare a separate campaign, check one
small case with ED, and explain the computing budget before production.
```

```text
Extend the implemented two-channel model to a nearest-neighbor kagome lattice at
t = 1, U1 = 1, U2 = 0, Delta = 0, beta = 4, mu = -5. Read the new-model skill, derive
the HS symmetry, implement the geometry and matching ED and observables,
and validate a small case. Report the computing budget before production.
```

```text
Add an observable for my research. Derive its estimator using the existing
Green-function convention, implement it in BAFQMC and ED, and verify it on a
small system.
```

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
- [Output observables](docs/observables.md): operator definitions and file conventions.
- [Extend to other models](docs/model-development.md): from Hamiltonian to tested implementation.
- [Research tasks for agents](docs/agent-workflows.md): new calculations and extensions.
- [Benchmark details](benchmarks/paper/README.md): parameters, data, and reproduction modes.
- [Contributing](CONTRIBUTING.md): develop and validate changes.

## Paper

arXiv: **2609.XXXXX** (forthcoming).
