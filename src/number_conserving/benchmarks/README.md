# Benchmarks

Commands below run from `src/number_conserving/` and exercise the solver during development. The live DQMC-vs-ED suite is:

```bash
make benchmark
```

`make benchmark` is an alias for `make benchmark-dqmc`. It runs the manifest in `benchmarks/dqmc_suite.json`, which includes the free-boson analytic case plus all four ED cases preserved in `benchmarks/dqmc_references/`. All cases must pass for `total_NE`, `total_kinetic`, `doubleOcc`, `squareOcc`, `numsquare_up`, and `numsquare_do`.

All live benchmark inputs use `dtau = beta / Ltrot = 0.01`. The interacting cases use `Nbin = 100000`, `Nsweep = 1`, `shiftLoc = 1.5`, and warm-up enabled. On the validation workstation under Linux/WSL with `MPI_NP=1`, the full live suite was observed at `real 604.50` seconds, about 10 minutes 5 seconds; budget at least 15 minutes, and longer if the machine is busy.

The fast benchmark is a real DQMC run for the no-interaction `U1=U2=0` case:

```bash
make benchmark-fast
```

It compares against the analytic reference in `benchmarks/dqmc_references/triangle_3x2_free_beta3_mu-2.5.json`. This is useful for quick executable, Green-function, particle-number, kinetic-normalization, and equal-time scalar observable checks, but it does not exercise the interacting Monte Carlo path.

`make check-fixtures` is a fast comparison-script check over checked-in fixture files. It discovers every `*.json` file under `benchmarks/references/`. Each reference points to a fixture DQMC-style output directory under `benchmarks/fixtures/`. This validates parsing and documented DQMC-to-ED normalization semantics, but it is not physics validation and must not replace `make benchmark` for algorithm changes.

For one-off debugging of a single case, override both `DQMC_BENCHMARK_INPUT_DIR` and `DQMC_BENCHMARK_REFERENCE`. Do not use a single-case override as the final benchmark evidence for substantive algorithm changes.

## Regression statistics

For an interacting live DQMC observable, `benchmarks/compare.py` first converts each bin to the ED scale. Live benchmark comparisons cover the `num_up`, `num_do`, `kinetic`, `doubleOcc`, `squareOcc`, `numsquare_up`, and `numsquare_do` conversions documented below. It then groups the 100000 samples into 10 blocks of 10000 samples. The reported `stderr` is the standard error of the Monte Carlo mean:

```text
stderr = std(block_means) / sqrt(number_of_blocks)
z = (DQMC_block_mean - ED_value) / stderr
```

The field `stderr_tolerance` is the allowed number of standard errors of the mean, not the standard deviation of individual Monte Carlo samples. The current references use `stderr_tolerance = 3.0`; a case passes when:

```text
abs(DQMC_block_mean - ED_value) <= max(atol, rtol * abs(ED_value), stderr_tolerance * stderr)
```

The free-boson case is deterministic and uses a tight absolute tolerance instead of block statistics.

The root production reproduction interface exports all valid means, standard
errors, and reference differences without a sigma acceptance threshold.
These configurable tolerances apply to development regression tests.

## Observable Conversions

For the current DQMC output:

```text
total_NE = last(num_up) + last(num_do)
total_kinetic = last(kinetic) * Lq
doubleOcc = last(doubleOcc)
squareOcc = last(squareOcc)
numsquare_up = last(numsquare_up)
numsquare_do = last(numsquare_do)
```

The ED script accumulates `NE = NE_b + NE_c` and its kinetic operator includes both b and c hopping layers. The DQMC `kinetic` file includes both flavors and is divided by `Lq` in `src/obser_equal.f90`. `doubleOcc` is compared per site. `squareOcc` is compared as the per-site half normal-ordered same-flavor onsite pair, `0.5 * <sum_i,s n_s,i (n_s,i - 1)> / Lq`. `numsquare_up` and `numsquare_do` are compared as total flavor number-square observables.

Interacting real DQMC benchmarks use sample means instead of the last bin:

```text
total_NE = mean(num_up) + mean(num_do)
total_kinetic = mean(kinetic) * Lq
doubleOcc = mean(doubleOcc)
squareOcc = mean(squareOcc)
numsquare_up = mean(numsquare_up)
numsquare_do = mean(numsquare_do)
```

## Custom 3x3 campaigns

For the full paper campaign, use the repository-root
[reproduction interface](../../../benchmarks/paper/README.md). It retains the
exact paper parameters, seeds, and ED policies. The separate example manifest
`benchmarks/campaigns/triangle_3x3_observables.json` is a starting point
for new research scans.

Create a working manifest and input directories under ignored `data/`:

```bash
python3 benchmarks/campaign_3x3.py write-manifest --output ../../runs/my-campaign/manifest.json
python3 benchmarks/campaign_3x3.py init-local --manifest ../../runs/my-campaign/manifest.json --output-dir ../../runs/my-campaign/inputs --seed 64464988
```

Each case has `paramC_sets.txt`, `confin.txt`, `seeds.txt`, and ED
`params.json`. Use a new directory for an independent chain.
The portable `run_paper.py` interface also accepts custom case manifests;
see [the solver guide](../README.md) for its schema and stage commands.

The relative-density interaction is `U1`; the relative-density scan uses
`U1=U` with `U2=0`. Both interaction channels remain available for new
calculations. For triangular hopping `t=1`, the single-particle minimum
is -3, so the noninteracting grand-canonical reference requires `mu < -3`.
Attractive `U2` scans use the stated finite-occupation reference prescription.

For lattice dimensions divisible by three, `K=(4*pi/3,0)` has code index
`(2*Lx/3+1,Ly/3+1)`, consistent with
`a1=(1,0)` and `a2=(1/2,sqrt(3)/2)`.

`benchmarks/ed/EDtriangle_quspin_3x3.py` resolves total-particle and
translation blocks, computes thermal expectation values, and saves a
checkpoint after each completed particle shell. The fields
`max_total_particles` and `convergence_policy` specify the shell range
and stopping rule. Use `--no-stop-when-reliable`, or set
`stop_when_reliable` to false, to complete the configured maximum.

The paper inputs retain their specified convergence policies. The attractive
scan uses six completed shells. Stable free points use the analytic Bose
distribution, including `<N^2>=<N>^2+sum_k n_k*(1+n_k)` for each flavor.
Each result records physical parameters, completed shells, and reference status.
See [the resource guide](../../../benchmarks/paper/RESOURCES.md) before
increasing particle cutoffs.

Campaign observables and DQMC file names:

| Observable | DQMC files used | Notes |
| --- | --- | --- |
| `density_total` | `density_total` | per-site total particle density |
| `energy_density` | `energy_density` | physical energy density without `-mu N` |
| `doubleOcc` | `doubleOcc` | per-site cross-flavor onsite density product |
| `num_up` | `num_up` | total `b`-flavor particle number |
| `num_do` | `num_do` | total `c`-flavor particle number |
| `numsquare_up` | `numsquare_up` | total `b`-flavor number-square moment `<N_b^2>` |
| `numsquare_do` | `numsquare_do` | total `c`-flavor number-square moment `<N_c^2>` |
| `IPR` | `density_site_total` | normalized total-density profile `sum_i rho_i^2/(sum_i rho_i)^2`; uniform 3x3 gives `1/9` |
| `S_SF_K` | `sf_K` | complex file, real/imag columns, K point; real part is the plotted structure factor |
| `S_PSF_Gamma` | `psf_Gamma` | complex file, real/imag columns, Gamma point; real part is the plotted structure factor |
| `S_DW_K` | `dw_K` | complex file, real/imag columns, K point; real part is the plotted structure factor |

Optional runtime smoke checks:

```bash
RUN_DQMC_SMOKE=1 python3 -m unittest \
  discover -s ../../tests/number_conserving -p 'test_obser_equal.py' -v

RUN_ED_SMOKE=1 python3 -m unittest \
  discover -s ../../tests/number_conserving -p 'test_ed_quspin_schema.py' -v
```

Use the Python interpreter from your QuSpin environment for ED work. Set
`ED_PYTHON=/path/to/python` to select that interpreter for optional runtime tests.
For paper cases, use `run_paper.py`; it accepts the exact distributed input
files and ED convergence policy through a portable manifest. The campaign
helper prepares broader custom sweeps; cluster scheduling is configured by the
user when needed.

## Pole Diagnostic Checks

Continuous pole diagnostics are structural and analysis outputs, not ED physics observables. They are not used in live benchmark pass/fail decisions.

For a run with `Nbin` bins and MPI size `ISIZE`, the diagnostic files contain `Nbin * ISIZE` configuration samples. One sample is written per bin per MPI rank, and rank 0 gathers rank-local samples and writes rows in rank order. These samples are not rank averaged.

Validate the diagnostic file structure and consistency with:

```bash
python3 benchmarks/check_pole_diagnostics.py <run_dir>
```

## Adding A Case

Add one JSON file under `benchmarks/references/` and one fixture directory under `benchmarks/fixtures/`. The reference JSON must include:

- `case`: stable case identifier, matching the JSON filename stem
- `dqmc_fixture`: repository-relative path to the fixture output directory
- `parameters`: `Lx`, `Ly`, `t`, `U1`, `U2`, `beta`, `mu`
- `observables`: the scalar observables needed by the case, with explicit operation names such as `last`, `mean`, `sum_mean`, or `mean_times_lq`

The fixture directory must contain scalar output files with DQMC names, currently `num_up`, `num_do`, `kinetic`, `doubleOcc`, `squareOcc`, `numsquare_up`, and `numsquare_do`.

Run the suite after adding a case:

```bash
python3 -m unittest discover -s ../../tests/number_conserving -p 'test_compare.py' -v
make check-fixtures
```

If the case is meant to exercise the live DQMC executable, add a JSON under `benchmarks/dqmc_references/`, add or reuse an input directory under `../../examples/number_conserving/benchmarks/`, wire it through `benchmarks/dqmc_suite.json`, and run `make benchmark`.

## Optional ED Recompute

Run:

```bash
make benchmark-ed
```

This runs `benchmarks/ed/EDtriangle_symm_NEblock.py` from inside `benchmarks/ed/`. It requires QuSpin and Numba and can take much longer than the default benchmark.
