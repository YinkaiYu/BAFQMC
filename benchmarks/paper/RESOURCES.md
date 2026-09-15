# Production resources and measured timings

`python3 reproduce.py` runs the complete 22-point BAFQMC/reference calculation
with the paper's parameters, then generates the processed data and figures.
The main-text figure comprises eight number-conserving U points and seven
paired Delta points; the Supplemental Material adds seven U1 points.
The default is one MPI rank, one numerical-library thread, and sequential cases.
Plan for **12–24 hours, 16 GiB RAM with at least 8 GiB available, and 8 GiB of
free disk** on a recent desktop. The time is a planning range assembled from
the measurements below; the whole integrated campaign has not been timed as
one uninterrupted run. A slower processor or filesystem can take longer.

Use `python3 reproduce.py --plan` to display this budget without starting a job.
The same estimate is printed before a production run. `--threads N` changes
the numerical-library thread count; the one-thread timing range should not be
divided by N because the stages scale differently.
Use `--scope main` or `--scope supplement` to select the corresponding figure,
and add `--model` to restrict the solver. The full campaign budget above includes
all three scans.

## Measured evidence

Calibration measurements used an AMD Ryzen 5 9600X, Linux/WSL2, Intel Fortran
2025.2.1/MKL, and one MPI rank/thread. Python ED used QuSpin 1.0.1, NumPy 2.4.4,
SciPy 1.17.1, and Numba 0.65.1.

| Calculation | Measurement | Peak resident memory |
|---|---|---|
| Paired BAFQMC, Delta=0.2, original beta=4 and Delta tau=0.01, 500 warmup iterations, 1000 bins | 46.38 s elapsed | 65 MiB |
| Paired ED, Delta=0.2, nmax=3, ncut=4, full 7297-state trace | 40.40 s elapsed | 3.01 GiB |
| Main-text number-conserving BAFQMC, U=0.25 (input U2), beta=4, original warmup, 1000 bins | 5.47 s elapsed | 61 MiB |
| Number-conserving BAFQMC, U1=-0.6, beta=1, original warmup, 1000 bins | 2.29 s elapsed | 61 MiB |
| Original seven paired 100000-bin production chains | 12.20 h summed case elapsed time | Not recorded |
| Historical supplemental attractive-density ED, six completed shells | 243.5–289.3 s per point | Not recorded |

Multiplying the paired 1000-bin calibration by 100 gives a conservative
single-point estimate of 1.29 h because it also multiplies warmup overhead.
Seven such points give about 9 h. Together with the historical complete runs,
**9–13 h** is a useful budget for the paired BAFQMC stage on comparable hardware.
The seven dense paired ED calculations add approximately **5–10 minutes**.
The two number-conserving calibrations give an approximately **1.3–1.7 hour**
budget for its 15 BAFQMC points. Number-conserving ED has varying block sizes;
the main-text U=0.25, seven-shell point dominates that reference stage.
Its historical successful checkpoint was written about 67.5 minutes after the
recorded start. That historical thread count was not retained. Allow roughly
**3–5 hours for the complete number-conserving part**, including its ED stage.
The recalculated paired ED point reproduced the archived four reference values
to 5.3e-18 absolute accuracy.

## Number-conserving ED memory

The archived interacting references stop after 4–7 particle-number shells.
The largest required momentum block is at U=0.25 (input U1=0, U2=0.25):
dimension 9075 in shell 7.
Four complex128 matrices of this size occupy about 4.91 GiB; sparse operators,
Python, and numerical-library workspace require additional memory. This is why
the full campaign should have at least 8 GiB available, even though BAFQMC
itself uses far less memory.

The free U1=U2=0 point is freshly evaluated with the analytic Bose distribution,
as in the paper. It does not launch an unnecessary truncated many-body ED.
The original convergence policies and particle cutoffs are preserved for the
interacting points. Memory checks stop before an oversized dense solve; they do
not silently lower the particle cutoff.

## Storage and continuation

BAFQMC repeatedly appends to many small files. On WSL, run those writes in the
Linux filesystem rather than `/mnt/c`. The production driver automatically
creates Linux temporary scratch space and records its path in `progress.json`.
Use `--work-dir /path/on/fast/linux/storage` to choose an empty scratch directory.
Each completed case/stage is copied into the selected output directory. Raw
outputs stay outside Git. Scratch and the copied outputs can coexist, so the
disk budget includes both.

Every case/stage has a log, elapsed time, and file hashes. To continue an
interrupted campaign:

```bash
python3 reproduce.py --resume --output benchmarks/paper/output/full-YYYYMMDD-HHMMSS
```

Use the same `--scope`, `--model`, `--python-ed`, and `--threads` as the original command.
Completed stages are checked and skipped. An interrupted chain is moved aside
in scratch and restarted from its original input; its partial samples are never
appended to a new chain. `--mode dqmc`, `--mode ed`, and `--mode analyze` can also
run the production stages separately with the same explicit `--output`.

After all results have been copied and checked, the scratch directory named in
`progress.json` may be removed. Keep the final `runs/`, tables, environment record,
and logs needed for your research.
