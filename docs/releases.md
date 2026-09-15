# Source releases

Each release identifies a computational source snapshot, its documented
benchmark inputs, and its validation results. Use a versioned Git tag and
GitHub release so researchers can return to the same code.

## Prepare the version

Update [CHANGELOG.md](../CHANGELOG.md) with the release date and concrete
changes. Keep the installation guide, agent instructions, and numerical
resource estimates consistent with the code. Check that the source package
contains the MIT license, solver source, ED code, curated inputs, processed
data, and reproduction scripts. Raw chains and generated build files remain
outside the release.

Run these checks in the documented Linux environment:

```bash
python3 scripts/doctor.py
make check
make physics
make smoke
python3 reproduce.py --mode check
git diff --check
```

The GNU Fortran GitHub Actions checks should also pass for the release commit.
When production orchestration changes, run the real continuation tests:

```bash
BAFQMC_RUN_MPI_TESTS=1 python3 -m unittest discover -s tests/reproduction -v
```

For physics changes, compare the Hamiltonian, estimator normalization, Trotter
step, sampling statistics, and ED occupation cutoffs with the previous
release. Run the corresponding small-system identities and live regression
calculations described in [the development guide](development.md). Keep
published benchmark inputs and processed data intact unless making an explicit,
documented update to the model, input settings, or reference calculation.

The default `python3 reproduce.py` command computes all 22 production points.
Its budget is approximately **12–24 hours, 16 GiB RAM, and 8 GiB free disk**
on the reference desktop; see [the measured resource guide](../benchmarks/paper/RESOURCES.md).
When a release changes production physics or benchmark inputs, rerun the
affected production scans and record their parameters and results. Distinguish
the completed checks from the full calculations in the release notes.

## Publish the source snapshot

Commit the completed changes and inspect `git status --short`. Tag the exact
validated commit with the new version, for example:

```bash
release_tag=v0.1.3
git tag -a "$release_tag" -m "BAFQMC ${release_tag}"
git push origin main
git push origin "$release_tag"
```

The Source release workflow reruns the computational CI for that tag, then
publishes a source archive and SHA-256 checksum as a GitHub release. Include
the changelog entry, relevant validation, and any changes to input or output
formats in the release notes. Keep the tag fixed once published; a correction
receives a new version.

GitHub supplies source archives for the tagged snapshot. Download one to a
fresh directory and check `python3 reproduce.py --plan` and
`python3 reproduce.py --mode check` before treating the release as complete.
The extracted package must work without a manuscript checkout or sibling
repository. The optional development container provides the same GNU/MPI
environment; opening it runs the environment doctor and leaves production
calculations under the researcher's control.
