# Contributing

Contributions to BAFQMC are welcome: new observables, algorithm improvements,
model extensions, computational workflows, and documentation. Contributions
are distributed under the repository's [MIT license](LICENSE).

Start with [AGENTS.md](AGENTS.md) and the
[development guide](docs/development.md). For a physical change, state the
Hamiltonian or estimator being implemented and its connection to the existing
conventions. Use a separate campaign for new research parameters and preserve
the published data package.

A useful pull request explains the problem, the resulting behavior, and the
checks performed. Include enough information to repeat a numerical check:
system size, couplings, temperature, Trotter step, sampling, seeds, and ED
cutoffs where relevant. Share small processed results; keep regenerable raw
chains and local build output out of Git.

Run targeted tests for the change and `git diff --check`. Changes to a
physical estimator should include a small-system comparison or identity test.
Documentation-only changes do not need a full production campaign.

Agent-assisted contributions follow the same scientific and computational
standards. Review the resulting code, data conventions, and reported tests
before opening a pull request.
