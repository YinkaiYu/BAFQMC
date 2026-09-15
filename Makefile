PYTHON ?= python3
PYTHON_ED ?= $(PYTHON)
SCOPE ?= all
MODEL ?= both
.DEFAULT_GOAL := help
.PHONY: help doctor plan reproduce build check physics numerics smoke docs

help:
	@echo 'Targets: doctor plan reproduce build check physics numerics smoke docs'
	@echo 'reproduce runs all 22 benchmark points (12-24 hours); SCOPE=main or supplement selects a subset.'
	@echo 'Set PYTHON and PYTHON_ED for separate analysis and QuSpin environments.'

doctor:
	$(PYTHON) scripts/doctor.py --python-ed "$(PYTHON_ED)"

plan:
	$(PYTHON) reproduce.py --plan --scope $(SCOPE) --model $(MODEL)

reproduce:
	$(PYTHON) reproduce.py --scope $(SCOPE) --model $(MODEL) --python-ed "$(PYTHON_ED)"

build:
	$(MAKE) -C src/number_conserving build
	$(MAKE) -C src/pairing build

check:
	$(PYTHON) -m unittest discover -s tests/reproduction -v
	$(PYTHON) -m unittest discover -s tests/number_conserving -p 'test_*.py' -v
	$(PYTHON) -m unittest discover -s tests/pairing -p 'test_*.py' -v
	$(PYTHON_ED) -m unittest discover -s tests/pairing/ed -p 'test_*.py' -v

physics:
	$(MAKE) numerics
	$(PYTHON_ED) -m unittest discover -s tests/physics -p 'test_ed_physics.py' -v
	$(MAKE) -C src/pairing benchmark-ed QU_SPIN_PYTHON="$(PYTHON_ED)"
	BAFQMC_RUN_MPI_TESTS=1 $(PYTHON) -m unittest discover -s tests/physics -p 'test_*solver*.py' -v

docs:
	$(PYTHON) scripts/build_docs.py

numerics:
	$(MAKE) -C src/common check

smoke:
	$(PYTHON) reproduce.py --mode smoke --python-ed "$(PYTHON_ED)" --output benchmarks/paper/output/smoke-$(shell date +%Y%m%d-%H%M%S)
