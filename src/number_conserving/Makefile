include ../common/compiler.mk

SRC_DIR := src
BUILD_DIR ?= build
OBJ_DIR := $(BUILD_DIR)/obj
MOD_DIR := $(BUILD_DIR)/mod
TARGET := $(BUILD_DIR)/bosonDQMC.out
PYTHON ?= python3
MPI_NP ?= 1
RUN_DIR ?= $(BUILD_DIR)/example
BENCHMARK_RUN_DIR ?=
BENCHMARK_REFERENCE ?= benchmarks/references
DQMC_FAST_BENCHMARK_INPUT_DIR ?= ../../examples/number_conserving/benchmarks/triangle_3x2_free_beta3_mu-2.5
DQMC_FAST_BENCHMARK_REFERENCE ?= benchmarks/dqmc_references/triangle_3x2_free_beta3_mu-2.5.json
DQMC_BENCHMARK_SUITE ?= benchmarks/dqmc_suite.json
DQMC_BENCHMARK_INPUT_DIR ?=
DQMC_BENCHMARK_REFERENCE ?=

SOURCES := \
  calc_basic.f90 \
  lattice.f90 \
  fields.f90 \
  non_interact.f90 \
  operator_Hubbard.f90 \
  model.f90 \
  process_matrix.f90 \
  multiply.f90 \
  obser_equal.f90 \
  obser_tau.f90 \
  stabilization.f90 \
  localU.f90 \
  dynamics.f90 \
  local_sweep.f90 \
  fourier_trans.f90 \
  pole_diagnostics.f90 \
  main.f90

OBJECTS := $(OBJ_DIR)/mymats.o $(OBJ_DIR)/random.o $(addprefix $(OBJ_DIR)/,$(SOURCES:.f90=.o))

.PHONY: all build clean print-config run-example benchmark benchmark-fast benchmark-dqmc check-fixtures benchmark-ed help
.DEFAULT_GOAL := build
all: build
build: $(TARGET)

$(TARGET): $(OBJECTS)
	$(FC) $(LDFLAGS) -o $@ $(OBJECTS) $(LDLIBS)

$(OBJ_DIR)/mymats.o: ../common/mymats.f90 Makefile ../common/compiler.mk | $(OBJ_DIR) $(MOD_DIR)
	$(FC) -c $(FFLAGS) $(MOD_FLAG) -I$(MOD_DIR) -o $@ $<

$(OBJ_DIR)/random.o: ../common/random.f90 Makefile ../common/compiler.mk | $(OBJ_DIR) $(MOD_DIR)
	$(FC) -c $(FFLAGS) $(MOD_FLAG) -I$(MOD_DIR) -o $@ $<

$(OBJ_DIR)/%.o: $(SRC_DIR)/%.f90 $(wildcard $(SRC_DIR)/*.f90) ../common/mymats.f90 ../common/random.f90 Makefile ../common/compiler.mk $(OBJ_DIR)/mymats.o | $(OBJ_DIR) $(MOD_DIR)
	$(FC) -c $(FFLAGS) $(MOD_FLAG) -I$(MOD_DIR) -o $@ $<

$(OBJ_DIR) $(MOD_DIR):
	mkdir -p $@

print-config:
	@echo "FC=$(FC)"
	@echo "COMPILER_VERSION=$(COMPILER_VERSION)"
	@echo "FFLAGS=$(FFLAGS)"
	@echo "LDLIBS=$(LDLIBS)"
	@echo "TARGET=$(TARGET)"

run-example: build
	$(PYTHON) scripts/run_example.py --output $(RUN_DIR) --np $(MPI_NP)

benchmark: benchmark-dqmc
benchmark-fast: build
	PYTHON=$(PYTHON) bash scripts/run_dqmc_benchmark.sh $(DQMC_FAST_BENCHMARK_INPUT_DIR) $(DQMC_FAST_BENCHMARK_REFERENCE) $(MPI_NP)

check-fixtures:
	$(PYTHON) benchmarks/compare.py --reference $(BENCHMARK_REFERENCE) $(if $(strip $(BENCHMARK_RUN_DIR)),--run-dir $(BENCHMARK_RUN_DIR),)

benchmark-dqmc: build
ifneq ($(strip $(DQMC_BENCHMARK_INPUT_DIR)$(DQMC_BENCHMARK_REFERENCE)),)
	PYTHON=$(PYTHON) bash scripts/run_dqmc_benchmark.sh $(DQMC_BENCHMARK_INPUT_DIR) $(DQMC_BENCHMARK_REFERENCE) $(MPI_NP)
else
	$(PYTHON) benchmarks/run_dqmc_suite.py --suite $(DQMC_BENCHMARK_SUITE) --np $(MPI_NP) --python $(PYTHON)
endif

benchmark-ed:
	mkdir -p $(BUILD_DIR)/ed-3x2
	cp benchmarks/ed/params_triangle_3x2.txt $(BUILD_DIR)/ed-3x2/params.txt
	cd $(BUILD_DIR)/ed-3x2 && $(PYTHON) $(abspath benchmarks/ed/EDtriangle_symm_NEblock.py)

clean:
	rm -rf $(BUILD_DIR)

help:
	@echo "Targets: build print-config run-example benchmark-fast benchmark-dqmc check-fixtures benchmark-ed clean"
