# Portable GNU MPI + BLAS/LAPACK default; Intel MPI/MKL is also supported.
# Source Intel's environment setup before building with Intel tools.
ifeq ($(filter command% environment,$(origin FC)),)
  ifneq ($(shell command -v mpifort 2>/dev/null),)
    FC := mpifort
  else ifneq ($(shell command -v mpiifx 2>/dev/null),)
    FC := mpiifx
  else ifneq ($(shell command -v mpiifort 2>/dev/null),)
    FC := mpiifort -fc=ifx
  else
    FC := mpifort
  endif
endif

COMPILER_VERSION := $(shell $(FC) --version 2>/dev/null | head -n 1)
ifneq ($(findstring GNU,$(COMPILER_VERSION)),)
  FFLAGS ?= -O3 -ffree-line-length-none -fallow-argument-mismatch
  MOD_FLAG = -J$(MOD_DIR)
  LDFLAGS ?=
  LDLIBS ?= -llapack -lblas
else
  FFLAGS ?= -O3
  MOD_FLAG = -module $(MOD_DIR)
  LDFLAGS ?=
  LDLIBS ?= -qmkl
endif
# The original source uses ordered Fortran module dependencies.
.NOTPARALLEL:
