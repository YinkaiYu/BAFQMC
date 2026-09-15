#!/usr/bin/env bash
set -euo pipefail

if [ "$#" -lt 1 ] || [ "$#" -gt 2 ]; then
  echo "usage: $0 RUN_DIR [MPI_NP]" >&2
  exit 2
fi

run_dir=$1
mpi_np=${2:-1}
repo_root=$(cd "$(dirname "$0")/.." && pwd -P)
target="$repo_root/build/bosonDQMC.out"

case "$run_dir" in
  /*) ;;
  *) run_dir="$repo_root/$run_dir" ;;
esac

for input in paramC_sets.txt confin.txt seeds.txt; do
  if [ ! -f "$run_dir/$input" ]; then
    echo "missing required input: $run_dir/$input" >&2
    exit 1
  fi
done

make -C "$repo_root" build

mpi_launcher=${MPIEXEC:-mpirun}
if ! command -v "$mpi_launcher" >/dev/null 2>&1; then
  echo "MPI launcher not found: $mpi_launcher" >&2
  exit 1
fi

cd "$run_dir"
export OMP_NUM_THREADS=${OMP_NUM_THREADS:-1}
export MKL_NUM_THREADS=${MKL_NUM_THREADS:-1}
export OPENBLAS_NUM_THREADS=${OPENBLAS_NUM_THREADS:-1}
"$mpi_launcher" -np "$mpi_np" "$target"
