#!/usr/bin/env bash
set -euo pipefail

run_dir="${1:-../../examples/number_conserving/examples/triangle_3x2}"
np="${2:-${MPI_NP:-1}}"

repo_root="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
exe="${repo_root}/build/bosonDQMC.out"

case "${run_dir}" in
  /*) ;;
  *) run_dir="${repo_root}/${run_dir}" ;;
esac

if [[ ! -d "${run_dir}" ]]; then
  echo "Run directory does not exist: ${run_dir}" >&2
  exit 2
fi

for required in paramC_sets.txt confin.txt seeds.txt; do
  if [[ ! -f "${run_dir}/${required}" ]]; then
    echo "Missing required input file: ${run_dir}/${required}" >&2
    exit 2
  fi
done

if ! command -v mpirun >/dev/null 2>&1; then
  echo "mpirun is not available on PATH" >&2
  exit 2
fi

make -C "${repo_root}" build

cd "${run_dir}"
exec mpirun -np "${np}" "${exe}"
