#!/usr/bin/env bash
set -euo pipefail

input_dir="${1:-../../examples/number_conserving/benchmarks/triangle_3x2_beta3_mu-2.5_u1_1_u2_0}"
reference="${2:-benchmarks/dqmc_references/triangle_3x2_beta3_mu-2.5_u1_1_u2_0.json}"
np="${3:-${MPI_NP:-1}}"
python="${PYTHON:-python3}"

repo_root="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"

case "${input_dir}" in
  /*) ;;
  *) input_dir="${repo_root}/${input_dir}" ;;
esac

case "${reference}" in
  /*) ;;
  *) reference="${repo_root}/${reference}" ;;
esac

if [[ ! -d "${input_dir}" ]]; then
  echo "Benchmark input directory does not exist: ${input_dir}" >&2
  exit 2
fi

run_dir="$(mktemp -d "${TMPDIR:-/tmp}/bosonDQMC-benchmark.XXXXXX")"

for required in paramC_sets.txt confin.txt seeds.txt; do
  if [[ ! -f "${input_dir}/${required}" ]]; then
    echo "Missing required input file: ${input_dir}/${required}" >&2
    exit 2
  fi
  cp "${input_dir}/${required}" "${run_dir}/${required}"
done

echo "DQMC benchmark run directory: ${run_dir}"
bash "${repo_root}/scripts/run_local.sh" "${run_dir}" "${np}"
cd "${repo_root}"
exec "${python}" benchmarks/compare.py --reference "${reference}" --run-dir "${run_dir}"
