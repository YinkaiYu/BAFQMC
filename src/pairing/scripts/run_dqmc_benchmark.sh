#!/usr/bin/env bash
set -euo pipefail

if [ "$#" -lt 2 ] || [ "$#" -gt 3 ]; then
  echo "usage: $0 INPUT_DIR REFERENCE_JSON [MPI_NP]" >&2
  exit 2
fi

input_dir=$1
reference=$2
mpi_np=${3:-1}
repo_root=$(cd "$(dirname "$0")/.." && pwd -P)
python_bin=${PYTHON:-python3}
tmp_parent="$repo_root/benchmarks/tmp"

case "$input_dir" in
  /*) ;;
  *) input_dir="$repo_root/$input_dir" ;;
esac
case "$reference" in
  /*) ;;
  *) reference="$repo_root/$reference" ;;
esac

mkdir -p "$tmp_parent"
tmp_run=$(mktemp -d "$tmp_parent/run.XXXXXX")
cp "$input_dir"/paramC_sets.txt "$tmp_run"/
cp "$input_dir"/confin.txt "$tmp_run"/
cp "$input_dir"/seeds.txt "$tmp_run"/

bash "$repo_root/scripts/run_local.sh" "$tmp_run" "$mpi_np"
"$python_bin" "$repo_root/benchmarks/compare.py" --reference "$reference" --run-dir "$tmp_run"
