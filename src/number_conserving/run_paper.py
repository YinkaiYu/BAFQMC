#!/usr/bin/env python3
"""Run number-conserving paper cases from exact distributed inputs.

Manifest paths are relative to the manifest file. Generated files are isolated in
OUTPUT/CASE/dqmc and OUTPUT/CASE/ed. Existing nonempty run directories are refused.
"""
from __future__ import annotations

import argparse
import json
import math
import os
from pathlib import Path
import re
import shlex
import shutil
import subprocess
import sys

ROOT = Path(__file__).resolve().parent
sys.path.insert(0, str(ROOT.parents[1]))
from src.number_conserving.benchmarks.campaign_analysis import compare_dqmc_ed_case, read_scalar_series, read_complex_series, read_vector_series, free_boson_reference_result
from src.number_conserving.benchmarks.ed.EDtriangle_quspin_3x3 import ED_OBSERVABLES

INPUTS = ('paramC_sets.txt', 'confin.txt', 'seeds.txt')


def resolve(base: Path, value: str) -> Path:
    path = Path(value)
    return path.resolve() if path.is_absolute() else (base / path).resolve()


def empty_output(path: Path) -> None:
    if path.exists() and any(path.iterdir()):
        raise ValueError(f'refusing to mix results in nonempty directory: {path}; choose a fresh --output')
    path.mkdir(parents=True, exist_ok=True)


def json_clean(value):
    if isinstance(value, float) and not math.isfinite(value):
        return None
    if isinstance(value, dict):
        return {key: json_clean(item) for key, item in value.items()}
    if isinstance(value, (list, tuple)):
        return [json_clean(item) for item in value]
    return value


def validate_completed_output(mode: str, directory: Path) -> None:
    """Check finite complete samples/schema, independently of statistical agreement."""
    if mode == 'dqmc':
        rows = (directory / 'paramC_sets.txt').read_text().splitlines()
        expected = int(rows[3].split()[1])
        lx, ly = map(int, rows[1].split()[:2])
        series = {
            name: read_scalar_series(directory / name)
            for name in ('density_total', 'energy_density', 'doubleOcc', 'num_up',
                         'num_do', 'onsite_n2_up', 'onsite_n2_do')
        }
        series.update({name: read_complex_series(directory / name)
                       for name in ('sf_K', 'psf_Gamma', 'dw_K')})
        series['density_site_total'] = read_vector_series(directory / 'density_site_total', lx*ly)
        for name, values in series.items():
            if len(values) != expected:
                raise ValueError(f'{directory / name}: expected {expected} bins, found {len(values)}')
    else:
        payload = json.loads((directory / 'results.json').read_text())
        if payload.get('status') not in {'converged', 'low_density_cutoff_accepted', 'incomplete', 'exact_free_boson'}:
            raise ValueError(f'{directory}: ED did not produce a physical result')
        for name in ED_OBSERVABLES:
            if not math.isfinite(float(payload['observables'][name])):
                raise ValueError(f'{directory}: nonfinite ED observable {name}')


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--manifest', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    parser.add_argument('--mode', choices=('dqmc', 'ed', 'analyze'), required=True)
    parser.add_argument('--case', action='append', default=[], help='run only this case ID; repeatable')
    parser.add_argument('--np', type=int, default=1, help='MPI ranks (paper production used 1)')
    parser.add_argument('--python', default=sys.executable, help='Python with QuSpin for ED')
    parser.add_argument('--dense-memory-cap-gib', type=float, default=12.0, help='ED dense-memory budget; <=0 disables guard')
    parser.add_argument('--dry-run', action='store_true', help='validate inputs and print commands without running')
    args = parser.parse_args()
    try:
        manifest_path = args.manifest.resolve()
        manifest = json.loads(manifest_path.read_text())
        cases = manifest['cases']
        if not isinstance(cases, list) or not cases:
            raise ValueError('manifest must define a nonempty cases list')
        ids = [case['id'] for case in cases]
        if len(ids) != len(set(ids)) or any(not re.fullmatch(r'[A-Za-z0-9_.-]+', name) or name in {'.', '..'} for name in ids):
            raise ValueError('case IDs must be unique simple directory names')
        if set(args.case) - set(ids):
            raise ValueError(f'unknown case IDs: {sorted(set(args.case)-set(ids))}')
        cases = [case for case in cases if not args.case or case['id'] in args.case]
        if args.np < 1:
            raise ValueError('--np must be positive')
        output = args.output.resolve()
        env = os.environ.copy()
        for key in ('OMP_NUM_THREADS', 'MKL_NUM_THREADS', 'OPENBLAS_NUM_THREADS'):
            env.setdefault(key, '1')
        prepared = []
        for case in cases:
            inputs = resolve(manifest_path.parent, case['input_dir'])
            ed_params = resolve(manifest_path.parent, case.get('ed_params', str(inputs / 'params.json')))
            if args.mode == 'dqmc':
                for name in INPUTS:
                    if not (inputs / name).is_file():
                        raise FileNotFoundError(inputs / name)
                target = output / case['id'] / 'dqmc'
                command = ['mpirun', '-np', str(args.np), str(ROOT / 'build' / 'bosonDQMC.out')]
            elif args.mode == 'ed':
                if not ed_params.is_file():
                    raise FileNotFoundError(ed_params)
                target = output / case['id'] / 'ed'
                parameters = json.loads(ed_params.read_text())['parameters']
                is_free = float(parameters['U1']) == 0.0 and float(parameters['U2']) == 0.0
                command = [] if is_free else [
                    args.python, str(ROOT / 'benchmarks' / 'ed' / 'EDtriangle_quspin_3x3.py'),
                    '--params', str(ed_params), '--output', str(target / 'results.json'),
                    '--dense-memory-cap-gib', str(args.dense_memory_cap_gib)]

            else:
                target = output / case['id']
                command = []
            if args.mode != 'analyze' and target.exists() and any(target.iterdir()):
                raise ValueError(f'nonempty run directory: {target}; choose a fresh --output')
            prepared.append((case, inputs, ed_params, target, command))
            if args.dry_run:
                print(f'{case["id"]}: cwd={target} ' + (shlex.join(command) if command else ('exact_free_boson' if args.mode == 'ed' else 'analyze')))
        if args.dry_run:
            return 0
        if args.mode == 'dqmc':
            subprocess.run(['make', '-C', str(ROOT), 'build'], env=env, check=True)
        records = []
        unavailable = False
        for case, inputs, ed_params, target, command in prepared:
            if args.mode in {'dqmc', 'ed'}:
                empty_output(target)
                if args.mode == 'dqmc':
                    for name in INPUTS:
                        shutil.copy2(inputs / name, target / name)
                else:
                    shutil.copy2(ed_params, target / 'params.json')
                print(f'Running {args.mode}: {case["id"]}', flush=True)
                (target / 'run_metadata.json').write_text(json.dumps({
                    'case': case, 'command': command, 'manifest': str(manifest_path),
                    'threads': {k: env[k] for k in ('OMP_NUM_THREADS', 'MKL_NUM_THREADS', 'OPENBLAS_NUM_THREADS')},
                }, indent=2) + '\n')
                with (target / 'run.log').open('w') as log:
                    if args.mode == 'ed' and not command:
                        # The manuscript's free point uses the exact grand-canonical
                        # reference, so no artificial many-body cutoff is introduced.
                        payload = json.loads(ed_params.read_text())
                        reference = free_boson_reference_result(
                            payload['parameters'], payload.get('convergence_policy'))
                        (target / 'results.json').write_text(json.dumps(reference, indent=2, allow_nan=False) + '\n')
                        log.write('Computed fresh exact grand-canonical free-boson reference.\n')
                    else:
                        subprocess.run(command, cwd=target, env=env, stdout=log, stderr=subprocess.STDOUT, check=True)
                validate_completed_output(args.mode, target)
            else:
                settings = {**manifest.get('analysis', {}), **case.get('analysis', {})}
                params_payload = json.loads(ed_params.read_text()) if ed_params.is_file() else {}
                parameters = case.get('parameters', params_payload.get('parameters'))
                if parameters is None:
                    raise ValueError(f'{case["id"]} must define physical parameters')
                record = compare_dqmc_ed_case(
                    target / 'dqmc', block_size=int(settings.get('block_size', 1000)),
                    lq=int(parameters['Lx']) * int(parameters['Ly']),
                    expected_parameters=parameters,
                    expected_policy=params_payload.get('convergence_policy'),
                    ed_result_path=target / 'ed' / 'results.json',
                    skip_samples=int(settings.get('skip_samples', 0)),
                    stderr_tolerance=float(settings.get('stderr_tolerance', 3.0)),
                    atol=float(settings.get('atol', 0.0)), rtol=float(settings.get('rtol', 0.0)),
                )
                record.update(id=case['id'], parameters=parameters)
                records.append(record)
                unavailable |= record['status'] in {'dqmc_unavailable', 'ed_unavailable'}
        if args.mode == 'analyze':
            output.mkdir(parents=True, exist_ok=True)
            path = output / 'analysis.json'
            path.write_text(json.dumps(json_clean({'manifest': str(manifest_path), 'records': records}), indent=2, allow_nan=False) + '\n')
            print(path)
            return 1 if unavailable else 0
        return 0
    except (KeyError, ValueError, OSError, subprocess.CalledProcessError) as exc:
        print(f'ERROR: {exc}', file=sys.stderr)
        return 2


if __name__ == '__main__':
    raise SystemExit(main())
