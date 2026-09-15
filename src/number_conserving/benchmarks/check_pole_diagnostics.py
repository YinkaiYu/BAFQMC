#!/usr/bin/env python3
"""Validate pole diagnostic output files from a DQMC run directory."""

from __future__ import annotations

import argparse
import math
import re
import sys
from pathlib import Path


SCALAR_FILES = (
    "pole_distance",
    "pole_x",
    "green_spectral_radius",
    "green_smax",
    "log_weight",
)


class DiagnosticError(Exception):
    """Raised when diagnostic files fail consistency checks."""


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Check DQMC pole diagnostic files in a run directory."
    )
    parser.add_argument("run_dir", type=Path)
    parser.add_argument("--ranks", type=int)
    parser.add_argument("--rtol", type=float, default=1e-9)
    parser.add_argument("--atol", type=float, default=1e-12)
    return parser.parse_args()


def first_numeric_rows(path: Path, count: int) -> list[list[float]]:
    rows: list[list[float]] = []
    for line_number, line in enumerate(path.read_text(encoding="utf-8").splitlines(), 1):
        stripped = line.strip()
        if not stripped:
            continue
        if len(rows) == count:
            break
        fields = stripped.split()
        try:
            rows.append([float(field) for field in fields])
        except ValueError as exc:
            row_number = len(rows) + 1
            raise DiagnosticError(
                f"{path.name} required row {row_number} (line {line_number}) is not numeric"
            ) from exc
    return rows


def parse_param_file(run_dir: Path) -> tuple[int, int, int]:
    path = run_dir / "paramC_sets.txt"
    if not path.is_file():
        raise DiagnosticError("missing paramC_sets.txt")

    rows = first_numeric_rows(path, 4)
    if len(rows) < 4:
        raise DiagnosticError("paramC_sets.txt must contain at least four nonblank rows")
    if len(rows[1]) < 2:
        raise DiagnosticError("paramC_sets.txt lattice row must contain Lx and Ly")
    if len(rows[3]) < 2:
        raise DiagnosticError("paramC_sets.txt update row must contain Nbin")

    lx = int(rows[1][0])
    ly = int(rows[1][1])
    nbin = int(rows[3][1])
    if lx <= 0 or ly <= 0 or nbin <= 0:
        raise DiagnosticError("Lx, Ly, and Nbin must be positive")
    return lx, ly, nbin


def parse_info_ranks(run_dir: Path) -> int | None:
    path = run_dir / "info.txt"
    if not path.is_file():
        return None

    pattern = re.compile(r"^\s*#\s*Cores\s*:\s*(\d+)\s*$")
    for line in path.read_text(encoding="utf-8").splitlines():
        match = pattern.match(line)
        if match:
            ranks = int(match.group(1))
            if ranks <= 0:
                raise DiagnosticError("info.txt # Cores must be positive")
            return ranks
    return None


def parse_ranks(run_dir: Path, cli_ranks: int | None) -> int:
    if cli_ranks is not None:
        if cli_ranks <= 0:
            raise DiagnosticError("--ranks must be positive")
        return cli_ranks
    return parse_info_ranks(run_dir) or 1


def read_scalar_file(path: Path) -> list[float]:
    if not path.is_file():
        raise DiagnosticError(f"missing {path.name}")

    values: list[float] = []
    for line_number, line in enumerate(path.read_text(encoding="utf-8").splitlines(), 1):
        stripped = line.strip()
        if not stripped:
            continue
        fields = stripped.split()
        if len(fields) != 1:
            raise DiagnosticError(f"{path.name}:{line_number} must contain one float")
        try:
            value = float(fields[0])
        except ValueError as exc:
            raise DiagnosticError(f"{path.name}:{line_number} is not a float") from exc
        if not math.isfinite(value):
            raise DiagnosticError(f"{path.name}:{line_number} is not finite")
        values.append(value)
    return values


def read_pole_z(path: Path, ndim: int) -> list[list[float]]:
    if not path.is_file():
        raise DiagnosticError("missing pole_z")

    expected_columns = 2 * ndim
    rows: list[list[float]] = []
    for line_number, line in enumerate(path.read_text(encoding="utf-8").splitlines(), 1):
        stripped = line.strip()
        if not stripped:
            continue
        fields = stripped.split()
        if len(fields) != expected_columns:
            raise DiagnosticError(
                f"pole_z:{line_number} has {len(fields)} columns; "
                f"expected {expected_columns}"
            )
        try:
            row = [float(field) for field in fields]
        except ValueError as exc:
            raise DiagnosticError(f"pole_z:{line_number} contains a non-float") from exc
        if not all(math.isfinite(value) for value in row):
            raise DiagnosticError(f"pole_z:{line_number} contains a non-finite value")
        rows.append(row)
    return rows


def require_close(name: str, sample: int, actual: float, expected: float, rtol: float, atol: float) -> None:
    if not math.isclose(actual, expected, rel_tol=rtol, abs_tol=atol):
        raise DiagnosticError(
            f"{name}[{sample}]={actual:.17g} does not match expected {expected:.17g}"
        )


def check_diagnostics(run_dir: Path, ranks_arg: int | None, rtol: float, atol: float) -> None:
    if not run_dir.is_dir():
        raise DiagnosticError(f"{run_dir} is not a directory")
    if rtol < 0 or atol < 0:
        raise DiagnosticError("--rtol and --atol must be non-negative")

    lx, ly, nbin = parse_param_file(run_dir)
    ranks = parse_ranks(run_dir, ranks_arg)
    ndim = lx * ly

    pole_z = read_pole_z(run_dir / "pole_z", ndim)
    scalar_values = {name: read_scalar_file(run_dir / name) for name in SCALAR_FILES}
    counts = {"pole_z": len(pole_z), **{name: len(values) for name, values in scalar_values.items()}}
    sample_count = counts["pole_z"]
    if any(count != sample_count for count in counts.values()):
        detail = ", ".join(f"{name}={count}" for name, count in sorted(counts.items()))
        raise DiagnosticError(f"diagnostic files have inconsistent sample counts: {detail}")
    if sample_count <= 0:
        raise DiagnosticError("diagnostic files contain no samples")

    expected_chunk = nbin * ranks
    if sample_count % expected_chunk != 0:
        raise DiagnosticError(
            f"sample count {sample_count} is not a positive multiple of "
            f"Nbin*ranks ({nbin}*{ranks}={expected_chunk})"
        )

    distances = scalar_values["pole_distance"]
    pole_x = scalar_values["pole_x"]
    radii = scalar_values["green_spectral_radius"]
    smax_values = scalar_values["green_smax"]
    for idx, distance in enumerate(distances, 1):
        if distance <= 0.0:
            raise DiagnosticError(f"pole_distance[{idx}] must be positive")
        require_close("pole_x", idx, pole_x[idx - 1], -math.log10(distance), rtol, atol)
        require_close(
            "green_spectral_radius",
            idx,
            radii[idx - 1],
            1.0 / distance,
            rtol,
            atol,
        )
        allowed_gap = max(atol, rtol * abs(radii[idx - 1]))
        if smax_values[idx - 1] + allowed_gap < radii[idx - 1]:
            raise DiagnosticError(
                f"green_smax[{idx}]={smax_values[idx - 1]:.17g} is below "
                f"green_spectral_radius[{idx}]={radii[idx - 1]:.17g}"
            )


def main() -> int:
    args = parse_args()
    try:
        check_diagnostics(args.run_dir, args.ranks, args.rtol, args.atol)
    except DiagnosticError as exc:
        print(f"ERROR: {exc}")
        return 1

    print(f"Pole diagnostic files passed: {args.run_dir}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
