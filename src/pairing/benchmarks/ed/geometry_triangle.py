#!/usr/bin/env python3
from __future__ import annotations

import argparse
import json
import math

A1 = (1.0, 0.0)
A2 = (0.5, math.sqrt(3.0) / 2.0)
B1 = (2.0 * math.pi, -2.0 * math.pi / math.sqrt(3.0))
B2 = (0.0, 4.0 * math.pi / math.sqrt(3.0))


def site_index(x: int, y: int, lx: int, ly: int) -> int:
    return (y % ly) * lx + (x % lx)


def site_coord(index: int, lx: int) -> tuple[int, int]:
    return index % lx, index // lx


def triangular_bonds(lx: int, ly: int) -> list[tuple[int, int]]:
    bonds: list[tuple[int, int]] = []
    for y in range(ly):
        for x in range(lx):
            i = site_index(x, y, lx, ly)
            bonds.append((i, site_index(x + 1, y, lx, ly)))
            bonds.append((i, site_index(x, y + 1, lx, ly)))
            bonds.append((i, site_index(x - 1, y + 1, lx, ly)))
    return bonds


def bond_table(lx: int, ly: int) -> dict[str, list[int]]:
    table: dict[str, list[int]] = {}
    bonds = triangular_bonds(lx, ly)
    for i in range(lx * ly):
        table[str(i)] = [j for start, j in bonds if start == i]
    return table


def k_points(lx: int, ly: int) -> list[dict[str, float | int]]:
    points: list[dict[str, float | int]] = []
    for iy in range(1, ly + 1):
        for ix in range(1, lx + 1):
            kx = (ix - 1) * B1[0] / lx + (iy - 1) * B2[0] / ly
            ky = (ix - 1) * B1[1] / lx + (iy - 1) * B2[1] / ly
            points.append({"index": len(points), "ix": ix, "iy": iy, "kx": kx, "ky": ky})
    return points


def payload(lx: int, ly: int) -> dict[str, object]:
    return {
        "Lx": lx,
        "Ly": ly,
        "a1": A1,
        "a2": A2,
        "b1": B1,
        "b2": B2,
        "bonds": bond_table(lx, ly),
        "k_points": k_points(lx, ly),
        "k_zero_index": 0,
    }


def assert_defaults() -> None:
    assert bond_table(2, 2) == {
        "0": [1, 2, 3],
        "1": [0, 3, 2],
        "2": [3, 0, 1],
        "3": [2, 1, 0],
    }
    assert bond_table(2, 1) == {
        "0": [1, 0, 1],
        "1": [0, 1, 0],
    }
    assert k_points(2, 2)[0]["kx"] == 0.0
    assert k_points(2, 2)[0]["ky"] == 0.0


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--Lx", type=int, default=2)
    parser.add_argument("--Ly", type=int, default=2)
    parser.add_argument("--assert-defaults", action="store_true")
    args = parser.parse_args()
    if args.assert_defaults:
        assert_defaults()
    print(json.dumps(payload(args.Lx, args.Ly), indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
