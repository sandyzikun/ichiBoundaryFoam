#!/usr/bin/env python3
"""Compare paired OpenFOAM ASCII fields at their latest times."""

from __future__ import annotations

import math
import pathlib
import re
import sys


def latest_time(case: pathlib.Path) -> str:
    values = []
    for path in case.iterdir():
        try:
            values.append((float(path.name), path.name))
        except (ValueError, TypeError):
            pass
    return max(values)[1]


def internal_values(path: pathlib.Path) -> list[tuple[float, ...]]:
    text = path.read_text()
    uniform = re.search(r"internalField\s+uniform\s+([^;]+);", text)
    if uniform:
        token = uniform.group(1).strip().strip("()")
        return [tuple(float(x) for x in token.split())]
    match = re.search(
        r"internalField\s+nonuniform\s+List<[^>]+>\s+"
        r"(\d+)\s*\((.*?)\)\s*;",
        text,
        re.S,
    )
    if not match:
        raise ValueError(f"Cannot parse {path}")
    body = match.group(2)
    tokens = re.findall(r"\([^()]+\)|[-+0-9.eE]+", body)
    return [
        tuple(float(x) for x in token.strip("()").split())
        for token in tokens
    ]


def metrics(a, b):
    if len(a) != len(b):
        raise ValueError(f"Field sizes differ: {len(a)} != {len(b)}")
    diffs = []
    refs = []
    for va, vb in zip(a, b):
        if len(va) != len(vb):
            raise ValueError("Field component counts differ")
        diffs.extend(x - y for x, y in zip(va, vb))
        refs.extend(x for x in va)
    maximum = max(abs(x) for x in diffs)
    rms = math.sqrt(sum(x * x for x in diffs) / len(diffs))
    ref_rms = math.sqrt(sum(x * x for x in refs) / len(refs))
    relative = rms / ref_rms if ref_rms else float("nan")
    return maximum, rms, relative


def main(root: pathlib.Path) -> None:
    cases = root / "cases"
    print("wave,model,field,time,maxAbs,rms,relativeRms")
    for wave in ("stokesI", "stokesII", "stokesV"):
        for model in ("laminar", "kEpsilon", "kOmegaSST"):
            lhs = cases / f"{wave}-{model}-interFoam"
            rhs = cases / f"{wave}-{model}-ichiBoundaryFoam"
            if not (lhs / "RUN_COMPLETE").exists() or not (
                rhs / "RUN_COMPLETE"
            ).exists():
                print(
                    f"{wave},{model},SKIPPED,incomplete,,,,"
                )
                continue
            ta = latest_time(lhs)
            tb = latest_time(rhs)
            if float(ta) != float(tb):
                raise ValueError(f"Latest times differ: {ta} != {tb}")
            for field in ("alpha.water", "U", "p_rgh"):
                result = metrics(
                    internal_values(lhs / ta / field),
                    internal_values(rhs / tb / field),
                )
                print(
                    f"{wave},{model},{field},{ta},"
                    f"{result[0]:.12g},{result[1]:.12g},"
                    f"{result[2]:.12g}"
                )


if __name__ == "__main__":
    main(pathlib.Path(sys.argv[1] if len(sys.argv) > 1 else "."))
