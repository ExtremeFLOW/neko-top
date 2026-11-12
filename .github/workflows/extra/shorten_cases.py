#!/usr/bin/env python3
"""Shorten Neko-TOP case files so the examples can run as CI smoke tests.

Each case keeps its own timestep, which is chosen for the stability of its
mesh, but its end time is capped to a few timesteps and its optimization loop
to a few iterations. Values that are already smaller are left unchanged, and
linear solver settings are not touched.

Neko reads case files with json-fortran, which accepts a single trailing comma
in objects and arrays. Those commas are removed before parsing, and the files
are written back as plain JSON.
"""

import argparse
import json
import pathlib
import re
import sys

TRAILING_COMMA = re.compile(r",(\s*[}\]])")


def shorten(path, steps, iterations):
    """Shorten one case file in place and return a description of changes."""
    data = json.loads(TRAILING_COMMA.sub(r"\1", path.read_text()))
    changes = []

    time = data.get("case", {}).get("time", {})
    if "end_time" in time and "timestep" in time:
        end_time = float(f"{steps * time['timestep']:.6g}")
        if end_time < time["end_time"]:
            changes.append(f"end_time {time['end_time']} -> {end_time}")
            time["end_time"] = end_time

    solver = data.get("optimization", {}).get("solver", {})
    if solver.get("max_iterations", 0) > iterations:
        changes.append(
            f"optimization iterations {solver['max_iterations']} -> {iterations}"
        )
        solver["max_iterations"] = iterations

    path.write_text(json.dumps(data, indent=4) + "\n")
    return changes


def main():
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument(
        "paths", nargs="+", type=pathlib.Path,
        help="case files, or directories to search for *.case files")
    parser.add_argument(
        "--steps", type=int, default=10,
        help="number of timesteps to run (default: %(default)s)")
    parser.add_argument(
        "--iterations", type=int, default=2,
        help="number of optimization iterations (default: %(default)s)")
    args = parser.parse_args()

    cases = []
    for path in args.paths:
        cases.extend(sorted(path.rglob("*.case")) if path.is_dir() else [path])
    if not cases:
        print("No case files found.", file=sys.stderr)
        return 1

    failed = False
    for case in cases:
        try:
            changes = shorten(case, args.steps, args.iterations)
        except (ValueError, TypeError) as error:
            print(f"{case}: could not shorten: {error}", file=sys.stderr)
            failed = True
            continue
        print(f"{case}: {', '.join(changes) or 'unchanged'}")

    return 1 if failed else 0


if __name__ == "__main__":
    sys.exit(main())
