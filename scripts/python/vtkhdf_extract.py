"""Inspect Neko's temporal VTKHDF output and extract single steps from it.

Neko and Neko-TOP write temporal VTKHDF output as a small main file,
<stem>.vtkhdf, holding the mesh and the step metadata, plus one HDF5 file
per output step, <stem>.data/<step>.h5, joined through virtual datasets.
This tool works on such output with bounded memory and never modifies it.

Usage:
    python3 vtkhdf_extract.py info FILE
    python3 vtkhdf_extract.py extract FILE (--step S | --time T)
                                      [-o OUT] [--overwrite]

info
    Print the type and version, the steps and their time range, the mesh
    sizes and cell types, every field with its dtype, components and
    storage (printf VDS, explicit VDS or inline), and the size of the
    main file and the per-step files, listing any that are missing. The
    per-step files are only stat'ed, never opened.

extract
    Write one step as a standalone, non-temporal VTKHDF file: the mesh, the
    step's point data with all data materialised (no virtual datasets) and
    the step's time as FieldData/Time. --step accepts negative values (-1
    is the last step); --time picks the nearest step and refuses a tie.
    The default output is <input_dir>/<stem>_step_<s>.vtkhdf, with <s>
    zero-padded to the width of the last step index. An existing output is
    only replaced with --overwrite, and never if it is the input or one of
    its per-step files.

Errors are reported as a single line on stderr with exit status 1.
Requires h5py and numpy (see requirements.txt).
"""

import sys
from argparse import ArgumentParser, Namespace
from pathlib import Path

import numpy as np

from vtkhdf_backend import (FileInfo, VtkhdfError, cell_types, find_step,
                            format_size, inspect, open_input, stat_sources,
                            step_label, write_extract)

VTK_CELL_NAMES = {
    9: "VTK_QUAD",
    12: "VTK_HEXAHEDRON",
    70: "VTK_LAGRANGE_QUADRILATERAL",
    72: "VTK_LAGRANGE_HEXAHEDRON",
}
STORAGE_NAMES = {
    "printf": "printf VDS",
    "explicit": "explicit VDS",
    "inline": "inline",
}
MAX_MISSING_LISTED = 10


def parse_args() -> Namespace:
    parser = ArgumentParser(
        description="Inspect Neko's temporal VTKHDF output and extract "
        "single steps from it.")
    commands = parser.add_subparsers(dest="command", required=True)

    info = commands.add_parser(
        "info", help="describe a VTKHDF file and its per-step files")
    info.add_argument("file", type=Path, help="main .vtkhdf file")

    extract = commands.add_parser(
        "extract", help="write one step as a static VTKHDF file")
    extract.add_argument("file", type=Path, help="main .vtkhdf file")
    which = extract.add_mutually_exclusive_group(required=True)
    which.add_argument("--step", type=int,
                       help="step index; negative counts from the end")
    which.add_argument("--time", type=float,
                       help="pick the step nearest to this time")
    extract.add_argument(
        "-o", "--output", type=Path, default=None,
        help="output file (default <input_dir>/<stem>_step_<s>.vtkhdf)")
    extract.add_argument("--overwrite", action="store_true",
                         help="replace an existing output file")
    return parser.parse_args()


def print_warnings(info: FileInfo) -> None:
    for message in info.warnings:
        print(f"warning: {message}", file=sys.stderr)


def describe_rows(row_shape: tuple[int, ...]) -> str:
    if not row_shape:
        return "scalar"
    if len(row_shape) == 1:
        return f"{row_shape[0]} components"
    return "x".join(str(n) for n in row_shape) + " components"


def describe_time(info: FileInfo) -> str:
    if not info.temporal:
        return "no (static)"
    times = info.times
    order = ("strictly increasing" if np.all(np.diff(times) > 0)
             else "NOT strictly increasing")
    return (f"yes: {info.n_steps} steps, t = {times[0]:g} .. "
            f"{times[-1]:g}, {order}")


def print_sources(info: FileInfo) -> None:
    sources = stat_sources(info)
    sizes = [nbytes for _, nbytes, missing in sources if not missing]
    missing = [(step, path) for step, _, paths in sources for path in paths]
    if not any(fld.kind != "inline" for fld in info.fields):
        print("step files  none (all point data inline)")
        return
    complete = len(sizes)
    line = (f"step files  {complete} of {info.n_steps} steps complete, "
            f"{format_size(sum(n for _, n, _ in sources))} total")
    if sizes:
        low, high = min(sizes), max(sizes)
        per_step = (format_size(low) if low == high
                    else f"{format_size(low)} .. {format_size(high)}")
        line += f", {per_step} per step"
    print(line)
    if not missing:
        print("missing     none")
        return
    print(f"missing     {len(missing)} file(s):")
    for step, path in missing[:MAX_MISSING_LISTED]:
        print(f"  step {step}: {path}")
    if len(missing) > MAX_MISSING_LISTED:
        print(f"  ... and {len(missing) - MAX_MISSING_LISTED} more")


def run_info(args: Namespace) -> None:
    with open_input(args.file) as f:
        info = inspect(f, args.file)
        types = cell_types(f)
    print_warnings(info)

    version = ".".join(str(v) for v in info.version) or "?"
    type_names = ", ".join(
        f"{t} ({VTK_CELL_NAMES.get(t, 'unknown')})" for t in types)
    print(f"file        {info.path} "
          f"({format_size(info.path.stat().st_size)})")
    print(f"type        {info.vtk_type}, VTKHDF version {version}")
    print(f"temporal    {describe_time(info)}")
    print(f"mesh        {info.n_points} points, {info.n_cells} cells, "
          f"{info.n_conn} connectivity ids, {info.n_parts} part(s)")
    print(f"cell types  {type_names or 'none'}")
    print(f"fields      {len(info.fields)}")
    width = max((len(fld.name) for fld in info.fields), default=0)
    for fld in info.fields:
        storage = STORAGE_NAMES[fld.kind]
        if fld.kind == "printf":
            storage += f" {fld.names[0]}"
        elif fld.kind == "explicit":
            storage += f" ({len(fld.names)} mappings)"
        print(f"  {fld.name:<{width}}  {str(fld.dtype):<8} "
              f"{describe_rows(fld.row_shape):<13} {storage}")
    print_sources(info)


def run_extract(args: Namespace) -> None:
    with open_input(args.file) as f:
        info = inspect(f, args.file)
        print_warnings(info)
        step = find_step(info, args.step, args.time)
        out = args.output
        if out is None:
            label = step_label(info.n_steps, step)
            out = info.path.parent / f"{info.path.stem}_step_{label}.vtkhdf"
        write_extract(f, info, step, out, args.overwrite)
    print(f"extracted step {step} (t = {info.times[step]:g}) of "
          f"{info.path} -> {out} ({format_size(out.stat().st_size)})")


def main() -> None:
    args = parse_args()
    commands = {"info": run_info, "extract": run_extract}
    try:
        commands[args.command](args)
    except (VtkhdfError, OSError) as error:
        print(f"error: {error}", file=sys.stderr)
        sys.exit(1)


if __name__ == "__main__":
    main()
