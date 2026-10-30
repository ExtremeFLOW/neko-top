"""Inspect, extract steps from and split Neko's temporal VTKHDF output.

Neko and Neko-TOP write temporal VTKHDF output as a small main file,
<stem>.vtkhdf, holding the mesh and the step metadata, plus one HDF5 file
per output step, <stem>.data/<step>.h5, joined through virtual datasets.
This tool works on such output with bounded memory and never modifies it.

Usage:
    python3 vtkhdf_extract.py info FILE
    python3 vtkhdf_extract.py extract FILE (--step S | --time T)
                                      [-o OUT] [--overwrite]
    python3 vtkhdf_extract.py split FILE (--steps-per-chunk N
                                    | --time-window W | --max-size SIZE)
                                    [--dry-run] [--overwrite]

info
    Print the type and version, the steps and their time range, the mesh
    sizes and cell types, and every field with its dtype, components and
    storage (printf VDS, explicit VDS or inline). Then the sizes: main
    file, per-step files and total, the size of one step's files (min,
    mean and max if they differ), any missing or unreadable per-step
    files (excluded from the sizes), and the footprint model split
    --max-size uses, with the smallest --max-size it accepts. The
    per-step files are only stat'ed, never opened.

extract
    Write one step as a standalone, non-temporal VTKHDF file: the mesh, the
    step's point data with all data materialised (no virtual datasets) and
    the step's time as FieldData/Time. --step accepts negative values (-1
    is the last step); --time picks the nearest step and refuses a tie.
    The default output is <input_dir>/<stem>_step_<s>.vtkhdf, with <s>
    zero-padded to the width of the last step index.

split
    Write chunks <stem>_steps_<a>-<b>.vtkhdf next to the input, each a
    temporal VTKHDF file for steps a..b (original times kept): a copy of
    the mesh plus point data as virtual datasets referencing the original
    per-step files, so no field data is copied. A manifest
    <stem>_steps_<a>-<b>.files lists the chunk and its per-step files
    relative to the input's directory, so that
        rsync -a --files-from=<manifest> <input_dir>/ HOST:DEST/
    ships exactly one self-contained chunk. Chunks hold N consecutive steps
    (--steps-per-chunk), the steps of each absolute time window
    [m W, (m+1) W) (--time-window, needs strictly increasing times), or as
    many steps as fit within SIZE bytes of footprint, chunk file plus its
    per-step files (--max-size, as estimated by the model info prints).
    SIZE is bytes, optionally with K, M, G, T or KiB, MiB, GiB, TiB (all
    powers of 1024, any case, e.g. 40G or "233.5 MiB"); every size the
    tool prints is accepted. KB, MB, GB, TB are refused as ambiguous.
    A 1-step chunk also gets FieldData/Time, as VTK reads it as static.
    a and b are step indices of FILE, so a chunk can itself be split or
    extracted from. Every per-step file and every output name is validated
    before anything is written; --dry-run does the same validation and
    prints the plan but writes nothing.

An existing output is only replaced with --overwrite, and never if it is
the input or one of its per-step files. Each output is written to a
temporary file and moved into place, so an error never leaves a partial
file behind.

Errors are reported as a single line on stderr with exit status 1;
SIGTERM, SIGHUP and Ctrl-C remove any partial output and exit with
128 + the signal number. Set VTKHDF_EXTRACT_TRACEBACK=1 to see the
traceback of an unexpected (internal) error.
Requires h5py and numpy (see requirements.txt).
"""

import os
import shlex
import signal
import sys
from argparse import ArgumentParser, Namespace
from pathlib import Path

import numpy as np

from vtkhdf_backend import (CHUNK_OVERHEAD, MAPPING_OVERHEAD, Chunk,
                            FileInfo, VtkhdfError, cell_types,
                            chunk_estimate, find_step, format_size, inspect,
                            mesh_bytes, open_input, parse_size,
                            plan_by_count, prepare_split, protected_paths,
                            size_model, smallest_max_size, stat_sources,
                            step_files, step_label, strictly_increasing,
                            write_chunk, write_extract, write_manifest)

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
EXAMPLE_STEPS = (10, 100, 1000, 10000)


class Terminated(BaseException):
    """Raised by the SIGTERM/SIGHUP handler, so that partial outputs are
    removed on the way out."""

    def __init__(self, signum: int):
        super().__init__(signum)
        self.signum = signum


def parse_args() -> Namespace:
    parser = ArgumentParser(
        description="Inspect, extract steps from and split Neko's temporal "
        "VTKHDF output.")
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

    split = commands.add_parser(
        "split", help="write chunks of steps that reference the originals")
    split.add_argument("file", type=Path, help="main .vtkhdf file")
    mode = split.add_mutually_exclusive_group(required=True)
    mode.add_argument("--steps-per-chunk", type=int, metavar="N",
                      help="N consecutive steps per chunk")
    mode.add_argument("--time-window", type=float, metavar="W",
                      help="one chunk per time window [m W, (m+1) W)")
    mode.add_argument(
        "--max-size", metavar="SIZE",
        help="largest footprint per chunk: bytes, or with K, M, G, T or "
        "KiB, MiB, GiB, TiB (powers of 1024), e.g. 40G or '233.5 MiB'")
    split.add_argument("--dry-run", action="store_true",
                       help="validate and print the plan, write nothing")
    split.add_argument("--overwrite", action="store_true",
                       help="replace existing chunks and manifests")
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
    finite = times[np.isfinite(times)]
    span = (f"t = {finite[0]:g} .. {finite[-1]:g}" if len(finite)
            else "no finite time")
    order = ("strictly increasing" if strictly_increasing(times)
             else "NOT strictly increasing")
    if len(finite) < len(times):
        order += f", {len(times) - len(finite)} non-finite"
    return f"yes: {info.n_steps} steps, {span}, {order}"


def sized(nbytes: int) -> str:
    return f"{format_size(nbytes)} ({nbytes} B)"


def print_sizes(f, info: FileInfo) -> None:
    """Sizes from stat only, and the split footprint model (the same
    size_model split --max-size plans with)."""
    main = info.path.stat().st_size
    print(f"main file   {sized(main)}")
    if not any(fld.kind != "inline" for fld in info.fields):
        print("step files  none (all point data inline)")
        print(f"total       {sized(main)}")
        print("split size  not possible: no per-step files to reference")
        return
    sizes, problems = stat_sources(info)
    bad = {path for _, path, _ in problems}
    complete = [s for s in range(info.n_steps)
                if not bad.intersection(step_files(info, s))]
    print(f"step files  {sized(sum(sizes.values()))} in {len(sizes)} "
          f"files, {len(complete)} of {info.n_steps} steps complete")
    print(f"total       {sized(main + sum(sizes.values()))}")
    per = [sum(sizes[p] for p in step_files(info, s)) for s in complete]
    if per and min(per) == max(per):
        print(f"per step    {sized(per[0])}")
    elif per:
        print(f"per step    min {sized(min(per))}, mean "
              f"{format_size(round(sum(per) / len(per)))},")
        print(f"            max {sized(max(per))}")
    if problems:
        print(f"missing     {len(problems)} file(s), excluded from the "
              "sizes above:")
        for step, path, why in problems[:MAX_MISSING_LISTED]:
            print(f"  step {step}: {path} ({why})")
        if len(problems) > MAX_MISSING_LISTED:
            print(f"  ... and {len(problems) - MAX_MISSING_LISTED} more")
    else:
        print("missing     none")
    if not info.temporal:
        return
    fixed, per_step = size_model(f, info, sizes)
    smallest = smallest_max_size(fixed, per_step)
    print(f"split size  per chunk: mesh {sized(mesh_bytes(f))}")
    print(f"            + overhead {sized(CHUNK_OVERHEAD)}")
    print(f"            per step: its files + {MAPPING_OVERHEAD} B per "
          f"field ({MAPPING_OVERHEAD * len(info.fields)} B)")
    print(f"            smallest --max-size "
          f"{format_size(smallest, round_up=True)} (exactly {smallest} B)")
    n = info.n_steps
    for k in [k for k in EXAMPLE_STEPS if k < n] + [n]:
        worst = max(chunk_estimate(fixed, per_step, steps)
                    for steps in plan_by_count(n, k))
        label = (f"all {n} steps in one chunk" if k == n
                 else f"{k} steps per chunk")
        print(f"            {label}: about {format_size(worst)}")


def run_info(args: Namespace) -> None:
    with open_input(args.file) as f:
        info = inspect(f, args.file)
        types = cell_types(f)
        print_warnings(info)
        version = ".".join(str(v) for v in info.version) or "?"
        type_names = ", ".join(
            f"{t} ({VTK_CELL_NAMES.get(t, 'unknown')})" for t in types)
        print(f"file        {info.path}")
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
        print_sizes(f, info)


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


def describe_chunk(info: FileInfo, chunk: Chunk, footprint: str) -> None:
    steps, times = chunk.steps, info.times
    print(f"{chunk.path.name}  steps {steps[0]}-{steps[-1]} "
          f"({len(steps)})  t = {times[steps[0]]:g} .. "
          f"{times[steps[-1]]:g}  footprint {footprint}")
    source_dir = shlex.quote(f"{info.path.parent}/")
    print(f"  rsync -a --files-from={shlex.quote(str(chunk.manifest))} "
          f"{source_dir} HOST:DEST/")


def run_split(args: Namespace) -> None:
    max_size = None if args.max_size is None else parse_size(args.max_size)
    with open_input(args.file) as f:
        info = inspect(f, args.file)
        print_warnings(info)
        chunks, warnings = prepare_split(
            f, info, steps_per_chunk=args.steps_per_chunk,
            time_window=args.time_window, max_size=max_size,
            overwrite=args.overwrite)
        for message in warnings:
            print(f"warning: {message}", file=sys.stderr)
        print(f"split {info.path}: {info.n_steps} steps into "
              f"{len(chunks)} chunk(s)")
        if args.dry_run:
            for chunk in chunks:
                describe_chunk(info, chunk,
                               f"~{format_size(chunk.estimate)} (estimate)")
            print("dry run: nothing written")
            return
        protected = protected_paths(info)
        for chunk in chunks:
            write_chunk(f, info, chunk.steps, chunk.path, args.overwrite,
                        protected)
            write_manifest(chunk.manifest, [chunk.path.name, *chunk.files],
                           args.overwrite)
            footprint = chunk.path.stat().st_size + chunk.source_bytes
            describe_chunk(info, chunk, format_size(footprint))
            if max_size is not None and footprint > max_size:
                print(f"warning: {chunk.path.name} footprint "
                      f"{format_size(footprint)} exceeds --max-size "
                      f"{format_size(max_size)}", file=sys.stderr)


def on_signal(signum: int, frame) -> None:
    raise Terminated(signum)


def main() -> None:
    args = parse_args()
    # Paths are printed as the file system's bytes, whatever the locale.
    sys.stdout.reconfigure(errors="surrogateescape")
    for name in ("SIGTERM", "SIGHUP"):
        if hasattr(signal, name):
            signal.signal(getattr(signal, name), on_signal)
    commands = {"info": run_info, "extract": run_extract,
                "split": run_split}
    try:
        commands[args.command](args)
    except (VtkhdfError, OSError) as error:
        print(f"error: {error}", file=sys.stderr)
        sys.exit(1)
    except Terminated as stop:
        print(f"error: terminated by {signal.Signals(stop.signum).name}",
              file=sys.stderr)
        sys.exit(128 + stop.signum)
    except KeyboardInterrupt:
        print("error: interrupted", file=sys.stderr)
        sys.exit(130)
    except Exception as error:  # last resort: no traceback for users
        if os.environ.get("VTKHDF_EXTRACT_TRACEBACK"):
            raise
        print(f"error: internal error: {type(error).__name__}: {error}",
              file=sys.stderr)
        sys.exit(1)


if __name__ == "__main__":
    main()
