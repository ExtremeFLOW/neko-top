"""HDF5 backend for vtkhdf_extract.py.

Neko writes temporal VTKHDF output as a main file ``<stem>.vtkhdf`` (mesh
and step metadata) whose PointData are virtual datasets (VDS) over one HDF5
file per output step, ``<stem>.data/<step>.h5``. This module reads such
files with bounded memory and writes derived files without ever modifying
the originals:

- inputs and per-step source files are only opened through open_input(),
  read-only, one source file at a time;
- an input PointData VDS is never sized, sliced or read: its shape, len(),
  ndim and repr() make HDF5 open every per-step file. Only its dtype,
  is_virtual and creation property list are used, and step data is read
  from the source file directly;
- bulk data only moves through copy_rows(), whose buffer is at most
  BLOCK_BYTES; metadata reads are O(NSteps + number of parts);
- a split chunk copies only the mesh: its point data are virtual datasets
  referencing the original per-step files, so it moves no field data.
"""

import contextlib
import math
import os
import re
from collections.abc import Iterator
from dataclasses import dataclass
from pathlib import Path

import h5py
import numpy as np
from h5py import h5a, h5d, h5p, h5s

BLOCK_BYTES = 64 * 1024**2
"""Largest buffer used for any bulk copy."""

TIME_EPS = 1e-9
"""Relative tolerance used to compare times."""

MESH_DATASETS = ("Points", "Connectivity", "Offsets", "Types")
COUNT_DATASETS = ("NumberOfPoints", "NumberOfCells",
                  "NumberOfConnectivityIds")
STEP_OFFSETS = ("PartOffsets", "PointOffsets", "CellOffsets",
                "ConnectivityIdOffsets")
SIZE_UNITS = "KMGT"
CHUNK_OVERHEAD = 1024**2
"""Fixed allowance for a chunk's HDF5 metadata in size estimates."""
MAPPING_OVERHEAD = 256
"""Allowance per step and field for a VDS mapping (about 90 B)."""

StepSources = tuple[int, list[Path], int]
"""(step, distinct source files, their total bytes)."""


class VtkhdfError(Exception):
    """A file or request the tool cannot handle; the message says why."""


@dataclass(frozen=True)
class FieldInfo:
    """One PointData array of an input file.

    kind is "printf" (one %b mapping, as Neko writes), "explicit" (one
    mapping per step) or "inline" (a regular dataset in the main file).
    names holds the stored VDS source file name(s): one for printf, one per
    step for explicit, none for inline. dset is the source dataset name, or
    for inline the dataset's path in the main file.
    """
    name: str
    dtype: np.dtype
    row_shape: tuple[int, ...]
    kind: str
    names: tuple[str, ...]
    dset: str


@dataclass(frozen=True)
class Chunk:
    """One planned split chunk.

    files are the per-step files it references, relative to the input's
    directory, in step order; source_bytes is their total size and
    estimate the planned footprint (chunk file plus files).
    """
    steps: range
    path: Path
    manifest: Path
    files: tuple[str, ...]
    source_bytes: int
    estimate: int


@dataclass
class FileInfo:
    """Validated metadata of an input file.

    A non-temporal file has n_steps == 1 and times None. part_counts maps
    each NumberOf* dataset name to its P per-part entries.
    """
    path: Path
    vtk_type: str
    version: tuple[int, ...]
    temporal: bool
    n_steps: int
    times: np.ndarray | None
    n_parts: int
    n_points: int
    n_cells: int
    n_conn: int
    part_counts: dict[str, np.ndarray]
    fields: list[FieldInfo]
    warnings: list[str]


# --------------------------------------------------------------------------
# Reading and validating


def open_input(path: Path) -> h5py.File:
    """Open an input or per-step source file; the only way either is
    opened."""
    if not path.is_file():
        raise VtkhdfError(f"{path}: no such file")
    try:
        return h5py.File(path, "r", locking="best-effort")
    except OSError as err:
        raise VtkhdfError(f"{path}: cannot open as HDF5: {err}") from None


def inspect(f: h5py.File, path: Path) -> FileInfo:
    """Read and validate the metadata of an open input file."""
    root = f.get("VTKHDF")
    if not isinstance(root, h5py.Group):
        raise VtkhdfError(f"{path}: no /VTKHDF group, not a VTKHDF file")
    vtk_type = attr_text(root, "Type")
    if vtk_type != "UnstructuredGrid":
        raise VtkhdfError(f"{path}: VTKHDF type {vtk_type!r} is not "
                          "supported, only UnstructuredGrid")
    version = tuple(int(v) for v in np.ravel(root.attrs.get("Version", [])))
    warnings: list[str] = []

    steps = root.get("Steps")
    temporal = steps is not None
    if temporal:
        n_steps = int(scalar_attr(steps, "NSteps", path))
        if n_steps < 1:
            raise VtkhdfError(f"{path}: Steps/NSteps is {n_steps}")
        times = head(steps, "Values", n_steps, path).astype(np.float64)
        parts = head(steps, "NumberOfParts", n_steps, path)
        if np.any(parts != parts[0]):
            raise VtkhdfError(f"{path}: NumberOfParts changes between "
                              "steps: time-varying mesh/partitioning not "
                              "supported")
        n_parts = int(parts[0])
        for name in STEP_OFFSETS:
            if np.any(head(steps, name, n_steps, path) != 0):
                raise VtkhdfError(f"{path}: Steps/{name} is non-zero: "
                                  "time-varying mesh/partitioning not "
                                  "supported")
    else:
        n_steps, times = 1, None
        n_parts = dataset(root, "NumberOfPoints", path).shape[0]
    if n_parts < 1:
        raise VtkhdfError(f"{path}: the file has {n_parts} parts")

    part_counts = {name: head(root, name, n_parts, path)
                   for name in COUNT_DATASETS}
    n_points, n_cells, n_conn = (int(part_counts[name].sum())
                                 for name in COUNT_DATASETS)
    expect_shape(root, "Points", (n_points, 3), path)
    expect_shape(root, "Types", (n_cells,), path)
    expect_shape(root, "Offsets", (n_cells + n_parts,), path)
    expect_shape(root, "Connectivity", (n_conn,), path)

    cell_data = root.get("CellData")
    if isinstance(cell_data, h5py.Group) and len(cell_data):
        raise VtkhdfError(f"{path}: CellData is not supported")
    field_data = root.get("FieldData")
    if isinstance(field_data, h5py.Group):
        extra = sorted(set(field_data) - {"Time"})
        if extra:
            warnings.append(f"{path}: input FieldData is ignored: "
                            + ", ".join(extra))

    fields = []
    point_data = root.get("PointData")
    if point_data is not None:
        offsets = steps.get("PointDataOffsets") if temporal else None
        for name in point_data:
            fields.append(inspect_field(point_data, name, offsets, n_steps,
                                        n_points, temporal, path))
    if (os.environ.get("HDF5_VDS_PREFIX")
            and any(fld.kind != "inline" for fld in fields)):
        warnings.append("HDF5_VDS_PREFIX is set: HDF5 and VTK resolve "
                        "relative VDS sources against it first, this tool "
                        "against the input's directory")

    return FileInfo(path, vtk_type, version, temporal, n_steps, times,
                    n_parts, n_points, n_cells, n_conn, part_counts, fields,
                    warnings)


def inspect_field(point_data: h5py.Group, name: str,
                  offsets: h5py.Group | None, n_steps: int, n_points: int,
                  temporal: bool, path: Path) -> FieldInfo:
    """Validate one PointData array and classify its storage."""
    ds = point_data[name]
    where = f"{path}: PointData/{name}"
    if not isinstance(ds, h5py.Dataset):
        raise VtkhdfError(f"{where} is not a dataset")
    if temporal:
        if offsets is None:
            raise VtkhdfError(f"{path}: Steps/PointDataOffsets is missing")
        got = head(offsets, name, n_steps, path)
        if not np.array_equal(got, np.arange(n_steps) * n_points):
            raise VtkhdfError(f"{path}: Steps/PointDataOffsets/{name} is "
                              "not step * NumberOfPoints: time-varying "
                              "mesh/partitioning not supported")

    if not ds.is_virtual:
        shape = ds.shape
        rows = n_steps * n_points
        if not shape or shape[0] < rows or (not temporal
                                           and shape[0] != rows):
            raise VtkhdfError(f"{where} has shape {shape}, expected "
                              f"{rows} rows")
        return FieldInfo(name, ds.dtype, shape[1:], "inline", (), ds.name)
    if not temporal:
        raise VtkhdfError(f"{where} is virtual in a non-temporal file: "
                          "not supported")

    dcpl = ds.id.get_create_plist()
    count = dcpl.get_virtual_count()
    if count == 0:
        raise VtkhdfError(f"{where} is a VDS without mappings")
    names = tuple(dcpl.get_virtual_filename(i) for i in range(count))
    dsets = {dcpl.get_virtual_dsetname(i) for i in range(count)}
    row_shape = tuple(dcpl.get_virtual_vspace(0).shape[1:])
    for stored in names:
        check_vds_name(stored, where)
    if len(dsets) != 1 or "%" in next(iter(dsets)):
        raise VtkhdfError(f"{where}: VDS source dataset names "
                          f"{sorted(dsets)} not supported")
    if any(dcpl.get_virtual_srcspace(i).get_select_type() != h5s.SEL_ALL
           for i in range(count)):
        raise VtkhdfError(f"{where}: VDS source selection is not ALL: "
                          "not supported")

    if count == 1 and "%b" in names[0]:
        if not is_printf_mapping(dcpl.get_virtual_vspace(0), n_points,
                                 row_shape):
            raise VtkhdfError(f"{where}: printf VDS mapping is not "
                              "N rows per step: not supported")
        kind = "printf"
    elif count == n_steps and not any("%b" in n for n in names):
        for i in range(count):
            if not is_step_rows(dcpl.get_virtual_vspace(i), i, n_points,
                                row_shape):
                raise VtkhdfError(f"{where}: VDS mapping {i} does not "
                                  f"cover exactly the rows of step {i}: "
                                  "not supported")
        kind = "explicit"
    else:
        raise VtkhdfError(f"{where}: VDS with {count} mappings for "
                          f"{n_steps} steps is not supported")
    return FieldInfo(name, ds.dtype, row_shape, kind, names, dsets.pop())


def check_vds_name(name: str, where: str) -> None:
    """Refuse a VDS source name with a '%' other than '%b'.

    HDF5 printf-parses every VDS source name, so any other '%' would be
    interpreted differently by HDF5 and by this tool.
    """
    if re.search(r"%(?!b)", name):
        raise VtkhdfError(f"{where}: VDS source name {name!r} contains a "
                          "'%' other than %b: not supported")


def is_printf_mapping(vspace: h5s.SpaceID, n_points: int,
                      row_shape: tuple[int, ...]) -> bool:
    """True for Neko's mapping: blocks of N rows from row 0, unlimited."""
    if (vspace.get_select_type() != h5s.SEL_HYPERSLABS
            or not vspace.is_regular_hyperslab()):
        return False
    start, stride, count, block = vspace.get_regular_hyperslab()
    return (all(s == 0 for s in start)
            and stride[0] == block[0] == n_points
            and count[0] == h5s.UNLIMITED
            and all(c == 1 for c in count[1:])
            and tuple(block[1:]) == row_shape)


def is_step_rows(vspace: h5s.SpaceID, step: int, n_points: int,
                 row_shape: tuple[int, ...]) -> bool:
    """True if the selection is exactly rows [step*N, (step+1)*N)."""
    if vspace.get_select_npoints() != n_points * math.prod(row_shape):
        return False
    low, high = vspace.get_select_bounds()
    return (tuple(vspace.shape[1:]) == row_shape
            and low == (step * n_points,) + (0,) * len(row_shape)
            and high == ((step + 1) * n_points - 1,)
            + tuple(r - 1 for r in row_shape))


def attr_text(obj: h5py.HLObject, name: str) -> str:
    """A string attribute as text, fixed- or variable-length."""
    value = obj.attrs.get(name)
    if isinstance(value, np.ndarray) and value.size == 1:
        value = value.reshape(-1)[0]
    if isinstance(value, bytes):
        value = value.decode("ascii", "replace")
    return "" if value is None else str(value).rstrip("\0")


def scalar_attr(obj: h5py.HLObject, name: str, path: Path):
    """A one-element attribute as a numpy scalar."""
    if name not in obj.attrs:
        raise VtkhdfError(f"{path}: attribute {obj.name}/{name} missing")
    value = np.asarray(obj.attrs[name])
    if value.size != 1:
        raise VtkhdfError(f"{path}: attribute {obj.name}/{name} has "
                          f"{value.size} entries, expected 1")
    return value.reshape(-1)[0]


def dataset(group: h5py.Group, name: str, path: Path) -> h5py.Dataset:
    """A required dataset of group."""
    ds = group.get(name)
    if not isinstance(ds, h5py.Dataset):
        raise VtkhdfError(f"{path}: dataset {group.name}/{name} missing")
    return ds


def head(group: h5py.Group, name: str, n: int, path: Path) -> np.ndarray:
    """The first n entries of a required 1-D metadata dataset."""
    ds = dataset(group, name, path)
    if len(ds.shape) != 1 or ds.shape[0] < n:
        raise VtkhdfError(f"{path}: {ds.name} has shape {ds.shape}, "
                          f"expected at least {n} entries")
    return ds[:n]


def expect_shape(group: h5py.Group, name: str, shape: tuple[int, ...],
                 path: Path) -> None:
    """Require a dataset of exactly the given shape."""
    ds = dataset(group, name, path)
    if ds.shape != shape:
        raise VtkhdfError(f"{path}: {ds.name} has shape {ds.shape}, "
                          f"expected {shape}")


def cell_types(f: h5py.File) -> list[int]:
    """The distinct VTK cell types, read in bounded blocks."""
    ds = f["VTKHDF/Types"]
    found: set[int] = set()
    block = max(1, BLOCK_BYTES // ds.dtype.itemsize)
    for row0 in range(0, ds.shape[0], block):
        found.update(np.unique(ds[row0:row0 + block]).tolist())
    return sorted(found)


# --------------------------------------------------------------------------
# Locating step data


def step_source(info: FileInfo, field: FieldInfo,
                step: int) -> tuple[Path, str, int]:
    """Where the rows of one step live: (file, dataset, first row).

    VDS source names are resolved like HDF5 does without HDF5_VDS_PREFIX:
    absolute as stored, relative against the input file's directory.
    """
    if field.kind == "inline":
        return info.path, field.dset, step * info.n_points
    name = field.names[0] if field.kind == "printf" else field.names[step]
    source = Path(name.replace("%b", str(step)))
    if not source.is_absolute():
        source = info.path.parent / source
    return source, field.dset, 0


def step_files(info: FileInfo, step: int) -> list[Path]:
    """The distinct per-step source files of one step."""
    files = dict.fromkeys(step_source(info, fld, step)[0]
                          for fld in info.fields if fld.kind != "inline")
    return list(files)


def stat_sources(info: FileInfo) -> list[tuple[int, int, list[Path]]]:
    """Per step: (step, bytes of present files, missing files).

    Uses stat only; no source file is opened.
    """
    result = []
    for step in range(info.n_steps):
        nbytes, missing = 0, []
        for source in step_files(info, step):
            try:
                nbytes += source.stat().st_size
            except FileNotFoundError:
                missing.append(source)
        result.append((step, nbytes, missing))
    return result


def check_sources(info: FileInfo,
                  steps: range | list[int]) -> list[StepSources]:
    """Validate the source data of each step.

    Each referenced file is opened read-only, one at a time, and every
    field's dataset in it must have exactly the step's N rows of the right
    row shape: a VDS mapping selects the whole source, and HDF5 refuses to
    read a mapping whose source has a different size. HDF5 would read a
    missing source silently as fill values, so this runs before anything
    is written.
    """
    result = []
    for step in steps:
        by_file: dict[Path, list[tuple[FieldInfo, str]]] = {}
        for fld in info.fields:
            if fld.kind != "inline":
                source, dset, _ = step_source(info, fld, step)
                by_file.setdefault(source, []).append((fld, dset))
        nbytes = 0
        for source, items in by_file.items():
            if not source.is_file():
                raise VtkhdfError(f"step {step}: per-step file {source} "
                                  "is missing")
            with open_input(source) as src:
                for fld, dset in items:
                    ds = src.get(dset)
                    want = (info.n_points,) + fld.row_shape
                    if not isinstance(ds, h5py.Dataset):
                        raise VtkhdfError(f"step {step}: {source} has no "
                                          f"dataset {dset!r}")
                    if ds.shape != want:
                        raise VtkhdfError(f"step {step}: {source}:{dset} "
                                          f"has shape {ds.shape}, expected "
                                          f"{want}")
            nbytes += source.stat().st_size
        result.append((step, list(by_file), nbytes))
    return result


def find_step(info: FileInfo, step: int | None = None,
              time: float | None = None) -> int:
    """Resolve --step (negative counts from the end) or --time (nearest
    step; a tie is an error naming the candidates)."""
    if not info.temporal:
        raise VtkhdfError(f"{info.path} is not temporal (no Steps group): "
                          "it already is a single snapshot")
    n = info.n_steps
    if step is not None:
        index = step + n if step < 0 else step
        if not 0 <= index < n:
            raise VtkhdfError(f"step {step} is out of range: "
                              f"{info.path} has {n} steps (0..{n - 1})")
        return index
    if time is None or not math.isfinite(time):
        raise VtkhdfError(f"invalid time {time}")
    times = info.times
    dist = np.abs(times - time)
    best = int(np.argmin(dist))
    # A tie is two steps equally far away up to float noise, measured
    # against their own spacing; equal times (restarts) always tie.
    ties = np.flatnonzero(dist - dist[best]
                          <= TIME_EPS * np.abs(times - times[best]))
    if len(ties) > 1:
        cands = ", ".join(f"{s} (t = {times[s]:g})" for s in ties)
        raise VtkhdfError(f"time {time:g} is equally close to steps "
                          f"{cands}; choose one with --step")
    return best


def step_label(n_steps: int, step: int) -> str:
    """The step index zero-padded to the width of n_steps - 1."""
    return f"{step:0{len(str(n_steps - 1))}d}"


def strictly_increasing(times: np.ndarray) -> bool:
    return bool(np.all(np.diff(times) > 0))


def time_order_warning(info: FileInfo) -> str | None:
    """A warning naming the first step whose time does not increase."""
    bad = np.flatnonzero(np.diff(info.times) <= 0)
    if not len(bad):
        return None
    s = int(bad[0]) + 1
    return (f"{info.path}: times are not strictly increasing (step {s} has "
            f"t = {info.times[s]:g} after t = {info.times[s - 1]:g}); "
            "each chunk keeps them as they are")


# --------------------------------------------------------------------------
# Chunk planning


def plan_by_count(n_steps: int, k: int) -> list[range]:
    """Consecutive runs of k steps; the last one may be shorter."""
    if k < 1:
        raise VtkhdfError(f"--steps-per-chunk must be at least 1, got {k}")
    return [range(a, min(a + k, n_steps)) for a in range(0, n_steps, k)]


def plan_by_time(times: np.ndarray, width: float) -> list[range]:
    """Group steps by absolute window [m W, (m + 1) W); empty windows are
    skipped. TIME_EPS absorbs float noise at window edges."""
    if not (math.isfinite(width) and width > 0):
        raise VtkhdfError(f"--time-window must be positive, got {width:g}")
    if not strictly_increasing(times):
        raise VtkhdfError("--time-window needs strictly increasing times, "
                          "and this file's are not; use --steps-per-chunk "
                          "or --max-size")
    index = np.floor(times / width + TIME_EPS).astype(np.int64)
    edges = [0, *(np.flatnonzero(np.diff(index)) + 1).tolist(), len(times)]
    return [range(a, b) for a, b in zip(edges[:-1], edges[1:])]


def plan_by_size(step_bytes: list[int], fixed_bytes: int,
                 limit: int) -> list[range]:
    """Greedy packing: add steps while fixed_bytes plus the steps' bytes
    stay within limit. A step that alone exceeds it is an error."""
    chunks, start, used = [], 0, fixed_bytes
    for step, nbytes in enumerate(step_bytes):
        need = fixed_bytes + nbytes
        if need > limit:
            raise VtkhdfError(f"step {step} alone needs about "
                              f"{format_size(need)} ({need} B: mesh, chunk "
                              "overhead and its per-step files), more than "
                              f"--max-size {format_size(limit)} ({limit} B)")
        if used + nbytes > limit:
            chunks.append(range(start, step))
            start, used = step, fixed_bytes
        used += nbytes
    chunks.append(range(start, len(step_bytes)))
    return chunks


def mesh_bytes(f: h5py.File) -> int:
    """Payload bytes of the mesh datasets a chunk copies."""
    root = f["VTKHDF"]
    return sum(root[name].dtype.itemsize * math.prod(root[name].shape)
               for name in MESH_DATASETS)


def source_relpath(info: FileInfo, source: Path) -> str:
    """A per-step file's path relative to the input's directory, as a
    chunk's VDS mapping and manifest store it.

    It must lie inside that directory, or an rsync --files-from manifest
    cannot ship it, and contain no '%', which HDF5 would printf-parse and
    Neko's reader would not.
    """
    rel = os.path.relpath(os.path.abspath(source),
                          os.path.abspath(info.path.parent))
    if ".." in Path(rel).parts:
        raise VtkhdfError(f"per-step file {source} is outside the input's "
                          "directory, so a chunk manifest cannot ship it")
    if "%" in rel:
        raise VtkhdfError(f"per-step file {rel!r} contains '%': not "
                          "supported in a chunk's VDS mapping")
    return rel


def chunk_stem(info: FileInfo, steps: range) -> str:
    first = step_label(info.n_steps, steps[0])
    last = step_label(info.n_steps, steps[-1])
    return f"{info.path.stem}_steps_{first}-{last}"


def prepare_split(f: h5py.File, info: FileInfo, *,
                  steps_per_chunk: int | None = None,
                  time_window: float | None = None,
                  max_size: int | None = None,
                  overwrite: bool = False) -> tuple[list[Chunk], list[str]]:
    """Validate everything a split needs and plan its chunks: returns the
    chunks and any warnings. Writes nothing.

    Exactly one of steps_per_chunk, time_window and max_size selects the
    planning mode. Every per-step file is checked with check_sources, and
    every planned output with check_output, before anything is written.
    """
    if not info.temporal:
        raise VtkhdfError(f"{info.path} is not temporal (no Steps group): "
                          "nothing to split")
    inline = [fld.name for fld in info.fields if fld.kind == "inline"]
    if inline:
        raise VtkhdfError(f"{info.path}: point data {', '.join(inline)} "
                          "are stored inline in the main file, so a chunk "
                          "has nothing to reference; use extract")
    warnings = []
    order = time_order_warning(info)
    if order and time_window is None:
        warnings.append(order)

    checked = check_sources(info, range(info.n_steps))
    sizes = {p: p.stat().st_size for _, files, _ in checked for p in files}
    rel = {p: source_relpath(info, p) for p in sizes}
    step_cost = [sum(sizes[p] for p in files)
                 + MAPPING_OVERHEAD * len(info.fields)
                 for _, files, _ in checked]
    fixed = mesh_bytes(f) + CHUNK_OVERHEAD
    if steps_per_chunk is not None:
        plan = plan_by_count(info.n_steps, steps_per_chunk)
    elif time_window is not None:
        plan = plan_by_time(info.times, time_window)
    else:
        plan = plan_by_size(step_cost, fixed, max_size)

    protected = protected_paths(info)
    chunks = []
    for steps in plan:
        stem = chunk_stem(info, steps)
        path = info.path.parent / f"{stem}.vtkhdf"
        manifest = info.path.parent / f"{stem}.files"
        for out in (path, manifest):
            check_output(out, overwrite, protected)
        files = list(dict.fromkeys(p for s in steps for p in checked[s][1]))
        chunks.append(Chunk(steps, path, manifest,
                            tuple(rel[p] for p in files),
                            sum(sizes[p] for p in files),
                            fixed + sum(step_cost[s] for s in steps)))
    return chunks, warnings


# --------------------------------------------------------------------------
# Writing


def protected_paths(info: FileInfo) -> set[str]:
    """Real paths of the input and of every per-step file it references."""
    sources = {info.path}
    for step in range(info.n_steps):
        sources.update(step_files(info, step))
    return {os.path.realpath(p) for p in sources}


def check_output(path: Path, overwrite: bool, protected: set[str]) -> None:
    """Refuse an output that is an original, exists without overwrite, or
    cannot be created."""
    if os.path.realpath(path) in protected:
        raise VtkhdfError(f"refusing to write {path}: it is the input or "
                          "one of its per-step files")
    if path.is_dir():
        raise VtkhdfError(f"output {path} is a directory")
    if not path.parent.is_dir():
        raise VtkhdfError(f"output directory {path.parent} does not exist")
    if os.path.lexists(path) and not overwrite:
        raise VtkhdfError(f"output {path} already exists (use --overwrite "
                          "to replace it)")


@contextlib.contextmanager
def atomic_output(path: Path, overwrite: bool) -> Iterator[Path]:
    """Yield a temporary path next to path; move it into place on success,
    delete it on any error."""
    tmp = path.with_name(f".{path.name}.tmp-{os.getpid()}")
    try:
        yield tmp
        if not overwrite and os.path.lexists(path):
            raise VtkhdfError(f"output {path} appeared while writing; "
                              "not replaced")
        os.replace(tmp, path)
    except BaseException:
        with contextlib.suppress(FileNotFoundError):
            tmp.unlink()
        raise


def copy_attrs(src: h5py.HLObject, dst: h5py.HLObject,
               skip: tuple[str, ...] = ()) -> None:
    """Copy attributes byte-exactly, keeping each one's stored type.

    h5py's high-level writes would turn fixed-length NULLTERM strings into
    NULLPAD or variable-length ones.
    """
    for name in src.attrs:
        if name in skip:
            continue
        key = name.encode()
        attr = h5a.open(src.id, key)
        ftype, space = attr.get_type(), attr.get_space()
        buf = np.empty(space.shape, dtype=ftype.dtype)
        attr.read(buf, mtype=ftype)
        h5a.create(dst.id, key, ftype, space).write(buf, mtype=ftype)


def copy_rows(src: h5py.Dataset, row0: int, nrows: int, dst: h5py.Dataset,
              dst_row0: int = 0) -> None:
    """Copy rows [row0, row0 + nrows) of src to dst from row dst_row0,
    through one reused buffer of at most BLOCK_BYTES."""
    row_shape = dst.shape[1:]
    row_bytes = dst.dtype.itemsize * math.prod(row_shape)
    block = max(1, BLOCK_BYTES // row_bytes)
    buf = np.empty((min(block, nrows),) + row_shape, dtype=dst.dtype)
    for done in range(0, nrows, block):
        n = min(block, nrows - done)
        src.read_direct(buf, np.s_[row0 + done:row0 + done + n], np.s_[:n])
        dst.write_direct(buf, np.s_[:n],
                         np.s_[dst_row0 + done:dst_row0 + done + n])


def copy_mesh(src_vtkhdf: h5py.Group, dst_vtkhdf: h5py.Group) -> None:
    """Copy the static mesh datasets as contiguous datasets."""
    for name in MESH_DATASETS:
        src = src_vtkhdf[name]
        dst = dst_vtkhdf.create_dataset(name, src.shape, dtype=src.dtype)
        copy_rows(src, 0, src.shape[0], dst)
        copy_attrs(src, dst)


def write_part_counts(info: FileInfo, dst_vtkhdf: h5py.Group,
                      repeats: int) -> None:
    """Write the P per-part NumberOf* entries repeats times, one write per
    repeat (VTK counts pieces as entries / steps)."""
    p = info.n_parts
    for name, counts in info.part_counts.items():
        dst = dst_vtkhdf.create_dataset(name, (repeats * p,),
                                        dtype=counts.dtype)
        for r in range(repeats):
            dst[r * p:(r + 1) * p] = counts


def fill_value(dtype: np.dtype):
    """Fill for a chunk's VDS, read wherever a source is missing: NaN for
    floats, the most negative value for signed and the largest for
    unsigned integers, HDF5's default (zero) otherwise."""
    if dtype.kind in "fc":
        return np.nan
    if dtype.kind == "i":
        return np.iinfo(dtype).min
    if dtype.kind == "u":
        return np.iinfo(dtype).max
    return None


def write_steps(info: FileInfo, steps: range, dst_vtkhdf: h5py.Group) -> None:
    """The chunk's Steps group: original times, P parts, zero offsets and
    PointDataOffsets j * N for its k steps."""
    k = len(steps)
    group = dst_vtkhdf.create_group("Steps")
    group.attrs.create("NSteps", k, dtype=np.int32)
    group.create_dataset("Values", data=info.times[steps.start:steps.stop])
    group.create_dataset("NumberOfParts",
                         data=np.full(k, info.n_parts, dtype=np.int64))
    for name in STEP_OFFSETS:
        group.create_dataset(name, data=np.zeros(k, dtype=np.int64))
    offsets = group.create_group("PointDataOffsets")
    for fld in info.fields:
        offsets.create_dataset(
            fld.name, data=np.arange(k, dtype=np.int64) * info.n_points)


def write_vds(info: FileInfo, field: FieldInfo, steps: range,
              src: h5py.Dataset, dst_group: h5py.Group) -> None:
    """A fixed-shape VDS whose mapping j covers rows [jN, (j+1)N) and
    selects all of the source of step steps[j]; mappings are added in
    ascending j because Neko's reader uses mapping index == step."""
    n, rows = info.n_points, field.row_shape
    shape = (len(steps) * n,) + rows
    dcpl = h5p.create(h5p.DATASET_CREATE)
    # No timestamps in the object header, as h5py's create_dataset does,
    # so that outputs are reproducible.
    dcpl.set_obj_track_times(False)
    fill = fill_value(field.dtype)
    if fill is not None:
        dcpl.set_fill_value(np.array(fill, dtype=field.dtype))
    vspace = h5s.create_simple(shape)
    src_space = h5s.create_simple((n,) + rows)
    for j, step in enumerate(steps):
        source, dset, _ = step_source(info, field, step)
        vspace.select_hyperslab((j * n,) + (0,) * len(rows),
                                (1,) * len(shape), block=(n,) + rows)
        dcpl.set_virtual(vspace, source_relpath(info, source).encode(),
                         dset.encode(), src_space)
    dsid = h5d.create(dst_group.id, field.name.encode(), src.id.get_type(),
                      h5s.create_simple(shape), dcpl=dcpl)
    copy_attrs(src, h5py.Dataset(dsid))


def write_chunk(f: h5py.File, info: FileInfo, steps: range, out: Path,
                overwrite: bool,
                protected: set[str] | None = None) -> list[Path]:
    """Write one split chunk: the mesh, Steps for its k steps and point
    data as VDS referencing the original per-step files; FieldData/Time
    only when k == 1, as VTK reads a 1-step file as static. Returns the
    referenced files. protected defaults to protected_paths(info).
    """
    if any(fld.kind == "inline" for fld in info.fields):
        raise VtkhdfError(f"{info.path}: inline point data cannot be "
                          "referenced by a chunk")
    check_output(out, overwrite, protected if protected is not None
                 else protected_paths(info))
    src_root = f["VTKHDF"]
    with atomic_output(out, overwrite) as tmp, h5py.File(tmp, "x") as dst:
        root = dst.create_group("VTKHDF")
        copy_attrs(src_root, root)
        copy_mesh(src_root, root)
        write_part_counts(info, root, len(steps))
        write_steps(info, steps, root)
        point_data = root.create_group("PointData")
        for fld in info.fields:
            write_vds(info, fld, steps, src_root["PointData"][fld.name],
                      point_data)
        if len(steps) == 1:
            root.create_group("FieldData").create_dataset(
                "Time", data=info.times[steps.start:steps.stop])
        dst.attrs["source_file"] = info.path.name
        dst.attrs["source_steps"] = np.array([steps[0], steps[-1]],
                                             dtype=np.int64)
    return list(dict.fromkeys(p for s in steps for p in step_files(info, s)))


def write_manifest(path: Path, entries: list[str], overwrite: bool) -> None:
    """Write one path per line, atomically."""
    with atomic_output(path, overwrite) as tmp:
        with open(tmp, "x", encoding="utf-8") as handle:
            handle.writelines(f"{entry}\n" for entry in entries)


def write_extract(f: h5py.File, info: FileInfo, step: int, out: Path,
                  overwrite: bool) -> None:
    """Write one step as a standalone non-temporal VTKHDF file: no Steps
    group, no VDS, the step's time kept as FieldData/Time."""
    if not info.temporal:
        raise VtkhdfError(f"{info.path} is not temporal: nothing to "
                          "extract")
    check_sources(info, [step])
    check_output(out, overwrite, protected_paths(info))
    src_root = f["VTKHDF"]
    with atomic_output(out, overwrite) as tmp, h5py.File(tmp, "x") as dst:
        root = dst.create_group("VTKHDF")
        copy_attrs(src_root, root)
        copy_mesh(src_root, root)
        write_part_counts(info, root, 1)
        point_data = root.create_group("PointData")
        for fld in info.fields:
            out_ds = point_data.create_dataset(
                fld.name, (info.n_points,) + fld.row_shape, dtype=fld.dtype)
            copy_attrs(src_root["PointData"][fld.name], out_ds)
            source, dset, row0 = step_source(info, fld, step)
            if fld.kind == "inline":
                copy_rows(f[dset], row0, info.n_points, out_ds)
            else:
                with open_input(source) as src:
                    copy_rows(src[dset], row0, info.n_points, out_ds)
        # A scalar dataspace is rejected by VTK, so Time has shape (1,).
        root.create_group("FieldData").create_dataset(
            "Time", data=np.array([info.times[step]], dtype=np.float64))
        dst.attrs["source_file"] = info.path.name
        dst.attrs["source_step"] = np.int64(step)


# --------------------------------------------------------------------------
# Sizes


def parse_size(text: str) -> int:
    """Bytes from e.g. '512', '40M' or '1.5G' (K/M/G/T are powers of
    1024)."""
    match = re.fullmatch(r"\s*(\d+(?:\.\d*)?)\s*([KMGT]?)\s*", text,
                         re.IGNORECASE)
    if match is None:
        raise VtkhdfError(f"invalid size {text!r}: expected a number with "
                          "an optional K, M, G or T suffix")
    number, unit = match.groups()
    power = SIZE_UNITS.index(unit.upper()) + 1 if unit else 0
    return int(float(number) * 1024**power)


def format_size(nbytes: int) -> str:
    """Bytes as a short human-readable string in powers of 1024."""
    if nbytes < 1024:
        return f"{nbytes} B"
    value = float(nbytes)
    for unit in SIZE_UNITS:
        value /= 1024
        if value < 1024:
            break
    return f"{value:.1f} {unit}iB"
