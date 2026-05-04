#!/usr/bin/env python3

import json
import os
import re
import sys
import time
from dataclasses import dataclass
from typing import Optional

import h5py
import numpy as np
from mpi4py import MPI

from pysemtools.datatypes.coef import Coef
from pysemtools.datatypes.field import FieldRegistry
from pysemtools.datatypes.msh import Mesh
from pysemtools.io.adios2.stream import DataStreamer
from pysemtools.io.ppymech.neksuite import pynekwrite
from pysemtools.io.utils import get_fld_from_ndarray
from pysemtools.rom.io_help import IoHelp
from pysemtools.rom.pod import POD

from neko_communicator import CtrlClient
from neko_communicator import MODE_ADJOINT
from neko_communicator import MODE_FORWARD
from neko_communicator import MODE_STOP
from neko_communicator import PHASE_ADJ_DONE
from neko_communicator import PHASE_ADJ_RUNNING
from neko_communicator import PHASE_FWD_DONE
from neko_communicator import PHASE_FWD_RUNNING
from neko_communicator import get_peer_root
from neko_communicator import make_local_comm

DEBUG = False


def log(comm: MPI.Comm, msg: str) -> None:
    if not DEBUG:
        return
    print(
        f"[py r={comm.Get_rank()}/{comm.Get_size()}] {msg}",
        flush=True,
    )


def as_bool(value, default: bool = False) -> bool:
    if value is None:
        return default
    if isinstance(value, str):
        return value.strip().lower() in ("1", "true", "yes", "y", "t")
    return bool(value)


def strip_json_comments(text: str) -> str:
    out = []
    idx = 0
    n_chars = len(text)
    in_string = False
    escaped = False

    while idx < n_chars:
        char = text[idx]

        if in_string:
            out.append(char)
            if escaped:
                escaped = False
            elif char == "\\":
                escaped = True
            elif char == '"':
                in_string = False
            idx += 1
            continue

        if char == '"':
            in_string = True
            out.append(char)
            idx += 1
            continue

        if char == "/" and idx + 1 < n_chars:
            nxt = text[idx + 1]
            if nxt == "/":
                idx += 2
                while idx < n_chars and text[idx] != "\n":
                    idx += 1
                continue
            if nxt == "*":
                idx += 2
                while idx + 1 < n_chars:
                    if text[idx] == "*" and text[idx + 1] == "/":
                        break
                    idx += 1
                idx += 2
                continue

        out.append(char)
        idx += 1

    return "".join(out)


def load_case_config(case_path: str) -> dict:
    with open(case_path, "r", encoding="utf-8") as handle:
        text = handle.read()
    cleaned = strip_json_comments(text)
    cleaned = re.sub(r",\s*([}\]])", r"\1", cleaned)
    return json.loads(cleaned)


def rotate_time_coeffs(path: str) -> Optional[str]:
    if not os.path.exists(path):
        return None

    stem, ext = os.path.splitext(path)
    idx = 1
    while True:
        candidate = f"{stem}_{idx:03d}{ext}"
        if not os.path.exists(candidate):
            os.rename(path, candidate)
            return candidate
        idx += 1


def count_enabled_scalars(case_cfg: dict) -> int:
    case_data = case_cfg.get("case", {})
    if not isinstance(case_data, dict):
        return 0

    scalar_cfg = case_data.get("scalar")
    if isinstance(scalar_cfg, dict):
        return int(as_bool(scalar_cfg.get("enabled"), True))

    scalar_cfgs = case_data.get("scalars")
    if not isinstance(scalar_cfgs, list):
        return 0

    return sum(
        1
        for entry in scalar_cfgs
        if not isinstance(entry, dict)
        or as_bool(entry.get("enabled"), True)
    )


@dataclass(frozen=True)
class PODConfig:
    batch_size: int
    keep_modes: int
    i_stream: int
    dtype: type
    write_modes: bool
    mode_output_dtype: type
    mode_output_wdsz: int
    mode_output_format: str
    mode_output_file_name: str
    include_scalar: bool
    debug: bool
    snapshot_dt: float

    @property
    def n_fields(self) -> int:
        return 4 if self.include_scalar else 3


@dataclass
class EnergyState:
    total: float = 0.0
    count: int = 0

    def add(
        self, comm: MPI.Comm, fields: list[np.ndarray], bm: np.ndarray
    ) -> None:
        local_energy = 0.0
        for field in fields:
            local_energy += np.sum(field * field * bm, dtype=np.float64)
        self.total += float(comm.allreduce(local_energy, op=MPI.SUM))
        self.count += 1


def load_pod_config(case_path: str) -> tuple[dict, PODConfig]:
    case = load_case_config(case_path)
    state_recovery = case.get("state_recovery", {})
    if not isinstance(state_recovery, dict):
        raise KeyError("case['state_recovery'] must be a dict")

    enabled_scalars = count_enabled_scalars(case)
    include_scalar = enabled_scalars > 0
    if enabled_scalars > 1:
        raise ValueError(
            "POD state recovery currently supports exactly one enabled "
            "scalar in the case."
        )

    dtype_name = str(state_recovery.get("dtype", "double")).strip().lower()
    if dtype_name == "single":
        dtype = np.float32
    elif dtype_name == "double":
        dtype = np.float64
    else:
        raise ValueError(f"Unsupported dtype '{dtype_name}'")

    output_precision = str(
        state_recovery.get("output_precision", "sp")
    ).strip().lower()
    if output_precision in ("sp", "single"):
        mode_output_dtype = np.float32
        mode_output_wdsz = 4
    elif output_precision in ("dp", "double"):
        mode_output_dtype = np.float64
        mode_output_wdsz = 8
    else:
        raise ValueError(
            f"Unsupported output_precision '{output_precision}'"
        )

    mode_output_format = str(
        state_recovery.get("output_format", "fld")
    ).strip().lower()
    mode_output_file_name = str(
        state_recovery.get("output_file_name", "POD_modes")
    ).strip()

    i_stream = int(state_recovery["i_stream"])
    timestep = float(case["case"]["time"]["timestep"])

    cfg = PODConfig(
        batch_size=int(state_recovery["batch_size"]),
        keep_modes=int(state_recovery["n_modes"]),
        i_stream=i_stream,
        dtype=dtype,
        write_modes=as_bool(state_recovery.get("write_modes")),
        mode_output_dtype=mode_output_dtype,
        mode_output_wdsz=mode_output_wdsz,
        mode_output_format=mode_output_format,
        mode_output_file_name=mode_output_file_name,
        include_scalar=include_scalar,
        debug=as_bool(state_recovery.get("debug")),
        snapshot_dt=timestep * i_stream,
    )
    return case, cfg

def recv_field(ds: DataStreamer, dtype: type) -> np.ndarray:
    field = get_fld_from_ndarray(ds.recieve(), ds.lx, ds.ly, ds.lz, ds.nelv)
    return field.astype(dtype)


def recv_fields(
    ds: DataStreamer, n_fields: int, dtype: type
) -> list[np.ndarray]:
    return [recv_field(ds, dtype) for _ in range(n_fields)]


def make_pod(
    comm: MPI.Comm,
    bm: np.ndarray,
    cfg: PODConfig,
) -> tuple[POD, IoHelp]:
    pod = POD(
        comm,
        number_of_modes_to_update=cfg.keep_modes,
        global_updates=True,
        auto_expand=False,
        bckend="numpy",
    )
    ioh = IoHelp(
        comm,
        number_of_fields=cfg.n_fields,
        batch_size=cfg.batch_size,
        field_size=bm.size,
        field_data_type=cfg.dtype,
        mass_matrix_data_type=cfg.dtype,
    )
    mass_fields = [np.copy(np.sqrt(bm)) for _ in range(cfg.n_fields)]
    ioh.copy_fieldlist_to_xi(mass_fields)
    ioh.bm1sqrt[:, :] = np.copy(ioh.xi[:, :])
    return pod, ioh


def flush_buffer(comm: MPI.Comm, pod: POD, ioh: IoHelp) -> None:
    if ioh.buffer_index > 0 and not ioh.update_from_buffer:
        log(comm, f"POD.update(buff) buffer_index={ioh.buffer_index}")
        pod.update(comm, buff=ioh.buff[:, :ioh.buffer_index])
        ioh.buffer_index = 0


def add_snapshot(
    comm: MPI.Comm,
    pod: POD,
    ioh: IoHelp,
    bm: np.ndarray,
    fields: list[np.ndarray],
    tcur: Optional[float],
    times: list[float],
    energy: EnergyState,
) -> None:
    if tcur is not None:
        times.append(float(tcur))

    energy.add(comm, fields, bm)
    ioh.copy_fieldlist_to_xi(fields)
    ioh.load_buffer(scale_snapshot=True)

    if ioh.update_from_buffer:
        log(comm, f"POD.update(buff) buffer_index={ioh.buffer_index}")
        pod.update(comm, buff=ioh.buff[:, :ioh.buffer_index])
        ioh.buffer_index = 0
        ioh.update_from_buffer = False


def available_modes(pod: POD, keep_modes: int) -> int:
    local_modes = getattr(pod, "u_1t", None)
    if local_modes is None:
        return 0
    return min(keep_modes, local_modes.shape[1])


def build_time_coefficients(
    comm: MPI.Comm,
    pod: POD,
    cfg: PODConfig,
    times: list[float],
) -> tuple[np.ndarray, str, int]:
    n_avail = available_modes(pod, cfg.keep_modes)
    if n_avail > 0:
        coeffs = (
            pod.d_1t[:n_avail, None] * pod.vt_1t[:n_avail, :]
        ).T
    else:
        coeffs = np.zeros((0, 0), dtype=np.float64)

    n_snaps = coeffs.shape[0]
    if len(times) == n_snaps:
        time_vec = np.asarray(times, dtype=np.float64)
    else:
        if times:
            log(
                comm,
                f"time list length {len(times)} != n_snaps {n_snaps}; "
                "falling back to dt",
            )
        time_vec = np.arange(n_snaps, dtype=np.float64) * cfg.snapshot_dt

    if coeffs.shape[1] < cfg.keep_modes:
        coeffs = np.hstack(
            [
                coeffs,
                np.zeros(
                    (n_snaps, cfg.keep_modes - coeffs.shape[1]),
                    dtype=coeffs.dtype,
                ),
            ]
        )

    header = "t," + ",".join(f"a{i + 1}" for i in range(cfg.keep_modes))
    return np.column_stack([time_vec, coeffs]), header, n_avail


def write_time_coefficients(
    comm: MPI.Comm,
    data: np.ndarray,
    header: str,
    rotate_existing: bool,
) -> None:
    if comm.Get_rank() != 0:
        return

    if rotate_existing:
        rotated = rotate_time_coeffs("pod_time_coeffs.csv")
        if rotated:
            log(comm, f"archived pod_time_coeffs.csv -> {rotated}")

    np.savetxt(
        "pod_time_coeffs.csv",
        data,
        delimiter=",",
        header=header,
        comments="",
    )


def report_energy_capture(
    pod: POD,
    energy: EnergyState,
    keep_modes: int,
) -> None:
    if energy.total <= 0.0 or energy.count == 0:
        print("POD energy capture: no snapshot energy available.", flush=True)
        return

    singular_values = getattr(pod, "d_1t", None)
    if singular_values is None or singular_values.size == 0:
        print("POD energy capture: no singular values available.", flush=True)
        return

    energies = np.asarray(singular_values, dtype=np.float64) ** 2
    fractions = 100.0 * energies / energy.total
    n_report = min(keep_modes, fractions.size)

    print("POD energy capture (% of snapshot energy):", flush=True)
    for idx in range(n_report):
        print(f"  mode {idx + 1}: {fractions[idx]:.6f}%", flush=True)
    print(
        f"  sum first {n_report}: "
        f"{np.sum(fractions[:n_report]):.6f}%",
        flush=True,
    )


def stream_mode_fields(
    ds: DataStreamer,
    ioh: IoHelp,
    pod: POD,
    cfg: PODConfig,
    zero_field: np.ndarray,
    n_avail: int,
) -> None:
    for idx in range(cfg.keep_modes):
        for field in mode_fields_1d(
            ioh, pod, cfg, zero_field, n_avail, idx
        ):
            ds.stream(field.astype(cfg.dtype, copy=False))


def mode_fields_1d(
    ioh: IoHelp,
    pod: POD,
    cfg: PODConfig,
    zero_field: np.ndarray,
    n_avail: int,
    idx: int,
) -> list[np.ndarray]:
    if idx < n_avail:
        return ioh.split_narray_to_1dfields(pod.u_1t[:, idx])
    return [zero_field for _ in range(cfg.n_fields)]


def mode_output_base(cfg: PODConfig) -> tuple[str, str, str]:
    file_name = os.path.normpath(cfg.mode_output_file_name)
    out_dir = os.path.dirname(file_name)
    out_stem = os.path.splitext(os.path.basename(file_name))[0]
    sample_base = os.path.join(out_dir, f"{out_stem}0")
    meta_path = os.path.join(out_dir, f"{out_stem}0.nek5000")
    return out_dir, sample_base, meta_path


def ensure_mode_output_dir(comm: MPI.Comm, out_dir: str) -> None:
    if comm.Get_rank() == 0 and out_dir:
        os.makedirs(out_dir, exist_ok=True)
    comm.Barrier()


def build_mode_output_fields(
    comm: MPI.Comm,
    ioh: IoHelp,
    pod: POD,
    cfg: PODConfig,
    zero_field: np.ndarray,
    n_avail: int,
    output_index: int,
) -> FieldRegistry:
    flat_fields = []
    for mode_idx in range(cfg.keep_modes):
        flat_fields.extend(
            mode_fields_1d(ioh, pod, cfg, zero_field, n_avail, mode_idx)
        )

    n_total = len(flat_fields)
    if n_total < 1 or n_total > 99:
        raise ValueError(
            f"Unsupported number of POD output fields: {n_total}"
        )

    if n_total == 1:
        field_names = ["p"]
    elif n_total == 2:
        field_names = ["p", "t"]
    elif n_total == 3:
        field_names = ["u", "v", "w"]
    elif n_total == 4:
        field_names = ["p", "u", "v", "w"]
    else:
        field_names = ["p", "u", "v", "w", "t"]
        field_names.extend(f"s{i}" for i in range(n_total - 5))

    fld = FieldRegistry(comm)
    fld.t = float(output_index)
    for name, field in zip(field_names, flat_fields):
        fld.add_field(
            comm,
            field_name=name,
            field=np.asarray(field, dtype=cfg.mode_output_dtype),
            dtype=cfg.mode_output_dtype,
        )

    return fld


def write_nek_index_file(
    comm: MPI.Comm,
    sample_base: str,
    meta_path: str,
    n_steps: int,
) -> None:
    if comm.Get_rank() != 0:
        return

    series_name = os.path.basename(sample_base[:-1])
    with open(meta_path, "w", encoding="utf-8") as handle:
        handle.write(
            f"filetemplate:         {series_name}%01d.f%05d\n"
        )
        handle.write(f"firsttimestep: {0:5d}\n")
        handle.write(f"numtimesteps: {n_steps:5d}\n")


def vtk_mode_field_names(cfg: PODConfig) -> list[str]:
    names = []
    for mode_idx in range(cfg.keep_modes):
        names.extend(
            [
                f"u_mode_{mode_idx + 1}",
                f"v_mode_{mode_idx + 1}",
                f"w_mode_{mode_idx + 1}",
            ]
        )
        if cfg.include_scalar:
            names.append(f"s_mode_{mode_idx + 1}")
    return names


def vtk_output_path(cfg: PODConfig) -> tuple[str, str]:
    file_name = os.path.normpath(cfg.mode_output_file_name)
    out_dir = os.path.dirname(file_name)
    out_stem = os.path.splitext(os.path.basename(file_name))[0]
    vtk_path = os.path.join(out_dir, f"{out_stem}_0.vtkhdf")
    return out_dir, vtk_path


def vtk_quad_ordering(lx: int, ly: int) -> np.ndarray:
    ordering = np.empty(lx * ly, dtype=np.int64)
    n_corners = 4
    n_edges = 2 * ((lx - 2) + (ly - 2))

    for j in range(1, ly + 1):
        for i in range(1, lx + 1):
            ibdy = int(i == 1 or i == lx)
            jbdy = int(j == 1 or j == ly)
            nbdy = ibdy + jbdy

            if nbdy == 2:
                if i == 1:
                    vtk_idx = 0 if j == 1 else 3
                else:
                    vtk_idx = 1 if j == 1 else 2
            elif nbdy == 1:
                offset = n_corners
                if ibdy == 0:
                    vtk_idx = (i - 2) + (lx - 2 + ly - 2 if j != 1 else 0)
                    vtk_idx += offset
                else:
                    vtk_idx = (j - 2)
                    vtk_idx += lx - 2 if i != 1 else 2 * (lx - 2) + (ly - 2)
                    vtk_idx += offset
            else:
                vtk_idx = n_corners + n_edges + (i - 2) + (lx - 2) * (j - 2)

            ordering[vtk_idx] = (i - 1) + lx * (j - 1)

    return ordering


def vtk_hex_ordering(lx: int, ly: int, lz: int) -> np.ndarray:
    ordering = np.empty(lx * ly * lz, dtype=np.int64)
    n_corners = 8
    n_edges = 4 * ((lx - 2) + (ly - 2) + (lz - 2))
    n_faces = 2 * (
        (lx - 2) * (ly - 2)
        + (lx - 2) * (lz - 2)
        + (ly - 2) * (lz - 2)
    )

    for k in range(1, lz + 1):
        for j in range(1, ly + 1):
            for i in range(1, lx + 1):
                ibdy = int(i == 1 or i == lx)
                jbdy = int(j == 1 or j == ly)
                kbdy = int(k == 1 or k == lz)
                nbdy = ibdy + jbdy + kbdy

                if nbdy == 3:
                    if i == 1:
                        vtk_idx = 0 if j == 1 else 3
                    else:
                        vtk_idx = 1 if j == 1 else 2
                    vtk_idx += 0 if k == 1 else 4
                elif nbdy == 2:
                    offset = n_corners
                    if ibdy == 0:
                        vtk_idx = (i - 2)
                        if j != 1:
                            vtk_idx += (lx - 2) + (ly - 2)
                        if k != 1:
                            vtk_idx += 2 * ((lx - 2) + (ly - 2))
                        vtk_idx += offset
                    elif jbdy == 0:
                        vtk_idx = (j - 2)
                        if i != 1:
                            vtk_idx += lx - 2
                        else:
                            vtk_idx += 2 * (lx - 2) + (ly - 2)
                        if k != 1:
                            vtk_idx += 2 * ((lx - 2) + (ly - 2))
                        vtk_idx += offset
                    else:
                        vtk_idx = (k - 2)
                        if i == 1:
                            vtk_idx += 0 if j == 1 else 3 * (lz - 2)
                        else:
                            vtk_idx += lz - 2 if j == 1 else 2 * (lz - 2)
                        vtk_idx += offset + 4 * ((lx - 2) + (ly - 2))
                elif nbdy == 1:
                    offset = n_corners + n_edges
                    if ibdy == 1:
                        vtk_idx = (j - 2) + (ly - 2) * (k - 2)
                        if i != 1:
                            vtk_idx += (ly - 2) * (lz - 2)
                        vtk_idx += offset
                    elif jbdy == 1:
                        offset += 2 * (ly - 2) * (lz - 2)
                        vtk_idx = (i - 2) + (lx - 2) * (k - 2)
                        if j != 1:
                            vtk_idx += (lx - 2) * (lz - 2)
                        vtk_idx += offset
                    else:
                        offset += 2 * (ly - 2) * (lz - 2)
                        offset += 2 * (lx - 2) * (lz - 2)
                        vtk_idx = (i - 2) + (lx - 2) * (j - 2)
                        if k != 1:
                            vtk_idx += (lx - 2) * (ly - 2)
                        vtk_idx += offset
                else:
                    offset = n_corners + n_edges + n_faces
                    vtk_idx = offset + (i - 2)
                    vtk_idx += (lx - 2) * ((j - 2) + (ly - 2) * (k - 2))

                ordering[vtk_idx] = (i - 1) + lx * ((j - 1) + ly * (k - 1))

    return ordering


def vtk_points_local(mesh: Mesh) -> np.ndarray:
    return np.column_stack(
        [
            np.asarray(mesh.x, dtype=np.float64).reshape(-1),
            np.asarray(mesh.y, dtype=np.float64).reshape(-1),
            np.asarray(mesh.z, dtype=np.float64).reshape(-1),
        ]
    )


def vtk_connectivity_local(
    mesh: Mesh,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    if mesh.gdim == 3:
        vtk_type = 72
        ordering = vtk_hex_ordering(mesh.lx, mesh.ly, mesh.lz)
    else:
        vtk_type = 70
        ordering = vtk_quad_ordering(mesh.lx, mesh.ly)

    points_per_elem = int(mesh.lx * mesh.ly * mesh.lz)
    conn_per_elem = int(ordering.size)
    local_cells = int(mesh.nelv)

    connectivity = np.empty(local_cells * conn_per_elem, dtype=np.int32)
    for elem in range(local_cells):
        start = elem * conn_per_elem
        connectivity[start : start + conn_per_elem] = (
            elem * points_per_elem + ordering
        ).astype(np.int32)

    offsets = np.arange(
        0,
        (local_cells + 1) * conn_per_elem,
        conn_per_elem,
        dtype=np.int32,
    )
    types = np.full(local_cells, vtk_type, dtype=np.uint8)
    return connectivity, offsets, types


def gather_vtk_static(
    comm: MPI.Comm,
    mesh: Mesh,
) -> tuple[
    Optional[np.ndarray],
    Optional[np.ndarray],
    Optional[np.ndarray],
    Optional[np.ndarray],
]:
    local_points = vtk_points_local(mesh)
    local_conn, local_offsets, local_types = vtk_connectivity_local(mesh)

    gathered_points = comm.gather(local_points, root=0)
    gathered_conn = comm.gather(local_conn, root=0)
    gathered_offsets = comm.gather(local_offsets, root=0)
    gathered_types = comm.gather(local_types, root=0)

    if comm.Get_rank() != 0:
        return None, None, None, None

    point_sizes = [pts.shape[0] for pts in gathered_points]
    point_offsets = np.cumsum([0] + point_sizes[:-1], dtype=np.int64)

    global_points = np.concatenate(gathered_points, axis=0)

    shifted_conn = []
    for conn_local, point_offset in zip(gathered_conn, point_offsets):
        shifted_conn.append(conn_local.astype(np.int64) + point_offset)
    global_conn = np.concatenate(shifted_conn).astype(np.int32)

    conn_per_elem = None
    for offsets_local in gathered_offsets:
        if offsets_local.size > 1:
            conn_per_elem = int(offsets_local[1] - offsets_local[0])
            break
    if conn_per_elem is None:
        raise ValueError("Could not determine VTK connectivity stride.")

    total_cells = sum(types.shape[0] for types in gathered_types)
    global_offsets = np.arange(
        0,
        (total_cells + 1) * conn_per_elem,
        conn_per_elem,
        dtype=np.int32,
    )
    global_types = np.concatenate(gathered_types).astype(np.uint8)
    return global_points, global_conn, global_offsets, global_types


def gather_vtk_step_fields(
    comm: MPI.Comm,
    field_names: list[str],
    ioh: IoHelp,
    pod: POD,
    cfg: PODConfig,
    zero_field: np.ndarray,
    n_avail: int,
) -> Optional[dict[str, np.ndarray]]:
    local_fields = {}
    flat_fields = []
    for mode_idx in range(cfg.keep_modes):
        flat_fields.extend(
            mode_fields_1d(ioh, pod, cfg, zero_field, n_avail, mode_idx)
        )

    for name, field in zip(field_names, flat_fields):
        local_fields[name] = np.asarray(
            field, dtype=cfg.mode_output_dtype
        ).reshape(-1)

    gathered = comm.gather(local_fields, root=0)
    if comm.Get_rank() != 0:
        return None

    merged = {}
    for name in field_names:
        merged[name] = np.concatenate(
            [rank_fields[name] for rank_fields in gathered]
        )
    return merged


def vtk_require_dataset(
    group: h5py.Group,
    name: str,
    dtype,
    shape: tuple[int, ...],
    maxshape: tuple[Optional[int], ...],
) -> h5py.Dataset:
    if name in group:
        return group[name]
    return group.create_dataset(
        name,
        shape=shape,
        maxshape=maxshape,
        chunks=True,
        dtype=dtype,
    )


def vtk_append_scalar(
    dataset: h5py.Dataset,
    values: np.ndarray,
) -> None:
    old_size = dataset.shape[0]
    dataset.resize((old_size + values.shape[0],))
    dataset[old_size:] = values


def vtk_append_steps(
    vtk_group: h5py.Group,
    field_names: list[str],
    output_index: int,
    n_points: int,
    n_cells: int,
    n_conn: int,
) -> None:
    steps = vtk_group.require_group("Steps")
    values = vtk_require_dataset(
        steps,
        "Values",
        np.float64,
        (0,),
        (None,),
    )
    n_parts = vtk_require_dataset(
        steps,
        "NumberOfParts",
        np.int64,
        (0,),
        (None,),
    )
    part_offsets = vtk_require_dataset(
        steps,
        "PartOffsets",
        np.int64,
        (0,),
        (None,),
    )
    point_offsets = vtk_require_dataset(
        steps,
        "PointOffsets",
        np.int64,
        (0,),
        (None,),
    )
    cell_offsets = vtk_require_dataset(
        steps,
        "CellOffsets",
        np.int64,
        (0,),
        (None,),
    )
    conn_offsets = vtk_require_dataset(
        steps,
        "ConnectivityIdOffsets",
        np.int64,
        (0,),
        (None,),
    )

    point_data_offsets = steps.require_group("PointDataOffsets")

    for ds, value in [
        (values, np.array([float(output_index)], dtype=np.float64)),
        (n_parts, np.array([1], dtype=np.int64)),
        (part_offsets, np.array([0], dtype=np.int64)),
        (point_offsets, np.array([0], dtype=np.int64)),
        (cell_offsets, np.array([0], dtype=np.int64)),
        (conn_offsets, np.array([0], dtype=np.int64)),
    ]:
        vtk_append_scalar(ds, value)

    base_offset = np.int64(output_index * n_points)
    for name in field_names:
        ds = vtk_require_dataset(
            point_data_offsets,
            name,
            np.int64,
            (0,),
            (None,),
        )
        vtk_append_scalar(ds, np.array([base_offset], dtype=np.int64))

    for name, value in [
        ("NumberOfPoints", np.array([n_points], dtype=np.int64)),
        ("NumberOfCells", np.array([n_cells], dtype=np.int64)),
        ("NumberOfConnectivityIds", np.array([n_conn], dtype=np.int64)),
    ]:
        ds = vtk_require_dataset(
            vtk_group,
            name,
            np.int64,
            (0,),
            (None,),
        )
        vtk_append_scalar(ds, value)


def write_modes_to_vtkhdf(
    comm: MPI.Comm,
    mesh: Mesh,
    ioh: IoHelp,
    pod: POD,
    cfg: PODConfig,
    zero_field: np.ndarray,
    n_avail: int,
    output_index: int,
) -> int:
    out_dir, vtk_path = vtk_output_path(cfg)
    ensure_mode_output_dir(comm, out_dir)

    field_names = vtk_mode_field_names(cfg)
    global_fields = gather_vtk_step_fields(
        comm,
        field_names,
        ioh,
        pod,
        cfg,
        zero_field,
        n_avail,
    )

    if output_index == 0:
        vtk_static = gather_vtk_static(comm, mesh)
    else:
        vtk_static = (None, None, None, None)

    if comm.Get_rank() == 0:
        if output_index == 0:
            points, conn, offsets, types = vtk_static
            with h5py.File(vtk_path, "w") as h5file:
                vtk_group = h5file.create_group("VTKHDF")
                vtk_group.attrs["Version"] = np.array([2, 6], dtype=np.int32)
                vtk_group.attrs["Type"] = np.bytes_("UnstructuredGrid")
                vtk_group.create_dataset(
                    "Points",
                    data=points,
                    dtype=np.float64,
                )
                vtk_group.create_dataset(
                    "Connectivity",
                    data=conn,
                    dtype=np.int32,
                )
                vtk_group.create_dataset(
                    "Offsets",
                    data=offsets,
                    dtype=np.int32,
                )
                vtk_group.create_dataset(
                    "Types",
                    data=types,
                    dtype=np.uint8,
                )
                point_data = vtk_group.create_group("PointData")
                for name in field_names:
                    point_data.create_dataset(
                        name,
                        shape=(0,),
                        maxshape=(None,),
                        chunks=True,
                        dtype=cfg.mode_output_dtype,
                    )

        with h5py.File(vtk_path, "r+") as h5file:
            vtk_group = h5file["VTKHDF"]
            point_data = vtk_group["PointData"]
            n_points = int(vtk_group["Points"].shape[0])
            n_cells = int(vtk_group["Types"].shape[0])
            n_conn = int(vtk_group["Connectivity"].shape[0])

            for name in field_names:
                vtk_append_scalar(point_data[name], global_fields[name])

            vtk_append_steps(
                vtk_group,
                field_names,
                output_index,
                n_points,
                n_cells,
                n_conn,
            )

    comm.Barrier()
    return output_index + 1


def write_modes_to_disk(
    comm: MPI.Comm,
    mesh: Mesh,
    ioh: IoHelp,
    pod: POD,
    cfg: PODConfig,
    zero_field: np.ndarray,
    n_avail: int,
    output_index: int,
) -> int:
    if not cfg.write_modes:
        return output_index

    fmt = cfg.mode_output_format
    if fmt in ("fld", "nek5000"):
        out_dir, sample_base, meta_path = mode_output_base(cfg)
        ensure_mode_output_dir(comm, out_dir)

        fld = build_mode_output_fields(
            comm,
            ioh,
            pod,
            cfg,
            zero_field,
            n_avail,
            output_index,
        )

        pynekwrite(
            f"{sample_base}.f{output_index:05d}",
            comm,
            msh=mesh,
            fld=fld,
            wdsz=cfg.mode_output_wdsz,
            istep=output_index,
            write_mesh=output_index == 0,
        )
        comm.Barrier()
        write_nek_index_file(
            comm,
            sample_base,
            meta_path,
            output_index + 1,
        )
        comm.Barrier()
        return output_index + 1

    if fmt == "vtkhdf":
        return write_modes_to_vtkhdf(
            comm,
            mesh,
            ioh,
            pod,
            cfg,
            zero_field,
            n_avail,
            output_index,
        )

    raise ValueError(f"Unsupported output_format '{fmt}'")


def init_runtime(
    comm: MPI.Comm,
    cfg: PODConfig,
) -> tuple[
    DataStreamer,
    Mesh,
    np.ndarray,
    POD,
    IoHelp,
    list[np.ndarray],
]:
    ds = DataStreamer(comm)

    x = recv_field(ds, cfg.dtype)
    y = recv_field(ds, cfg.dtype)
    z = recv_field(ds, cfg.dtype)

    initial_fields = recv_fields(ds, cfg.n_fields, cfg.dtype)

    mesh = Mesh(comm, x=x, y=y, z=z, create_connectivity=False)
    coef = Coef(mesh, comm)
    pod, ioh = make_pod(comm, coef.B, cfg)
    return ds, mesh, coef.B, pod, ioh, initial_fields


def main() -> None:
    global DEBUG

    world = MPI.COMM_WORLD
    comm = make_local_comm(world)
    rank = comm.Get_rank()
    peer_root = get_peer_root()

    case_path = sys.argv[1] if len(sys.argv) > 1 else "POD_rugby_ball.case"
    case_path = os.path.abspath(case_path)
    case, cfg = load_pod_config(case_path)
    DEBUG = cfg.debug

    ds, mesh, bm, pod, ioh, initial_fields = init_runtime(comm, cfg)
    times = []
    energy = EnergyState()
    add_snapshot(comm, pod, ioh, bm, initial_fields, 0.0, times, energy)
    mode_output_index = 0

    ctrl = None
    if rank == 0:
        ctrl = CtrlClient(world, peer_root, cfg.debug)

    my_count = getattr(ds, "py2f_field_my_count", None)
    if my_count is None:
        my_count = initial_fields[0].size
    zero_field = np.zeros(my_count, dtype=cfg.dtype)

    try:
        while True:
            if rank == 0:
                state = ctrl.read_state()
            else:
                state = (None, None, None, None)

            mode, phase, step, tcur = comm.bcast(state, root=0)

            if mode is None or mode == MODE_STOP:
                break

            if mode == MODE_FORWARD and phase == PHASE_FWD_RUNNING:
                fields = recv_fields(ds, cfg.n_fields, cfg.dtype)
                add_snapshot(comm, pod, ioh, bm, fields, tcur, times, energy)
                continue

            if mode == MODE_FORWARD and phase == PHASE_FWD_DONE:
                flush_buffer(comm, pod, ioh)
                pod.scale_modes(comm, bm1sqrt=ioh.bm1sqrt, op="div")
                pod.rotate_local_modes_to_global(comm)

                coeffs, header, n_avail = build_time_coefficients(
                    comm, pod, cfg, times
                )
                write_time_coefficients(
                    comm, coeffs, header, cfg.write_modes
                )
                if rank == 0:
                    report_energy_capture(pod, energy, cfg.keep_modes)
                    ctrl.send_cmd(MODE_ADJOINT, PHASE_ADJ_RUNNING)

                stream_mode_fields(ds, ioh, pod, cfg, zero_field, n_avail)
                mode_output_index = write_modes_to_disk(
                    comm,
                    mesh,
                    ioh,
                    pod,
                    cfg,
                    zero_field,
                    n_avail,
                    mode_output_index,
                )
                continue

            if mode == MODE_ADJOINT and phase == PHASE_ADJ_DONE:
                pod, ioh = make_pod(comm, bm, cfg)
                times = []
                energy = EnergyState()
                add_snapshot(
                    comm,
                    pod,
                    ioh,
                    bm,
                    initial_fields,
                    0.0,
                    times,
                    energy,
                )
                continue

            time.sleep(0.001)

    finally:
        if rank == 0 and ctrl is not None:
            ctrl.close()
        comm.Barrier()
        ds.finalize()


if __name__ == "__main__":
    main()
