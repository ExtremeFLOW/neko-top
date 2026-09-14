#!/usr/bin/env python3

import json
import os
import re
from dataclasses import dataclass
from typing import Optional

import numpy as np
from mpi4py import MPI

from pysemtools.datatypes.msh import Mesh
from pysemtools.io.adios2.stream import DataStreamer
from pysemtools.io.hdf.vtkhdf import VTKHDFFile
from pysemtools.io.utils import get_fld_from_ndarray
from pysemtools.io.wrappers import write_data
from pysemtools.rom.io_help import IoHelp
from pysemtools.rom.pod import POD

DEBUG = False


def set_debug(enabled: bool) -> None:
    global DEBUG
    DEBUG = enabled


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


def rotate_existing_file(path: str) -> Optional[str]:
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
    elif output_precision in ("dp", "double"):
        mode_output_dtype = np.float64
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
        rotated = rotate_existing_file("pod_time_coeffs.csv")
        if rotated:
            log(comm, f"archived pod_time_coeffs.csv -> {rotated}")

    np.savetxt(
        "pod_time_coeffs.csv",
        data,
        delimiter=",",
        header=header,
        comments="",
    )


def singular_value_header(keep_modes: int) -> str:
    return "iteration,total_energy," + ",".join(
        f"sigma{i + 1}" for i in range(keep_modes)
    )


def singular_value_row(
    pod: POD,
    keep_modes: int,
    iteration: int,
) -> np.ndarray:
    singular_values = np.asarray(
        getattr(pod, "d_1t", np.zeros((0,), dtype=np.float64)),
        dtype=np.float64,
    ).reshape(-1)

    if singular_values.size < keep_modes:
        singular_values = np.pad(
            singular_values,
            (0, keep_modes - singular_values.size),
            mode="constant",
        )
    else:
        singular_values = singular_values[:keep_modes]

    total_energy = float(np.sum(singular_values * singular_values))
    return np.concatenate(
        [
            np.array([float(iteration), total_energy], dtype=np.float64),
            singular_values,
        ]
    )


def write_singular_values(
    comm: MPI.Comm,
    pod: POD,
    keep_modes: int,
    iteration: int,
    rotate_existing: bool,
) -> None:
    if comm.Get_rank() != 0:
        return

    path = "pod_singular_values.csv"
    if rotate_existing:
        rotated = rotate_existing_file(path)
        if rotated:
            log(comm, f"archived {path} -> {rotated}")

    mode = "w" if rotate_existing else "a"
    row = singular_value_row(pod, keep_modes, iteration)
    fmt = ["%d", "%.18e"] + ["%.18e"] * keep_modes

    with open(path, mode, encoding="utf-8") as handle:
        if rotate_existing:
            handle.write(singular_value_header(keep_modes) + "\n")
        np.savetxt(
            handle,
            row.reshape(1, -1),
            delimiter=",",
            fmt=fmt,
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


def mode_output_base(cfg: PODConfig) -> tuple[str, str]:
    file_name = os.path.normpath(cfg.mode_output_file_name)
    out_dir = os.path.dirname(file_name)
    out_stem = os.path.splitext(os.path.basename(file_name))[0]
    return out_dir, out_stem


def fld_mode_output_paths(cfg: PODConfig) -> tuple[str, str, str]:
    out_dir, out_stem = mode_output_base(cfg)
    sample_base = os.path.join(out_dir, f"{out_stem}0")
    meta_path = os.path.join(out_dir, f"{out_stem}0.nek5000")
    return out_dir, sample_base, meta_path


def ensure_mode_output_dir(comm: MPI.Comm, out_dir: str) -> None:
    if comm.Get_rank() == 0 and out_dir:
        os.makedirs(out_dir, exist_ok=True)
    comm.Barrier()


def fld_mode_field_names(n_total: int) -> list[str]:
    if n_total < 1 or n_total > 99:
        raise ValueError(
            f"Unsupported number of POD output fields: {n_total}"
        )

    if n_total == 1:
        return ["p"]
    if n_total == 2:
        return ["p", "t"]
    if n_total == 3:
        return ["u", "v", "w"]
    if n_total == 4:
        return ["p", "u", "v", "w"]

    names = ["p", "u", "v", "w", "t"]
    names.extend(f"s{i}" for i in range(n_total - 5))
    return names


def build_mode_output_data(
    ioh: IoHelp,
    pod: POD,
    cfg: PODConfig,
    zero_field: np.ndarray,
    n_avail: int,
    field_names: list[str],
) -> dict[str, np.ndarray]:
    flat_fields = []
    for mode_idx in range(cfg.keep_modes):
        flat_fields.extend(
            mode_fields_1d(ioh, pod, cfg, zero_field, n_avail, mode_idx)
        )

    return {
        name: np.asarray(field, dtype=cfg.mode_output_dtype)
        for name, field in zip(field_names, flat_fields)
    }


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


def vtk_mode_output_path(
    cfg: PODConfig,
    output_index: int,
) -> tuple[str, str, str]:
    out_dir, out_stem = mode_output_base(cfg)
    vtk_path = os.path.join(out_dir, f"{out_stem}_{output_index:05d}.vtkhdf")
    mesh_link = f"{out_stem}_{0:05d}.vtkhdf"
    return out_dir, vtk_path, mesh_link


def build_vtk_mode_output_data(
    ioh: IoHelp,
    pod: POD,
    cfg: PODConfig,
    zero_field: np.ndarray,
    n_avail: int,
) -> dict[str, np.ndarray]:
    return build_mode_output_data(
        ioh,
        pod,
        cfg,
        zero_field,
        n_avail,
        vtk_mode_field_names(cfg),
    )


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
    out_dir, vtk_path, mesh_link = vtk_mode_output_path(
        cfg, output_index
    )
    ensure_mode_output_dir(comm, out_dir)
    mode_data = build_vtk_mode_output_data(
        ioh,
        pod,
        cfg,
        zero_field,
        n_avail,
    )

    vtk_file = VTKHDFFile(
        comm,
        vtk_path,
        "w",
        parallel=comm.Get_size() > 1,
    )

    if output_index == 0:
        vtk_file.write_mesh_data(mesh.x, mesh.y, mesh.z)
    else:
        vtk_file.link_to_existing_mesh(mesh_link)

    for name, field in mode_data.items():
        vtk_file.write_point_data(
            name,
            field.astype(cfg.mode_output_dtype),
        )

    vtk_file.close()
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
        out_dir, sample_base, meta_path = fld_mode_output_paths(cfg)
        ensure_mode_output_dir(comm, out_dir)

        mode_data = build_mode_output_data(
            ioh,
            pod,
            cfg,
            zero_field,
            n_avail,
            fld_mode_field_names(cfg.keep_modes * cfg.n_fields),
        )

        write_data(
            comm,
            fname=f"{sample_base}.f{output_index:05d}",
            data_dict=mode_data,
            parallel_io=False,
            dtype=cfg.mode_output_dtype,
            msh=[mesh.x, mesh.y, mesh.z],
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
