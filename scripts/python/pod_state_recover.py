#!/usr/bin/env python3

import json
import os
import re
import sys
import time
from dataclasses import dataclass
from typing import Optional

import numpy as np
from mpi4py import MPI

from pysemtools.datatypes.coef import Coef
from pysemtools.datatypes.msh import Mesh
from pysemtools.io.adios2.stream import DataStreamer
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

    i_stream = int(state_recovery["i_stream"])
    timestep = float(case["case"]["time"]["timestep"])

    cfg = PODConfig(
        batch_size=int(state_recovery["batch_size"]),
        keep_modes=int(state_recovery["n_modes"]),
        i_stream=i_stream,
        dtype=dtype,
        write_modes=as_bool(state_recovery.get("write_modes")),
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
            log0(f"archived pod_time_coeffs.csv -> {rotated}")

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
        if idx < n_avail:
            fields_1d = ioh.split_narray_to_1dfields(pod.u_1t[:, idx])
            for field in fields_1d:
                ds.stream(field.astype(cfg.dtype, copy=False))
            continue

        for _ in range(cfg.n_fields):
            ds.stream(zero_field)


def init_runtime(
    comm: MPI.Comm,
    cfg: PODConfig,
) -> tuple[DataStreamer, np.ndarray, POD, IoHelp, list[np.ndarray]]:
    ds = DataStreamer(comm)

    x = recv_field(ds, cfg.dtype)
    y = recv_field(ds, cfg.dtype)
    z = recv_field(ds, cfg.dtype)

    initial_fields = recv_fields(ds, cfg.n_fields, cfg.dtype)

    mesh = Mesh(comm, x=x, y=y, z=z, create_connectivity=False)
    coef = Coef(mesh, comm)
    pod, ioh = make_pod(comm, coef.B, cfg)
    return ds, coef.B, pod, ioh, initial_fields


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

    ds, bm, pod, ioh, initial_fields = init_runtime(comm, cfg)
    times = []
    energy = EnergyState()
    add_snapshot(comm, pod, ioh, bm, initial_fields, 0.0, times, energy)

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
