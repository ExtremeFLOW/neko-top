#!/usr/bin/env python3
# pod_state_recover.py
#
# Debug-first in-situ driver for Neko <-> pySEMTools using ADIOS2 SST.
# - Uses ONE pySEMTools DataStreamer for the whole lifetime (mesh + forward snapshots + mode return).
# - Uses a rank-0-only control client (COMM_SELF) that:
#     reads  "neko_ctrl_state"
#     writes "neko_ctrl_cmd"
# - Adds very verbose, rank-tagged prints around every blocking ADIOS/MPI operation.
#
# IMPORTANT: This file is written for `import adios2.bindings as adios2`.
# In that binding, DefineVariable expects a numpy ndarray (including 0-D), not np.int32 scalars.

import os
import sys
import time
import json
import re
from typing import Optional, Tuple

import numpy as np
from mpi4py import MPI
import adios2.bindings as adios2

from pysemtools.io.adios2.stream import DataStreamer
from pysemtools.io.utils import get_fld_from_ndarray
from pysemtools.datatypes.msh import Mesh
from pysemtools.datatypes.coef import Coef
from pysemtools.rom.pod import POD
from pysemtools.rom.io_help import IoHelp


# -------------------------
# Control enums (keep consistent with Fortran/C++)
# -------------------------
MODE_IDLE    = 0
MODE_FORWARD = 1
MODE_ADJOINT = 2
MODE_STOP    = 9

PHASE_INIT        = 0
PHASE_FWD_RUNNING = 10
PHASE_FWD_DONE    = 11
PHASE_ADJ_RUNNING = 20
PHASE_ADJ_DONE    = 21


# -------------------------
# Logging helpers
# -------------------------
DEBUG = False

def log(comm: MPI.Comm, msg: str) -> None:
    if not DEBUG:
        return
    r = comm.Get_rank()
    s = comm.Get_size()
    print(f"[py r={r}/{s}] {msg}", flush=True)


def log0(msg: str) -> None:
    if not DEBUG:
        return
    print(f"[py ctrl r=0/1] {msg}", flush=True)


def load_case_config(case_path: str) -> dict:
    with open(case_path, "r", encoding="utf-8") as handle:
        raw = handle.read()
    # Allow trailing commas in .case files
    cleaned = re.sub(r",\s*([}\]])", r"\1", raw)
    return json.loads(cleaned)


def rotate_time_coeffs(path: str) -> Optional[str]:
    if not os.path.exists(path):
        return None
    base, ext = os.path.splitext(path)
    idx = 1
    while True:
        candidate = f"{base}_{idx:03d}{ext}"
        if not os.path.exists(candidate):
            os.rename(path, candidate)
            return candidate
        idx += 1


def wait_for_sst_contact(stem: str, timeout: float = 120.0, poll: float = 0.05) -> str:
    """
    SST rendezvous is file-based. Depending on ADIOS2 version/config,
    the contact file may be created as 'stem' or 'stem.sst'.
    We wait until one exists in the current working directory.
    """
    candidates = [stem, stem + ".sst"]
    t0 = time.time()
    cwd = os.getcwd()
    while time.time() - t0 < timeout:
        for c in candidates:
            if os.path.exists(c):
                return os.path.join(cwd, c)
        time.sleep(poll)
    raise RuntimeError(
        f"Timed out waiting for SST contact file for '{stem}' in cwd={cwd}. "
        f"Checked {candidates} for {timeout} seconds."
    )


# -------------------------
# Control client (rank-0 only, COMM_SELF)
# -------------------------
class CtrlClient:
    """
    Rank-0-only control channel.
    Uses COMM_SELF so it does NOT depend on python MPI ranks matching neko ranks.

    Reads:
      - stream name: "neko_ctrl_state"   (Neko writes)
    Writes:
      - stream name: "neko_ctrl_cmd"     (Python writes)
    """

    def __init__(self, debug: bool = True):
        self.debug = debug
        self.comm = MPI.COMM_SELF
        self.rank = 0  # COMM_SELF

        if self.debug:
            log0(f"CtrlClient: cwd={os.getcwd()}")
            log0("CtrlClient: creating ADIOS(COMM_SELF)")
        self.ad = adios2.ADIOS(self.comm)

        # Reader IO
        if self.debug:
            log0("CtrlClient: DeclareIO read (py_ctrl_read) + SetEngine(SST)")
        self.io_r = self.ad.DeclareIO("py_ctrl_read")
        self.io_r.SetEngine("SST")

        # IMPORTANT: wait for contact file before Open(Read) to avoid startup race.
        if self.debug:
            log0("CtrlClient: waiting for SST contact file for 'neko_ctrl_state'")
        found = wait_for_sst_contact("neko_ctrl_state", timeout=120.0)
        if self.debug:
            log0(f"CtrlClient: found contact file: {found}")

        if self.debug:
            log0("CtrlClient: Open('neko_ctrl_state', Read) (this may still block briefly)")
        self.rd = self.io_r.Open("neko_ctrl_state", adios2.Mode.Read)

        # Writer IO
        if self.debug:
            log0("CtrlClient: DeclareIO write (py_ctrl_write) + SetEngine(SST)")
        self.io_w = self.ad.DeclareIO("py_ctrl_write")
        self.io_w.SetEngine("SST")

        # Define scalar variables as 0-D ndarrays (required by adios2.bindings)
        self._mode_cmd_buf = np.zeros((), dtype=np.int32)
        self._phase_cmd_buf = np.zeros((), dtype=np.int32)

        if self.debug:
            log0("CtrlClient: DefineVariable(mode_cmd, phase_cmd) as 0-D arrays")
        self.v_mode_cmd = self.io_w.DefineVariable("mode_cmd", self._mode_cmd_buf)
        self.v_phase_cmd = self.io_w.DefineVariable("phase_cmd", self._phase_cmd_buf)

        # Open writer
        if self.debug:
            log0("CtrlClient: Open('neko_ctrl_cmd', Write)")
        self.wr = self.io_w.Open("neko_ctrl_cmd", adios2.Mode.Write)

        if self.debug:
            log0("CtrlClient: ready")

    def read_state(self, poll_sleep: float = 0.001) -> Tuple[Optional[int], Optional[int], Optional[int], Optional[float]]:
        """
        Poll BeginStep until OK or EndOfStream.
        Returns (mode, phase, step, time) or (None, None, None, None) on EOS.
        """
        while True:
            if self.debug:
                log0("CtrlClient: BeginStep(Read state)")
            st = self.rd.BeginStep()
            if st == adios2.StepStatus.OK:
                v_mode = self.io_r.InquireVariable("mode")
                v_phase = self.io_r.InquireVariable("phase")
                v_step = self.io_r.InquireVariable("step")
                v_time = self.io_r.InquireVariable("time")

                mode_buf = np.zeros((), dtype=np.int32)
                phase_buf = np.zeros((), dtype=np.int32)
                step_buf = np.zeros((), dtype=np.int32)
                time_buf = np.zeros((), dtype=np.float64)

                self.rd.Get(v_mode, mode_buf)
                self.rd.Get(v_phase, phase_buf)
                self.rd.Get(v_step, step_buf)
                self.rd.Get(v_time, time_buf)
                self.rd.EndStep()

                m = int(mode_buf)
                p = int(phase_buf)
                s = int(step_buf)
                t = float(time_buf)
                if self.debug:
                    log0(f"CtrlClient: got state mode={m} phase={p} step={s} t={t}")
                return m, p, s, t

            if st == adios2.StepStatus.NotReady:
                time.sleep(poll_sleep)
                continue

            # EndOfStream or other status
            if self.debug:
                log0(f"CtrlClient: state stream ended / not OK (status={st})")
            return None, None, None, None

    def send_cmd(self, mode_cmd: int, phase_cmd: int) -> None:
        if self.debug:
            log0(f"CtrlClient: send_cmd mode={mode_cmd} phase={phase_cmd}")
        self._mode_cmd_buf[...] = np.int32(mode_cmd)
        self._phase_cmd_buf[...] = np.int32(phase_cmd)

        if self.debug:
            log0("CtrlClient: BeginStep(Write cmd)")
        self.wr.BeginStep()
        self.wr.Put(self.v_mode_cmd, self._mode_cmd_buf)
        self.wr.Put(self.v_phase_cmd, self._phase_cmd_buf)
        self.wr.EndStep()
        if self.debug:
            log0("CtrlClient: EndStep(Write cmd)")

    def close(self) -> None:
        if self.debug:
            log0("CtrlClient: Close reader/writer")
        # Always close to avoid SST warnings
        self.rd.Close()
        self.wr.Close()


# -------------------------
# POD helper
# -------------------------
def make_pod(comm: MPI.Comm, bm: np.ndarray, n_fields: int, batch_size: int, keep_modes: int, dtype) -> Tuple[POD, IoHelp]:
    pod = POD(comm, number_of_modes_to_update=keep_modes, global_updates=True, auto_expand=False, bckend="numpy")
    ioh = IoHelp(comm,
                 number_of_fields=n_fields,
                 batch_size=batch_size,
                 field_size=bm.size,
                 field_data_type=dtype,
                 mass_matrix_data_type=dtype)
    mass_list = [np.copy(np.sqrt(bm)) for _ in range(n_fields)]
    ioh.copy_fieldlist_to_xi(mass_list)
    ioh.bm1sqrt[:, :] = np.copy(ioh.xi[:, :])
    return pod, ioh


# -------------------------
# Main
# -------------------------
def main() -> None:
    global DEBUG
    world = MPI.COMM_WORLD
    # keep your original "split for MPMD" approach (col=1)
    comm = world.Split(1, world.Get_rank())
    rank = comm.Get_rank()

    case_path = sys.argv[1] if len(sys.argv) > 1 else "POD_rugby_ball.case"
    case_path = os.path.abspath(case_path)
    case = load_case_config(case_path)

    sr = case.get("state_recovery", {})
    if not isinstance(sr, dict):
        raise KeyError("case['state_recovery'] must be a dict")

    debug = sr.get("debug", False)
    if isinstance(debug, str):
        debug = debug.strip().lower() in ("1", "true", "yes", "y", "t")
    else:
        debug = bool(debug)
    DEBUG = debug

    batch_size = int(sr["batch_size"])
    keep_modes = int(sr["n_modes"])
    i_stream = int(sr["i_stream"])
    dtype_str = str(sr["dtype"]).strip().lower()
    write_modes = sr.get("write_modes", False)
    if isinstance(write_modes, str):
        write_modes = write_modes.strip().lower() in ("1", "true", "yes", "y", "t")
    else:
        write_modes = bool(write_modes)

    # time settings (used for CSV)
    dt = float(case["case"]["time"]["timestep"])
    snapshot_dt = dt * i_stream

    if dtype_str == "single":
        dtype = np.float32
    elif dtype_str == "double":
        dtype = np.float64
    else:
        raise ValueError(f"Unsupported dtype '{dtype_str}'")

    log(comm, "Python: init DataStreamer (single lifetime)")
    ds = DataStreamer(comm)

    log(comm, "Python: receive mesh x,y,z")
    x = get_fld_from_ndarray(ds.recieve(), ds.lx, ds.ly, ds.lz, ds.nelv).astype(dtype)
    y = get_fld_from_ndarray(ds.recieve(), ds.lx, ds.ly, ds.lz, ds.nelv).astype(dtype)
    z = get_fld_from_ndarray(ds.recieve(), ds.lx, ds.ly, ds.lz, ds.nelv).astype(dtype)

    # Build mesh/coefs on all python ranks (pySEMTools expects this)
    msh = Mesh(comm, x=x, y=y, z=z, create_connectivity=False)
    coef = Coef(msh, comm)
    bm = coef.B

    n_fields = 3
    pod, ioh = make_pod(comm, bm, n_fields, batch_size, keep_modes, dtype)

    # Rank0 control client only (COMM_SELF)
    ctrl = CtrlClient(debug=DEBUG) if rank == 0 else None

    log(comm, "Python: Barrier after ctrl init")
    comm.Barrier()

    log(comm, "Python: entering control loop")

    # For padding when sending modes back
    # Use local field count if available, else fall back to x.size.
    my_count = getattr(ds, "py2f_field_my_count", None)
    if my_count is None:
        my_count = x.size
    zero_field = np.zeros(my_count, dtype=dtype)

    try:
        while True:
            # ---- read control state on rank0, broadcast to all python ranks ----
            if rank == 0:
                mode, phase, step, tcur = ctrl.read_state()
            else:
                mode, phase, step, tcur = None, None, None, None

            mode, phase, step, tcur = comm.bcast((mode, phase, step, tcur), root=0)
            log(comm, f"Loop: state mode={mode} phase={phase} step={step} t={tcur}")

            if mode is None:
                log(comm, "Loop: ctrl EOS -> exit")
                break
            if mode == MODE_STOP:
                log(comm, "Loop: MODE_STOP -> exit")
                break

            # ---- forward streaming snapshots ----
            if mode == MODE_FORWARD and phase == PHASE_FWD_RUNNING:
                log(comm, "Forward: about to recieve u,v,w (3x)")
                u = get_fld_from_ndarray(ds.recieve(), ds.lx, ds.ly, ds.lz, ds.nelv).astype(dtype)
                v = get_fld_from_ndarray(ds.recieve(), ds.lx, ds.ly, ds.lz, ds.nelv).astype(dtype)
                w = get_fld_from_ndarray(ds.recieve(), ds.lx, ds.ly, ds.lz, ds.nelv).astype(dtype)

                ioh.copy_fieldlist_to_xi([u, v, w])
                ioh.load_buffer(scale_snapshot=True)
                if ioh.update_from_buffer:
                    log(comm, f"POD.update(buff) with buffer_index={ioh.buffer_index}")
                    pod.update(comm, buff=ioh.buff[:, :ioh.buffer_index])
                continue

            # ---- forward finished: compute/write + send modes back + tell neko to switch ----
            if mode == MODE_FORWARD and phase == PHASE_FWD_DONE:
                log(comm, "Forward done: flush buffer -> update POD")
                if ioh.buffer_index > 0:
                    pod.update(comm, buff=ioh.buff[:, :ioh.buffer_index])

                log(comm, "Forward done: scale/rotate modes")
                pod.scale_modes(comm, bm1sqrt=ioh.bm1sqrt, op="div")
                pod.rotate_local_modes_to_global(comm)

                # Build time coeffs (rank0 writes)
                n_avail = min(keep_modes, pod.u_1t.shape[1]) if hasattr(pod, "u_1t") else 0
                if n_avail > 0:
                    A = (pod.d_1t[:n_avail, None] * pod.vt_1t[:n_avail, :]).T
                else:
                    A = np.zeros((0, 0), dtype=np.float64)

                nsnaps = A.shape[0]
                tvec = np.arange(nsnaps, dtype=np.float64) * snapshot_dt
                # pad to keep_modes columns
                if A.shape[1] < keep_modes:
                    A = np.hstack([A, np.zeros((nsnaps, keep_modes - A.shape[1]), dtype=A.dtype)])
                out = np.column_stack([tvec, A])
                header = "t," + ",".join([f"a{i+1}" for i in range(keep_modes)])

                if rank == 0:
                    if write_modes:
                        rotated = rotate_time_coeffs("pod_time_coeffs.csv")
                        if rotated:
                            log(comm, f"Rank0: archived pod_time_coeffs.csv -> {rotated}")
                    log(comm, "Rank0: writing pod_time_coeffs.csv")
                    np.savetxt("pod_time_coeffs.csv", out, delimiter=",", header=header, comments="")

                log(comm, "Barrier after CSV write")
                comm.Barrier()

                # IMPORTANT ordering:
                # 1) send modes back (so Neko can receive immediately after ctrl_wait_cmd returns)
                # 2) then send command that tells Neko to switch to adjoint
                #
                # If your Fortran side waits for cmd BEFORE receiving modes, flip this ordering.
                log(comm, f"Streaming modes back: expect 3*{keep_modes} streams")
                for j in range(keep_modes):
                    if j < n_avail:
                        fields1d = ioh.split_narray_to_1dfields(pod.u_1t[:, j])
                        log(comm, f"stream mode {j+1}/{keep_modes} (real)")
                        ds.stream(fields1d[0].astype(dtype, copy=False))
                        ds.stream(fields1d[1].astype(dtype, copy=False))
                        ds.stream(fields1d[2].astype(dtype, copy=False))
                    else:
                        log(comm, f"stream mode {j+1}/{keep_modes} (zero padded)")
                        ds.stream(zero_field); ds.stream(zero_field); ds.stream(zero_field)

                log(comm, "Barrier after streaming modes")
                comm.Barrier()

                if rank == 0:
                    log0("Rank0: sending ctrl cmd MODE_ADJOINT/PHASE_ADJ_RUNNING")
                    ctrl.send_cmd(MODE_ADJOINT, PHASE_ADJ_RUNNING)

                log(comm, "Barrier after sending ctrl cmd")
                comm.Barrier()
                continue

            # ---- adjoint done: reset POD buffers ----
            if mode == MODE_ADJOINT and phase == PHASE_ADJ_DONE:
                log(comm, "Adjoint done: reset POD/ioh")
                pod, ioh = make_pod(comm, bm, n_fields, batch_size, keep_modes, dtype)
                continue

            # if we get here, we're in an unhandled state
            log(comm, f"Unhandled state: mode={mode} phase={phase} -> idle spin")
            time.sleep(0.001)

    finally:
        # Clean shutdown
        if rank == 0 and ctrl is not None:
            ctrl.close()
        log(comm, "Barrier before ds.finalize()")
        comm.Barrier()
        log(comm, "Finalize DataStreamer")
        ds.finalize()


if __name__ == "__main__":
    main()
