#!/usr/bin/env python3

import os
import sys
import time

import numpy as np
from mpi4py import MPI

from pysemtools.datatypes.coef import Coef
from pysemtools.datatypes.msh import Mesh
from pysemtools.io.adios2.stream import DataStreamer
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
from pod_state_recover_tools import EnergyState
from pod_state_recover_tools import PODConfig
from pod_state_recover_tools import add_snapshot
from pod_state_recover_tools import build_time_coefficients
from pod_state_recover_tools import flush_buffer
from pod_state_recover_tools import load_pod_config
from pod_state_recover_tools import make_pod
from pod_state_recover_tools import recv_field
from pod_state_recover_tools import recv_fields
from pod_state_recover_tools import report_energy_capture
from pod_state_recover_tools import set_debug
from pod_state_recover_tools import stream_mode_fields
from pod_state_recover_tools import write_modes_to_disk
from pod_state_recover_tools import write_singular_values
from pod_state_recover_tools import write_time_coefficients


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
    world = MPI.COMM_WORLD
    comm = make_local_comm(world)
    rank = comm.Get_rank()
    peer_root = get_peer_root()

    case_path = sys.argv[1] if len(sys.argv) > 1 else "POD_rugby_ball.case"
    case_path = os.path.abspath(case_path)
    _, cfg = load_pod_config(case_path)
    set_debug(cfg.debug)

    ds, mesh, bm, pod, ioh, initial_fields = init_runtime(comm, cfg)
    times = []
    energy = EnergyState()
    add_snapshot(comm, pod, ioh, bm, initial_fields, 0.0, times, energy)
    mode_output_index = 0
    pod_iteration = 0
    singular_value_history_new = True

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
                pod_iteration += 1
                write_singular_values(
                    comm,
                    pod,
                    cfg.keep_modes,
                    pod_iteration,
                    singular_value_history_new,
                )
                singular_value_history_new = False
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
