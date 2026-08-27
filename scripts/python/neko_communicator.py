#!/usr/bin/env python3

import os
from typing import Optional

import numpy as np
from mpi4py import MPI


# Control protocol:
# - mode describes the coarse solver stage, while phase describes the
#   synchronization point inside that stage.
# - MODE_IDLE + PHASE_INIT is a one-off startup tick used to prove that the
#   control channel is alive. The Python driver may ignore it.
# - MODE_FORWARD + PHASE_FWD_RUNNING means Neko is in the forward solve and is
#   currently streaming snapshots to Python.
# - MODE_FORWARD + PHASE_FWD_DONE means the forward window is complete. Neko is
#   now blocked at the forward-to-adjoint boundary while Python finalizes the
#   POD basis, writes outputs, and prepares the adjoint data.
# - MODE_ADJOINT + PHASE_ADJ_RUNNING is Python's reply that POD modes and time
#   coefficients are ready, so Neko may enter adjoint restore/reconstruction.
# - MODE_ADJOINT + PHASE_ADJ_DONE means the adjoint window is complete and the
#   Python side should discard the old POD state and wait for a new forward
#   window.
# - MODE_STOP is a terminal shutdown message. The current implementation pairs
#   it with PHASE_ADJ_DONE as a placeholder phase, but the phase is not
#   semantically important once MODE_STOP is seen.
MODE_IDLE = 0
MODE_FORWARD = 1
MODE_ADJOINT = 2
MODE_STOP = 9

PHASE_INIT = 0
PHASE_FWD_RUNNING = 10
PHASE_FWD_DONE = 11
PHASE_ADJ_RUNNING = 20
PHASE_ADJ_DONE = 21

CTRL_TAG_STATE_INT = 4101
CTRL_TAG_STATE_REAL = 4102
CTRL_TAG_CMD = 4103


def get_comm_color() -> Optional[int]:
    value = os.getenv("NEKO_COMM_ID")
    if value is None:
        return None
    try:
        return int(value)
    except ValueError as exc:
        raise ValueError("Invalid NEKO_COMM_ID") from exc


def get_peer_root() -> int:
    value = os.getenv("NEKO_CTRL_PEER_ROOT")
    if value is None:
        raise ValueError("NEKO_CTRL_PEER_ROOT must be set for POD control.")
    try:
        peer_root = int(value)
    except ValueError as exc:
        raise ValueError("Invalid NEKO_CTRL_PEER_ROOT") from exc
    if peer_root < 0:
        raise ValueError("NEKO_CTRL_PEER_ROOT must be non-negative.")
    return peer_root


def make_local_comm(world: MPI.Comm) -> MPI.Comm:
    """Create this application's communicator in an ADIOS2-enabled MPMD run.

    This matches Neko's first MPI_COMM_WORLD split when Neko is compiled with
    ADIOS2. A non-ADIOS2 Neko build duplicates MPI_COMM_WORLD before splitting,
    so its Python peer must make that matching Dup call first.
    """
    comm_color = get_comm_color()
    if comm_color is None:
        return world.Dup()
    return world.Split(comm_color, world.Get_rank())


class CtrlClient:
    def __init__(
        self,
        world: MPI.Comm,
        peer_root: int,
        debug: bool,
    ):
        self.debug = debug
        self.world = world
        self.peer_root = peer_root

    def read_state(
        self,
    ) -> tuple[Optional[int], Optional[int], Optional[int], Optional[float]]:
        state_i = np.zeros(3, dtype=np.int32)
        state_t = np.zeros(1, dtype=np.float64)

        self.world.Recv(
            [state_i, MPI.INT],
            source=self.peer_root,
            tag=CTRL_TAG_STATE_INT,
        )
        self.world.Recv(
            [state_t, MPI.DOUBLE],
            source=self.peer_root,
            tag=CTRL_TAG_STATE_REAL,
        )

        result = (
            int(state_i[0]),
            int(state_i[1]),
            int(state_i[2]),
            float(state_t[0]),
        )
        if self.debug:
            print(
                "recv_state mode="
                f"{result[0]} phase={result[1]} "
                f"step={result[2]} t={result[3]}",
                flush=True,
            )
        return result

    def send_cmd(self, mode_cmd: int, phase_cmd: int) -> None:
        if self.debug:
            print(
                f"send_cmd mode={mode_cmd} phase={phase_cmd}",
                flush=True,
            )
        cmd = np.asarray([mode_cmd, phase_cmd], dtype=np.int32)
        self.world.Send(
            [cmd, MPI.INT],
            dest=self.peer_root,
            tag=CTRL_TAG_CMD,
        )

    def close(self) -> None:
        return
