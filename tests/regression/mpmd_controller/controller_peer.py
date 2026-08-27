#!/usr/bin/env python3

from pathlib import Path
import sys

from mpi4py import MPI


REPO_ROOT = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(REPO_ROOT / "scripts" / "python"))

from neko_communicator import (  # noqa: E402
    MODE_ADJOINT,
    MODE_FORWARD,
    PHASE_ADJ_RUNNING,
    PHASE_FWD_DONE,
    CtrlClient,
    get_peer_root,
)


def main() -> int:
    world = MPI.COMM_WORLD
    if world.Get_size() != 2 or world.Get_rank() != 0:
        raise RuntimeError("Expected one Python rank at world rank 0.")

    # mpmd_launch_shared assigns Python color 1 and Neko color 0. The launcher
    # does not perform these collectives on behalf of either application.
    # Match Neko's non-ADIOS2 initialization: duplicate the world communicator
    # before splitting it into the Neko and Python color groups.
    neko_global_comm = world.Dup()
    neko_local_comm = neko_global_comm.Split(1, world.Get_rank())

    # Control traffic intentionally stays on world so it can cross from the
    # Neko group to this Python group. POD data uses local communicators later.
    client = CtrlClient(world, get_peer_root(), debug=False)
    mode, phase, step, time = client.read_state()
    if (mode, phase, step) != (MODE_FORWARD, PHASE_FWD_DONE, 7):
        raise RuntimeError(f"Unexpected Neko state: {(mode, phase, step)}")
    if abs(time - 1.25) > 1.0e-14:
        raise RuntimeError(f"Unexpected Neko time: {time}")

    client.send_cmd(MODE_ADJOINT, PHASE_ADJ_RUNNING)
    print("MPMD controller regression passed", flush=True)

    # Keep both communicators alive until Neko has completed its matching setup.
    _ = (neko_global_comm, neko_local_comm)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
