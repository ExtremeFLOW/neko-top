#!/usr/bin/env python3

"""Validate the Python packages required by the ADIOS2 POD MPMD runtime."""

import os
import sys
import traceback

# Validation is run before mpirun. Do not let importing mpi4py initialize or
# finalize a singleton MPI world in this preflight process.
os.environ.setdefault("MPI4PY_RC_INITIALIZE", "0")
os.environ.setdefault("MPI4PY_RC_FINALIZE", "0")

try:
    import mpi4py

    mpi4py.rc.initialize = False
    mpi4py.rc.finalize = False
except Exception:
    pass


CHECKS = (
    ("numpy", "import numpy"),
    ("mpi4py.MPI", "from mpi4py import MPI"),
    ("adios2.bindings", "import adios2.bindings"),
    ("pysemtools.datatypes.coef", "from pysemtools.datatypes.coef import Coef"),
    ("pysemtools.datatypes.msh", "from pysemtools.datatypes.msh import Mesh"),
    (
        "pysemtools.io.adios2.stream",
        "from pysemtools.io.adios2.stream import DataStreamer",
    ),
    (
        "pysemtools.io.utils",
        "from pysemtools.io.utils import get_fld_from_ndarray",
    ),
    ("pysemtools.rom.io_help", "from pysemtools.rom.io_help import IoHelp"),
    ("pysemtools.rom.pod", "from pysemtools.rom.pod import POD"),
)


def main() -> int:
    namespace: dict[str, object] = {}
    failed = False

    for label, statement in CHECKS:
        try:
            exec(statement, namespace, namespace)
        except Exception:
            print(f"Failed to import {label}", file=sys.stderr)
            traceback.print_exc()
            failed = True

    return 1 if failed else 0


if __name__ == "__main__":
    raise SystemExit(main())
