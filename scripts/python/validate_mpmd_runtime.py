#!/usr/bin/env python3

"""Validate the Python packages required by the ADIOS2/POD MPMD runtime."""

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


def module_path(namespace: dict[str, object], name: str) -> str:
    module = namespace.get(name)
    if module is None:
        return "<not imported>"
    return str(getattr(module, "__file__", "<built-in>"))


def main() -> int:
    namespace: dict[str, object] = {}
    failed = False

    print(f"python: {sys.executable}")
    print(f"prefix: {sys.prefix}")

    for label, statement in CHECKS:
        try:
            exec(statement, namespace, namespace)
        except Exception:
            print(f"Failed to import {label}", file=sys.stderr)
            traceback.print_exc()
            failed = True

    if not failed:
        exec("import numpy, mpi4py, adios2, pysemtools", namespace, namespace)
        print(f"numpy: {module_path(namespace, 'numpy')}")
        print(f"mpi4py: {module_path(namespace, 'mpi4py')}")
        print(f"adios2: {module_path(namespace, 'adios2')}")
        print(f"pysemtools: {module_path(namespace, 'pysemtools')}")

    return 1 if failed else 0


if __name__ == "__main__":
    raise SystemExit(main())
