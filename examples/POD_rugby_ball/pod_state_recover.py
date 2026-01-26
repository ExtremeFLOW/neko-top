#!/usr/bin/env python3
import runpy
import sys
from pathlib import Path

shared = (
    Path(__file__).resolve().parents[2]
    / "sources/state_recovery/POD_state_recover/pod_state_recover.py"
)

if not shared.exists():
    raise SystemExit(f"Shared pod_state_recover.py not found: {shared}")

sys.argv[0] = str(shared)
runpy.run_path(str(shared), run_name="__main__")
