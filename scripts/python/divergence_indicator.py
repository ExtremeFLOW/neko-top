#!/usr/bin/env python3
"""Report a divergence indicator for Neko/Nek5000 field output.

The rugby_ball forward output currently stores fields as generic scalar slots:
``s0`` is pressure-like, while ``s1``, ``s2`` and ``s3`` are velocity-like.
Those are therefore the default velocity fields.
"""

from __future__ import annotations

import argparse
import glob
import os
import re
from pathlib import Path

import numpy as np
from mpi4py import MPI
from pysemtools.datatypes.coef import Coef
from pysemtools.datatypes.field import FieldRegistry
from pysemtools.datatypes.msh import Mesh
from pysemtools.io.ppymech.neksuite import pynekread


def _field_step(path: str) -> int:
    match = re.search(r"\.f(\d+)$", path)
    return int(match.group(1)) if match else -1


def _latest_field(log_dir: Path, prefix: str) -> Path:
    candidates = glob.glob(str(log_dir / f"{prefix}.f[0-9][0-9][0-9][0-9][0-9]"))
    if not candidates:
        raise FileNotFoundError(f"No field files found for {log_dir / (prefix + '.fXXXXX')}")
    return Path(max(candidates, key=_field_step))


def _all_fields(log_dir: Path, prefix: str) -> list[Path]:
    candidates = glob.glob(str(log_dir / f"{prefix}.f[0-9][0-9][0-9][0-9][0-9]"))
    if not candidates:
        raise FileNotFoundError(f"No field files found for {log_dir / (prefix + '.fXXXXX')}")
    return [Path(p) for p in sorted(candidates, key=_field_step)]


def _global_sum(comm: MPI.Comm, value: float) -> float:
    return float(comm.allreduce(float(value), op=MPI.SUM))


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Compute SEM divergence indicators from pySEMTools/Neko field files."
    )
    parser.add_argument(
        "--log-dir",
        default="logs/rugby_ball",
        type=Path,
        help="Directory containing forward_fields output.",
    )
    parser.add_argument(
        "--prefix",
        default="forward_fields0",
        help="Field output prefix, without .fXXXXX.",
    )
    parser.add_argument(
        "--field",
        type=Path,
        default=None,
        help="Specific field file to analyze. If omitted, all prefix.fXXXXX files are processed.",
    )
    parser.add_argument(
        "--mesh-field",
        type=Path,
        default=None,
        help="Field file containing coordinates. Defaults to prefix.f00000 in log-dir.",
    )
    parser.add_argument(
        "--velocity-fields",
        nargs=3,
        default=("s1", "s2", "s3"),
        metavar=("U", "V", "W"),
        help="FieldRegistry names for velocity components.",
    )
    parser.add_argument(
        "--dtype",
        choices=("single", "double"),
        default="single",
        help="Data precision to request from pynekread.",
    )
    parser.add_argument(
        "--output",
        type=Path,
        default=None,
        help="Plot output path. Defaults to log-dir/divergence_vs_iteration.png.",
    )
    parser.add_argument(
        "--csv",
        type=Path,
        default=None,
        help="CSV output path. Defaults to log-dir/divergence_vs_iteration.csv.",
    )
    return parser.parse_args()


def compute_divergence_metrics(
    field_file: Path,
    comm: MPI.Comm,
    data_dtype: np.dtype,
    coef: Coef,
    velocity_fields: tuple[str, str, str],
) -> dict[str, float]:
    fld = FieldRegistry(comm)
    pynekread(str(field_file), comm, data_dtype=data_dtype, fld=fld, overwrite_fld=True)

    missing = [name for name in velocity_fields if name not in fld.registry]
    if missing:
        available = ", ".join(sorted(fld.registry.keys()))
        raise KeyError(
            f"Missing velocity field(s): {missing}. Available fields are: {available}"
        )

    u = fld.registry[velocity_fields[0]]
    v = fld.registry[velocity_fields[1]]
    w = fld.registry[velocity_fields[2]]

    dudx = coef.dudxyz(u, coef.drdx, coef.dsdx, coef.dtdx)
    dvdy = coef.dudxyz(v, coef.drdy, coef.dsdy, coef.dtdy)
    dwdz = coef.dudxyz(w, coef.drdz, coef.dsdz, coef.dtdz)
    div = dudx + dvdy + dwdz

    weights = coef.B.astype(div.dtype, copy=False)
    volume = _global_sum(comm, np.sum(weights))
    div_l2_sq = _global_sum(comm, np.sum(weights * div * div))
    vel_l2_sq = _global_sum(comm, np.sum(weights * (u * u + v * v + w * w)))
    div_abs_int = _global_sum(comm, np.sum(weights * np.abs(div)))
    div_int = _global_sum(comm, np.sum(weights * div))
    local_linf = float(np.max(np.abs(div))) if div.size else 0.0
    div_linf = float(comm.allreduce(local_linf, op=MPI.MAX))

    div_l2 = np.sqrt(div_l2_sq)
    div_rms = np.sqrt(div_l2_sq / volume)
    vel_l2 = np.sqrt(vel_l2_sq)
    rel_div_to_vel = div_l2 / vel_l2 if vel_l2 > 0 else np.nan

    return {
        "file_step": _field_step(str(field_file)),
        "volume": volume,
        "integral_div": div_int,
        "mean_abs_div": div_abs_int / volume,
        "rms_div": div_rms,
        "l2_div": div_l2,
        "linf_abs_div": div_linf,
        "l2_velocity": vel_l2,
        "l2_div_over_l2_velocity": rel_div_to_vel,
    }


def write_csv(path: Path, rows: list[dict[str, float]]) -> None:
    headers = [
        "file_step",
        "volume",
        "integral_div",
        "mean_abs_div",
        "rms_div",
        "l2_div",
        "linf_abs_div",
        "l2_velocity",
        "l2_div_over_l2_velocity",
    ]
    with path.open("w", encoding="utf-8") as handle:
        handle.write(",".join(headers) + "\n")
        for row in rows:
            handle.write(",".join(f"{row[key]:.16e}" for key in headers) + "\n")


def plot_series(path: Path, rows: list[dict[str, float]]) -> None:
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    steps = np.asarray([row["file_step"] for row in rows])
    rms = np.asarray([row["rms_div"] for row in rows])
    rel = np.asarray([row["l2_div_over_l2_velocity"] for row in rows])
    linf = np.asarray([row["linf_abs_div"] for row in rows])

    fig, ax = plt.subplots(figsize=(7.0, 4.2), constrained_layout=True)
    ax.plot(steps, rms, marker="o", markersize=3.5, linewidth=1.4, label="RMS(div)")
    ax.plot(steps, rel, marker="s", markersize=3.0, linewidth=1.2, label="L2(div)/L2(velocity)")
    ax.plot(steps, linf, marker="^", markersize=3.0, linewidth=1.0, label="Linf(abs(div))")
    ax.set_xlabel("Field iteration")
    ax.set_ylabel("Divergence indicator")
    ax.set_yscale("log")
    ax.grid(True, which="both", alpha=0.3)
    ax.legend()
    fig.savefig(path, dpi=200)
    plt.close(fig)


def main() -> None:
    args = parse_args()
    comm = MPI.COMM_WORLD
    rank = comm.Get_rank()

    data_dtype = np.single if args.dtype == "single" else np.double
    target_files = [args.field] if args.field is not None else _all_fields(args.log_dir, args.prefix)
    mesh_file = args.mesh_field if args.mesh_field is not None else args.log_dir / f"{args.prefix}.f00000"
    output = args.output if args.output is not None else args.log_dir / "divergence_vs_iteration.png"
    csv = args.csv if args.csv is not None else args.log_dir / "divergence_vs_iteration.csv"

    if not mesh_file.exists():
        raise FileNotFoundError(f"Mesh field file does not exist: {mesh_file}")
    for target_file in target_files:
        if not target_file.exists():
            raise FileNotFoundError(f"Target field file does not exist: {target_file}")

    msh = Mesh(comm, create_connectivity=False)
    pynekread(str(mesh_file), comm, data_dtype=data_dtype, msh=msh)
    coef = Coef(msh, comm)

    rows = [
        compute_divergence_metrics(
            target_file,
            comm,
            data_dtype,
            coef,
            tuple(args.velocity_fields),
        )
        for target_file in target_files
    ]

    if rank == 0:
        latest = rows[-1]
        print("Divergence indicator")
        print(f"  mesh field:        {mesh_file}")
        if len(target_files) == 1:
            print(f"  analyzed field:    {target_files[0]}")
        else:
            print(f"  analyzed fields:   {target_files[0]} ... {target_files[-1]}")
        print(f"  velocity fields:   {', '.join(args.velocity_fields)}")
        print(f"  samples:           {len(rows)}")
        print(f"  latest iteration:  {latest['file_step']:.0f}")
        print(f"  volume:            {latest['volume']:.8e}")
        print(f"  integral(div):     {latest['integral_div']:.8e}")
        print(f"  mean(abs(div)):    {latest['mean_abs_div']:.8e}")
        print(f"  rms(div):          {latest['rms_div']:.8e}")
        print(f"  L2(div):           {latest['l2_div']:.8e}")
        print(f"  Linf(abs(div)):    {latest['linf_abs_div']:.8e}")
        print(f"  L2(velocity):      {latest['l2_velocity']:.8e}")
        print(f"  L2(div)/L2(vel):   {latest['l2_div_over_l2_velocity']:.8e}")
        write_csv(csv, rows)
        print(f"  wrote CSV:         {csv}")
        if len(rows) > 1:
            plot_series(output, rows)
            print(f"  wrote plot:        {output}")


if __name__ == "__main__":
    main()
