#!/usr/bin/env python3
"""Make first/last rugby_ball design and velocity panels.

The figure has three rows and two columns:
  1. z-slice of the selected design scalar, black/white
  2. z-slice of velocity magnitude with velocity vectors and a sample line
  3. velocity magnitude along that physical sample line

The 2D panels are rotated with physical x pointing upward.  In plot
coordinates this is

    x_plot = -y_physical, y_plot = x_physical.
"""

from __future__ import annotations

import argparse
import glob
import re
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from mpi4py import MPI
from pysemtools.datatypes.field import FieldRegistry
from pysemtools.datatypes.msh import Mesh
from pysemtools.interpolation.probes import Probes
from pysemtools.io.ppymech.neksuite import pynekread


def field_step(path: str | Path) -> int:
    match = re.search(r"\.f(\d+)$", str(path))
    return int(match.group(1)) if match else -1


def fields_for_prefix(log_dir: Path, prefix: str) -> dict[int, Path]:
    files = glob.glob(str(log_dir / f"{prefix}.f[0-9][0-9][0-9][0-9][0-9]"))
    return {field_step(path): Path(path) for path in files}


def read_mesh_and_fields(path: Path, comm: MPI.Comm, dtype: np.dtype) -> tuple[Mesh, FieldRegistry]:
    msh = Mesh(comm, create_connectivity=False)
    fld = FieldRegistry(comm)
    pynekread(str(path), comm, data_dtype=dtype, msh=msh, fld=fld, overwrite_fld=True)
    return msh, fld


def read_fields(path: Path, comm: MPI.Comm, dtype: np.dtype) -> FieldRegistry:
    fld = FieldRegistry(comm)
    pynekread(str(path), comm, data_dtype=dtype, fld=fld, overwrite_fld=True)
    return fld


def slice_to_grid(
    msh: Mesh,
    values: np.ndarray,
    z_value: float,
    z_tol: float,
    decimals: int = 10,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    mask = np.abs(msh.z - z_value) <= z_tol
    if not np.any(mask):
        z_flat = msh.z.ravel()
        nearest = float(z_flat[np.argmin(np.abs(z_flat - z_value))])
        mask = np.isclose(msh.z, nearest)

    x = np.round(msh.x[mask].astype(np.float64), decimals)
    y = np.round(msh.y[mask].astype(np.float64), decimals)
    v = values[mask].astype(np.float64)

    xs = np.unique(x)
    ys = np.unique(y)
    xi = np.searchsorted(xs, x)
    yi = np.searchsorted(ys, y)

    accum = np.zeros((ys.size, xs.size), dtype=np.float64)
    count = np.zeros_like(accum)
    np.add.at(accum, (yi, xi), v)
    np.add.at(count, (yi, xi), 1.0)

    grid = np.divide(accum, count, out=np.full_like(accum, np.nan), where=count > 0)
    x_grid, y_grid = np.meshgrid(xs, ys)
    return x_grid, y_grid, grid


def interpolate_line_with_probes(
    msh: Mesh,
    fields: list[np.ndarray],
    comm: MPI.Comm,
    n_points: int,
    line_x: float,
    line_y0: float,
    line_y1: float,
    line_z: float,
) -> tuple[np.ndarray, list[np.ndarray]]:
    rank = comm.Get_rank()
    y = np.linspace(line_y0, line_y1, n_points)
    if rank == 0:
        probes_xyz = np.column_stack(
            (
                np.full(n_points, line_x),
                y,
                np.full(n_points, line_z),
            )
        )
    else:
        probes_xyz = None

    probes = Probes(
        comm,
        output_fname="/tmp/rugby_ball_line_probes.csv",
        probes=probes_xyz,
        msh=msh,
        write_coords=False,
        point_interpolator_type="multiple_point_legendre_numpy",
        max_pts=max(n_points, 128),
        find_points_comm_pattern="point_to_point",
    )
    probes.interpolate_from_field_list(0.0, fields, comm, write_data=False)

    if rank == 0:
        sampled = [probes.interpolated_fields[:, i + 1].copy() for i in range(len(fields))]
    else:
        sampled = [np.empty(0) for _ in fields]
    return y, sampled


def rotated_coords(x_grid: np.ndarray, y_grid: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    return -y_grid, x_grid


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Plot selected rugby_ball design and velocity panels.")
    parser.add_argument("--log-dir", type=Path, default=Path("logs/rugby_ball"))
    parser.add_argument("--design-prefix", default="design0")
    parser.add_argument("--forward-prefix", default="forward_fields0")
    parser.add_argument(
        "--design-field",
        default="s1",
        help="Design scalar to plot. pySEMTools uses zero-based names; use s0 for raw [0,1] design, s1 for second scalar.",
    )
    parser.add_argument("--velocity-fields", nargs=3, default=("s1", "s2", "s3"))
    parser.add_argument("--z", type=float, default=0.0)
    parser.add_argument("--z-tol", type=float, default=1.0e-7)
    parser.add_argument("--line-x", type=float, default=0.5)
    parser.add_argument("--line-y0", type=float, default=0.0)
    parser.add_argument("--line-y1", type=float, default=1.0)
    parser.add_argument("--line-z", type=float, default=0.0)
    parser.add_argument("--line-points", type=int, default=300)
    parser.add_argument("--quiver-stride", type=int, default=5)
    parser.add_argument("--dtype", choices=("single", "double"), default="single")
    parser.add_argument("--output", type=Path, default=Path("logs/rugby_ball/first_last_panels.png"))
    parser.add_argument(
        "--steps",
        nargs="+",
        type=int,
        default=None,
        help="Iterations to plot as columns, e.g. --steps 1 10 20. Defaults to first and last common iteration.",
    )
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    comm = MPI.COMM_WORLD
    rank = comm.Get_rank()
    dtype = np.single if args.dtype == "single" else np.double

    design_files = fields_for_prefix(args.log_dir, args.design_prefix)
    forward_files = fields_for_prefix(args.log_dir, args.forward_prefix)
    common_steps = sorted(set(design_files).intersection(forward_files))
    if not common_steps:
        raise FileNotFoundError("No common design/forward field iterations found.")
    if args.steps is None:
        steps = [common_steps[0], common_steps[-1]]
    else:
        missing = [step for step in args.steps if step not in common_steps]
        if missing:
            available = f"{common_steps[0]}..{common_steps[-1]}"
            raise ValueError(f"Requested iteration(s) not available: {missing}. Available range: {available}")
        steps = args.steps

    msh, _ = read_mesh_and_fields(design_files[common_steps[0]], comm, dtype)
    cases = []
    for step in steps:
        design = read_fields(design_files[step], comm, dtype)
        forward = read_fields(forward_files[step], comm, dtype)

        if args.design_field not in design.registry:
            raise KeyError(f"{args.design_field} not found in {design_files[step]}; available {list(design.registry)}")
        missing_vel = [name for name in args.velocity_fields if name not in forward.registry]
        if missing_vel:
            raise KeyError(f"{missing_vel} not found in {forward_files[step]}; available {list(forward.registry)}")

        u = forward.registry[args.velocity_fields[0]]
        v = forward.registry[args.velocity_fields[1]]
        w = forward.registry[args.velocity_fields[2]]
        speed = np.sqrt(u * u + v * v + w * w)

        xg, yg, design_grid = slice_to_grid(msh, design.registry[args.design_field], args.z, args.z_tol)
        _, _, speed_grid = slice_to_grid(msh, speed, args.z, args.z_tol)
        _, _, u_grid = slice_to_grid(msh, u, args.z, args.z_tol)
        _, _, v_grid = slice_to_grid(msh, v, args.z, args.z_tol)

        line_y, line_fields = interpolate_line_with_probes(
            msh,
            [u, v, w],
            comm,
            args.line_points,
            args.line_x,
            args.line_y0,
            args.line_y1,
            args.line_z,
        )
        if rank == 0:
            line_speed = np.sqrt(line_fields[0] ** 2 + line_fields[1] ** 2 + line_fields[2] ** 2)
        else:
            line_speed = np.empty(0)

        cases.append(
            {
                "step": step,
                "xg": xg,
                "yg": yg,
                "design": design_grid,
                "speed": speed_grid,
                "u": u_grid,
                "v": v_grid,
                "line_y": line_y,
                "line_speed": line_speed,
            }
        )

    if rank != 0:
        return

    xrot, yrot = rotated_coords(cases[0]["xg"], cases[0]["yg"])
    x_limits = (float(np.nanmin(xrot)), float(np.nanmax(xrot)))
    y_limits = (float(np.nanmin(yrot)), float(np.nanmax(yrot)))
    line_xrot = -np.linspace(args.line_y0, args.line_y1, args.line_points)
    line_yrot = np.full(args.line_points, args.line_x)

    design_min = min(float(np.nanmin(case["design"])) for case in cases)
    design_max = max(float(np.nanmax(case["design"])) for case in cases)
    speed_max = max(float(np.nanmax(case["speed"])) for case in cases)
    line_max = max(float(np.nanmax(case["line_speed"])) for case in cases)
    line_min = min(float(np.nanmin(case["line_speed"])) for case in cases)
    line_pad = 0.05 * max(line_max - line_min, line_max, 1.0)
    linear_ylim = (max(0.0, line_min - line_pad), line_max + line_pad)
    positive_line_min = min(
        float(np.nanmin(case["line_speed"][case["line_speed"] > 0.0]))
        for case in cases
        if np.any(case["line_speed"] > 0.0)
    )
    log_ylim = (positive_line_min * 0.5, line_max * 1.5)
    n_cols = len(cases)

    fig, axes = plt.subplots(
        4,
        n_cols,
        figsize=(5.25 * n_cols, 14.4),
        gridspec_kw={"height_ratios": [1.0, 1.0, 0.42, 0.42]},
        constrained_layout=True,
        squeeze=False,
    )

    for col, case in enumerate(cases):
        xr, yr = rotated_coords(case["xg"], case["yg"])
        title_suffix = f"iteration {case['step']}"

        ax = axes[0, col]
        im_design = ax.pcolormesh(
            xr,
            yr,
            case["design"],
            shading="auto",
            cmap="gray",
            vmin=design_min,
            vmax=design_max,
        )
        ax.set_title(f"Design {title_suffix}")
        ax.set_xlim(x_limits)
        ax.set_ylim(y_limits)
        ax.set_aspect("equal", adjustable="box")
        ax.set_ylabel("physical x")
        fig.colorbar(im_design, ax=ax, fraction=0.046, pad=0.02)

        ax = axes[1, col]
        im_speed = ax.pcolormesh(
            xr,
            yr,
            case["speed"],
            shading="auto",
            cmap="viridis",
            vmin=0.0,
            vmax=speed_max,
        )
        stride = max(args.quiver_stride, 1)
        qx = xr[::stride, ::stride]
        qy = yr[::stride, ::stride]
        qu = -case["v"][::stride, ::stride]
        qv = case["u"][::stride, ::stride]
        ax.quiver(qx, qy, qu, qv, color="white", angles="xy", scale_units="xy", scale=18.0, width=0.002)
        ax.plot(line_xrot, line_yrot, color="red", linewidth=1.5)
        ax.set_title(f"Velocity magnitude {title_suffix}")
        ax.set_xlim(x_limits)
        ax.set_ylim(y_limits)
        ax.set_aspect("equal", adjustable="box")
        ax.set_ylabel("physical x")
        fig.colorbar(im_speed, ax=ax, fraction=0.046, pad=0.02)

        ax = axes[2, col]
        ax.plot(-case["line_y"], case["line_speed"], color="black", linewidth=1.5)
        ax.set_xlim(x_limits)
        ax.set_ylim(linear_ylim)
        ax.grid(True, alpha=0.3)
        ax.set_xlabel("- physical y")
        ax.set_ylabel("|u|")
        ax.set_title(f"Line x={args.line_x:g}, z={args.line_z:g}")

        ax = axes[3, col]
        positive_speed = np.maximum(case["line_speed"], np.finfo(float).tiny)
        ax.plot(-case["line_y"], positive_speed, color="black", linewidth=1.5)
        ax.set_xlim(x_limits)
        ax.set_ylim(log_ylim)
        ax.set_yscale("log")
        ax.grid(True, which="both", alpha=0.3)
        ax.set_xlabel("- physical y")
        ax.set_ylabel("|u|")
        ax.set_title("Line magnitude, log scale")

    for ax in axes[:2, :].ravel():
        ax.set_xlabel("- physical y")

    args.output.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(args.output, dpi=220)
    plt.close(fig)
    print(f"Wrote {args.output}")


if __name__ == "__main__":
    main()
