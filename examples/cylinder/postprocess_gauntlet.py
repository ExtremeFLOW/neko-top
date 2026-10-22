#!/usr/bin/env python3
"""Post-process the Re40, Re60, Re180, and Re200 cylinder gauntlet results.

This script must be run with a Python environment containing pySEMTools.  It
uses pySEMTools to read the Neko field files and extracts one spectral-element
plane in the homogeneous z direction for all field plots.

For the steady Re40/Re60 cases it produces:

* an accumulating baseflow/Newton/GMRES convergence PNG sequence;
* one overlaid direct/adjoint stability spectrum; and
* separate PNG-frame sequences of the leading direct and adjoint modes; and
* transient-growth singular-vector input/output pairs, when available.

For the periodic Re180/Re200 cases it produces:

* accumulating periodic-orbit Newton/GMRES and period-convergence PNG sequences;
* numbered spanwise-vorticity shedding frames; and
* overlaid direct-LNS/adjoint Floquet multiplier and exponent spectra.
"""

from __future__ import annotations

import argparse
import cmath
import csv
import json
import logging
import math
import re
from dataclasses import dataclass
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.tri as mtri
import numpy as np
from matplotlib.colors import Normalize, TwoSlopeNorm
from matplotlib.patches import Circle
from mpi4py import MPI
from pysemtools.datatypes.coef import Coef
from pysemtools.datatypes.field import FieldRegistry
from pysemtools.datatypes.msh import Mesh
from pysemtools.datatypes.msh_connectivity import MeshConnectivity
from pysemtools.io.ppymech.neksuite import pynekread


FLOAT = r"[+\-]?(?:\d+(?:\.\d*)?|\.\d+)(?:[Ee][+\-]?\d+)?"
NEWTON_RE = re.compile(
    rf"Start step\s+(?P<step>\d+):\s+rnorm=\s*(?P<res>{FLOAT}),"
    rf"\s*tol=\s*(?P<tol>{FLOAT})"
)
GMRES_RE = re.compile(
    rf"GMRES\(k\)\s+(?P<kind>init|inner|outer) step(?:\s+\d+)?\s*:"
    rf"\s*\|res\|=\s*(?P<res>{FLOAT}),\s*tol=\s*(?P<tol>{FLOAT})"
)
FIELD_INDEX_RE = re.compile(r"\.f(?P<index>\d+)$")


@dataclass
class NewtonStep:
    index: int
    residual: float
    tolerance: float
    gmres_residuals: list[float]


@dataclass
class SliceMesh:
    x: np.ndarray
    y: np.ndarray
    triangles: np.ndarray
    z: float


@dataclass
class EigenvalueRow:
    iteration: int
    value: complex
    residual: float
    converged: bool


@dataclass
class FloquetRow:
    index: int
    multiplier: complex
    exponent: complex
    residual: float


def parse_args() -> argparse.Namespace:
    script_dir = Path(__file__).resolve().parent
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--results",
        type=Path,
        default=script_dir / "gauntlet" / "results",
        help="Gauntlet results root (default: %(default)s)",
    )
    parser.add_argument(
        "--output",
        type=Path,
        default=script_dir / "gauntlet" / "postprocessing",
        help="Output directory (default: %(default)s)",
    )
    parser.add_argument(
        "--reynolds",
        nargs="+",
        choices=("Re40", "Re60", "Re180", "Re200"),
        default=("Re40", "Re60", "Re180", "Re200"),
        help="Reynolds-number cases to process",
    )
    parser.add_argument(
        "--component",
        choices=("u", "v", "w"),
        default="u",
        help="Field component rendered in convergence and mode plots",
    )
    parser.add_argument("--mode-frames", type=int, default=48)
    parser.add_argument("--xlim", nargs=2, type=float, default=(-2.0, 15.0))
    parser.add_argument("--ylim", nargs=2, type=float, default=(-5.0, 5.0))
    parser.add_argument("--skip-baseflow", action="store_true")
    parser.add_argument("--skip-spectra", action="store_true")
    parser.add_argument("--skip-modes", action="store_true")
    parser.add_argument("--skip-transient-growth", action="store_true")
    parser.add_argument("--skip-shedding", action="store_true")
    parser.add_argument(
        "--vorticity-limit",
        type=float,
        default=None,
        help="Symmetric shedding-frame colour limit (default: robust value from the data)",
    )
    return parser.parse_args()


def require_file(path: Path) -> Path:
    if not path.is_file():
        raise FileNotFoundError(f"Required gauntlet result is missing: {path}")
    return path


def field_sort_key(path: Path) -> int:
    match = FIELD_INDEX_RE.search(path.name)
    if match is None:
        raise ValueError(f"Cannot extract field index from {path.name}")
    return int(match.group("index"))


def parse_newton_log(path: Path) -> list[NewtonStep]:
    steps: list[NewtonStep] = []
    current: NewtonStep | None = None

    for line in path.read_text(encoding="utf-8", errors="replace").splitlines():
        newton_match = NEWTON_RE.search(line)
        if newton_match:
            current = NewtonStep(
                index=int(newton_match.group("step")),
                residual=float(newton_match.group("res")),
                tolerance=float(newton_match.group("tol")),
                gmres_residuals=[],
            )
            steps.append(current)
            continue

        gmres_match = GMRES_RE.search(line)
        if current is None or gmres_match is None:
            continue
        # The outer-step value repeats the last inner-step value.  Keeping only
        # init/inner entries gives a continuous Krylov-iteration index even
        # when restarted GMRES resets its printed local counter.
        if gmres_match.group("kind") != "outer":
            current.gmres_residuals.append(float(gmres_match.group("res")))

    if not steps:
        raise ValueError(f"No Newton convergence records found in {path}")
    if any(not step.gmres_residuals for step in steps):
        raise ValueError(f"At least one Newton step has no GMRES history in {path}")
    return steps


def build_slice_mesh(mesh: Mesh) -> tuple[SliceMesh, int]:
    if mesh.x.ndim != 4:
        raise ValueError(f"Expected SEM arrays with four dimensions, got {mesh.x.shape}")

    z_by_plane = np.mean(mesh.z, axis=(0, 2, 3))
    target_z = 0.5 * (float(np.min(mesh.z)) + float(np.max(mesh.z)))
    plane = int(np.argmin(np.abs(z_by_plane - target_z)))
    x = np.asarray(mesh.x[:, plane, :, :], dtype=np.float64).reshape(-1)
    y = np.asarray(mesh.y[:, plane, :, :], dtype=np.float64).reshape(-1)

    nel, ly, lx = mesh.x[:, plane, :, :].shape
    triangles = np.empty((2 * nel * (ly - 1) * (lx - 1), 3), dtype=np.int64)
    n = 0
    points_per_element = ly * lx
    for element in range(nel):
        base = element * points_per_element
        for j in range(ly - 1):
            for i in range(lx - 1):
                lower_left = base + j * lx + i
                lower_right = lower_left + 1
                upper_left = base + (j + 1) * lx + i
                upper_right = upper_left + 1
                triangles[n] = (lower_left, lower_right, upper_right)
                triangles[n + 1] = (lower_left, upper_right, upper_left)
                n += 2

    return SliceMesh(x=x, y=y, triangles=triangles, z=float(z_by_plane[plane])), plane


def read_slice_sequence(
    files: list[Path], component: str
) -> tuple[SliceMesh, list[np.ndarray], list[np.ndarray]]:
    if not files:
        raise ValueError("No field files were supplied")

    comm = MPI.COMM_SELF
    mesh = Mesh(comm)
    field = FieldRegistry(comm)
    geometry: SliceMesh | None = None
    plane = 0
    values: list[np.ndarray] = []
    signatures: list[np.ndarray] = []

    for file_index, path in enumerate(files):
        pynekread(
            str(path),
            comm,
            data_dtype=np.float64,
            msh=mesh if file_index == 0 else None,
            fld=field,
            overwrite_fld=True,
        )
        if file_index == 0:
            geometry, plane = build_slice_mesh(mesh)
        if component not in field.registry:
            raise KeyError(f"Field {component!r} is not present in {path}")

        values.append(
            np.asarray(field.registry[component][:, plane, :, :], dtype=np.float64)
            .reshape(-1)
            .copy()
        )
        velocity_parts = [
            np.asarray(field.registry[name][:, plane, :, :], dtype=np.float64).reshape(-1)
            for name in ("u", "v", "w")
            if name in field.registry
        ]
        signatures.append(np.concatenate(velocity_parts).copy())

    assert geometry is not None
    return geometry, values, signatures


def unique_state_indices(signatures: list[np.ndarray]) -> list[int]:
    indices: list[int] = []
    previous: np.ndarray | None = None
    for index, signature in enumerate(signatures):
        scale = max(float(np.linalg.norm(signature)), 1.0)
        if previous is None or float(np.linalg.norm(signature - previous)) > 1.0e-12 * scale:
            indices.append(index)
            previous = signature
    return indices


def add_cylinder(ax: plt.Axes) -> None:
    ax.add_patch(Circle((0.0, 0.0), 0.5, facecolor="white", edgecolor="black", lw=1.0, zorder=5))


def draw_slice(
    ax: plt.Axes,
    geometry: SliceMesh,
    values: np.ndarray,
    *,
    cmap: str,
    norm: Normalize,
    xlim: tuple[float, float],
    ylim: tuple[float, float],
):
    triangulation = mtri.Triangulation(geometry.x, geometry.y, geometry.triangles)
    image = ax.tripcolor(triangulation, values, shading="gouraud", cmap=cmap, norm=norm)
    add_cylinder(ax)
    ax.set_aspect("equal", adjustable="box")
    ax.set_xlim(xlim)
    ax.set_ylim(ylim)
    ax.set_xlabel(r"$x/D$")
    ax.set_ylabel(r"$y/D$")
    return image


def prepare_frames_directory(path: Path) -> Path:
    """Create a frame directory and remove PNGs left by an earlier render."""
    path.mkdir(parents=True, exist_ok=True)
    for frame in path.glob("*.png"):
        frame.unlink()
    return path


def make_baseflow_convergence(
    case: str,
    case_results: Path,
    output_dir: Path,
    component: str,
    xlim: tuple[float, float],
    ylim: tuple[float, float],
) -> None:
    baseflow_dir = case_results / "baseflow"
    steps = parse_newton_log(require_file(baseflow_dir / "lightkrylov.log"))
    outer_dir = baseflow_dir / "viz" / f"outer_baseflow_newton_{case}"
    fields = sorted(outer_dir.glob("*.f[0-9]*"), key=field_sort_key)
    if not fields:
        raise FileNotFoundError(f"No Newton outer-state fields found in {outer_dir}")

    geometry, raw_values, signatures = read_slice_sequence(fields, component)
    states = [raw_values[index] for index in unique_state_indices(signatures)]
    if len(states) < len(steps) + 1:
        raise ValueError(
            f"{case}: found {len(states)} distinct outer states for "
            f"{len(steps)} Newton iterations; expected an initial state plus every update"
        )
    # State zero is the initial guess.  The selected states are the updates at
    # the end of each Newton iteration; linspace also keeps this robust if a
    # solver version emits extra distinct line-search evaluations.
    state_indices = np.rint(np.linspace(1, len(states) - 1, len(steps))).astype(int)
    iteration_states = [states[index] for index in state_indices]

    vmin = min(float(np.nanmin(values)) for values in iteration_states)
    vmax = max(float(np.nanmax(values)) for values in iteration_states)
    norm = Normalize(vmin=vmin, vmax=vmax)
    colors = plt.get_cmap("turbo")(np.linspace(0.05, 0.95, len(steps)))
    all_positive = [step.residual for step in steps]
    all_positive.extend(value for step in steps for value in step.gmres_residuals)
    residual_min = min(value for value in all_positive if value > 0.0)
    residual_max = max(all_positive)
    max_gmres_iterations = max(len(step.gmres_residuals) for step in steps)

    frames_dir = prepare_frames_directory(output_dir / "baseflow_convergence_frames")
    for last_index, (current_step, current_state) in enumerate(
        zip(steps, iteration_states)
    ):
        visible_steps = steps[: last_index + 1]
        fig = plt.figure(figsize=(14, 9), constrained_layout=True)
        grid = fig.add_gridspec(2, 2, height_ratios=(1.18, 1.0))
        field_ax = fig.add_subplot(grid[0, :])
        newton_ax = fig.add_subplot(grid[1, 0])
        gmres_ax = fig.add_subplot(grid[1, 1])

        image = draw_slice(
            field_ax,
            geometry,
            current_state,
            cmap="viridis",
            norm=norm,
            xlim=xlim,
            ylim=ylim,
        )
        field_ax.set_title(
            f"{case} baseflow after Newton iteration {current_step.index} "
            f"({component}, z={geometry.z:.3g})"
        )
        fig.colorbar(image, ax=field_ax, pad=0.015, label=component)

        newton_ax.semilogy(
            [item.index for item in visible_steps],
            [item.residual for item in visible_steps],
            "o-",
            color="black",
            lw=1.8,
        )
        newton_ax.set_xlim(0.5, len(steps) + 0.5)
        newton_ax.set_ylim(residual_min * 0.5, residual_max * 2.0)
        newton_ax.set_xlabel("Newton iteration")
        newton_ax.set_ylabel(r"$\|F(q)\|$")
        newton_ax.set_title("Newton convergence")
        newton_ax.grid(True, which="both", alpha=0.3)

        for color, item in zip(colors[: last_index + 1], visible_steps):
            gmres_ax.semilogy(
                np.arange(len(item.gmres_residuals)),
                item.gmres_residuals,
                color=color,
                lw=1.3,
                label=f"Newton {item.index}",
            )
            gmres_ax.axhline(
                item.tolerance, color=color, ls="--", lw=0.7, alpha=0.45
            )
        gmres_ax.set_xlim(0, max_gmres_iterations - 1)
        gmres_ax.set_ylim(residual_min * 0.5, residual_max * 2.0)
        gmres_ax.set_xlabel("Inner Krylov iteration")
        gmres_ax.set_ylabel(r"$\|r\|$")
        gmres_ax.set_title("GMRES convergence by Newton iteration")
        gmres_ax.grid(True, which="both", alpha=0.3)
        gmres_ax.legend(ncol=2, fontsize=8, loc="best")

        fig.suptitle(f"Cylinder gauntlet: {case}", fontsize=15)
        output = frames_dir / f"baseflow_convergence_{last_index + 1:04d}.png"
        fig.savefig(output, dpi=200)
        plt.close(fig)
    print(f"[{case}] wrote {len(steps)} convergence frames to {frames_dir}")


def read_period_history(path: Path) -> tuple[np.ndarray, np.ndarray]:
    iterations: list[int] = []
    periods: list[float] = []
    with require_file(path).open(newline="", encoding="utf-8") as stream:
        for row in csv.DictReader(stream):
            try:
                iterations.append(int(row["iteration"]))
                periods.append(float(row["T"]))
            except (KeyError, TypeError, ValueError) as error:
                raise ValueError(f"Malformed period-history row in {path}: {row}") from error
    if not periods:
        raise ValueError(f"No period estimates found in {path}")
    return np.asarray(iterations), np.asarray(periods)


def make_periodic_orbit_convergence(
    case: str,
    case_results: Path,
    output_dir: Path,
    component: str,
    xlim: tuple[float, float],
    ylim: tuple[float, float],
) -> None:
    orbit_dir = case_results / "periodic_orbit"
    steps = parse_newton_log(require_file(orbit_dir / "lightkrylov.log"))
    outer_dir = orbit_dir / "viz" / f"outer_periodic_orbit_{case}"
    fields = sorted(outer_dir.glob("*.f[0-9]*"), key=field_sort_key)
    if not fields:
        raise FileNotFoundError(f"No periodic-orbit Newton states found in {outer_dir}")

    geometry, values, signatures = read_slice_sequence(fields, component)
    period_iterations, periods = read_period_history(
        orbit_dir / "viz" / f"outer_periodic_orbit_{case}_period.csv"
    )
    distinct_indices = unique_state_indices(signatures)
    if len(distinct_indices) < len(steps) + 1:
        raise ValueError(
            f"{case}: found {len(distinct_indices)} distinct periodic-orbit states for "
            f"{len(steps)} Newton iterations; expected an initial state plus every update"
        )
    selected_indices = np.rint(
        np.linspace(1, len(distinct_indices) - 1, len(steps))
    ).astype(int)
    selected_field_indices = [distinct_indices[index] for index in selected_indices]
    iteration_states = [values[index] for index in selected_field_indices]
    if len(periods) == len(values):
        iteration_periods = [periods[index] for index in selected_field_indices]
    else:
        period_selection = np.rint(
            np.linspace(1, len(periods) - 1, len(steps))
        ).astype(int)
        iteration_periods = [periods[index] for index in period_selection]
    norm = Normalize(
        vmin=min(float(np.nanmin(item)) for item in values),
        vmax=max(float(np.nanmax(item)) for item in values),
    )
    positive_residuals = [step.residual for step in steps]
    positive_residuals.extend(value for step in steps for value in step.gmres_residuals)
    residual_min = min(value for value in positive_residuals if value > 0.0)
    residual_max = max(positive_residuals)
    colors = plt.get_cmap("turbo")(np.linspace(0.05, 0.95, len(steps)))
    max_gmres_iterations = max(len(step.gmres_residuals) for step in steps)

    convergence_dir = prepare_frames_directory(
        output_dir / "periodic_orbit_convergence_frames"
    )
    for last_index, (current_step, current_state, current_period) in enumerate(
        zip(steps, iteration_states, iteration_periods)
    ):
        visible_steps = steps[: last_index + 1]
        fig = plt.figure(figsize=(15, 9), constrained_layout=True)
        grid = fig.add_gridspec(2, 2, height_ratios=(1.12, 1.0))
        field_ax = fig.add_subplot(grid[0, :])
        newton_ax = fig.add_subplot(grid[1, 0])
        gmres_ax = fig.add_subplot(grid[1, 1])

        image = draw_slice(
            field_ax,
            geometry,
            current_state,
            cmap="viridis",
            norm=norm,
            xlim=xlim,
            ylim=ylim,
        )
        field_ax.set_title(
            f"{case} periodic orbit at phase zero after Newton iteration "
            f"{current_step.index} "
            f"({component}, T={current_period:.9g}, z={geometry.z:.3g})"
        )
        fig.colorbar(image, ax=field_ax, pad=0.015, label=component)

        newton_ax.semilogy(
            [step.index for step in visible_steps],
            [step.residual for step in visible_steps],
            "o-",
            color="black",
            lw=1.8,
        )
        newton_ax.set_xlim(0.5, len(steps) + 0.5)
        newton_ax.set_xlabel("Newton iteration")
        newton_ax.set_ylabel(r"$\|F(q,T)\|$")
        newton_ax.set_ylim(residual_min * 0.5, residual_max * 2.0)
        newton_ax.set_title("Periodic-orbit Newton convergence")
        newton_ax.grid(True, which="both", alpha=0.3)

        for color, step in zip(colors[: last_index + 1], visible_steps):
            gmres_ax.semilogy(
                np.arange(len(step.gmres_residuals)),
                step.gmres_residuals,
                color=color,
                lw=1.3,
                label=f"Newton {step.index}",
            )
            gmres_ax.axhline(
                step.tolerance, color=color, ls="--", lw=0.7, alpha=0.45
            )
        gmres_ax.set_xlim(0, max_gmres_iterations - 1)
        gmres_ax.set_xlabel("Inner Krylov iteration")
        gmres_ax.set_ylabel(r"$\|r\|$")
        gmres_ax.set_ylim(residual_min * 0.5, residual_max * 2.0)
        gmres_ax.set_title("GMRES convergence by Newton iteration")
        gmres_ax.grid(True, which="both", alpha=0.3)
        gmres_ax.legend(ncol=2, fontsize=8, loc="best")

        fig.suptitle(f"Cylinder gauntlet: {case}", fontsize=15)
        output = convergence_dir / f"periodic_orbit_convergence_{last_index + 1:04d}.png"
        fig.savefig(output, dpi=200)
        plt.close(fig)
    print(
        f"[{case}] wrote {len(steps)} periodic-orbit convergence frames "
        f"to {convergence_dir}"
    )

    final_period = periods[-1]
    errors = np.abs(periods - final_period)
    positive_errors = errors[errors > 0.0]
    error_min = float(np.min(positive_errors)) if positive_errors.size else 1.0e-16
    error_max = float(np.max(positive_errors)) if positive_errors.size else 1.0
    period_min = float(np.min(periods))
    period_max = float(np.max(periods))
    period_padding = max(0.05 * (period_max - period_min), 1.0e-10)
    period_dir = prepare_frames_directory(output_dir / "period_convergence_frames")
    for last_index in range(len(periods)):
        visible = slice(0, last_index + 1)
        visible_errors = errors[visible]
        visible_positive = visible_errors > 0.0
        fig, ax = plt.subplots(figsize=(8, 5.5), constrained_layout=True)
        ax.plot(
            period_iterations[visible],
            periods[visible],
            "o-",
            color="tab:blue",
            label=r"Period $T_n$",
        )
        ax.axhline(
            final_period,
            color="black",
            ls="--",
            lw=1.0,
            label=rf"$T={final_period:.9g}$",
        )
        ax.set_xlim(float(period_iterations[0]) - 0.5, float(period_iterations[-1]) + 0.5)
        ax.set_ylim(period_min - period_padding, period_max + period_padding)
        ax.set_xlabel("Saved outer iteration")
        ax.set_ylabel(r"Period estimate $T_n$", color="tab:blue")
        ax.tick_params(axis="y", labelcolor="tab:blue")
        ax.grid(True, alpha=0.3)
        error_ax = ax.twinx()
        error_ax.semilogy(
            period_iterations[visible][visible_positive],
            visible_errors[visible_positive],
            "s-",
            color="tab:red",
            label=r"$|T_n-T_{final}|$",
        )
        error_ax.set_ylim(error_min * 0.5, error_max * 2.0)
        error_ax.set_ylabel(r"Absolute period error $|T_n-T_{final}|$", color="tab:red")
        error_ax.tick_params(axis="y", labelcolor="tab:red")
        handles, labels = ax.get_legend_handles_labels()
        error_handles, error_labels = error_ax.get_legend_handles_labels()
        ax.legend(handles + error_handles, labels + error_labels, loc="best")
        ax.set_title(f"{case} periodic-orbit period convergence")
        output = period_dir / f"period_convergence_{last_index + 1:04d}.png"
        fig.savefig(output, dpi=220)
        plt.close(fig)
    print(f"[{case}] wrote {len(periods)} period frames to {period_dir}")


def read_eigenvalues(path: Path) -> list[EigenvalueRow]:
    rows: list[EigenvalueRow] = []
    for line in require_file(path).read_text(encoding="utf-8", errors="replace").splitlines():
        tokens = line.split()
        if len(tokens) != 6:
            continue
        try:
            row = EigenvalueRow(
                iteration=int(tokens[0]),
                value=complex(float(tokens[1]), float(tokens[2])),
                residual=float(tokens[4]),
                converged=tokens[5].upper() == "T",
            )
        except ValueError:
            continue
        rows.append(row)
    if not rows:
        raise ValueError(f"No eigenvalues found in {path}")
    final_iteration = max(row.iteration for row in rows)
    return [row for row in rows if row.iteration == final_iteration]


def read_end_time(case_dir: Path) -> float:
    case_files = sorted(case_dir.glob("*.case"))
    if len(case_files) != 1:
        raise ValueError(f"Expected one .case file in {case_dir}, found {len(case_files)}")
    data = json.loads(case_files[0].read_text(encoding="utf-8"))
    return float(data["case"]["time"]["end_time"])


def continuous_spectrum(rows: list[EigenvalueRow], interval: float) -> list[tuple[EigenvalueRow, complex]]:
    return [(row, cmath.log(row.value) / interval) for row in rows if abs(row.value) > 0.0]


def make_spectrum(case: str, case_results: Path, output_dir: Path) -> None:
    direct_dir = case_results / "LNS"
    adjoint_dir = case_results / "adjoint"
    direct = continuous_spectrum(
        read_eigenvalues(direct_dir / "eigs_output.txt"), read_end_time(direct_dir)
    )
    adjoint = continuous_spectrum(
        read_eigenvalues(adjoint_dir / "eigs_output.txt"), read_end_time(adjoint_dir)
    )

    fig, ax = plt.subplots(figsize=(8, 6), constrained_layout=True)
    for label, spectrum, color, marker in (
        ("Forward LNS", direct, "tab:blue", "o"),
        ("Adjoint LNS", adjoint, "tab:red", "x"),
    ):
        for index, (row, eigenvalue) in enumerate(spectrum):
            kwargs = {
                "marker": marker,
                "s": 65 if row.converged else 42,
                "color": color,
                "alpha": 1.0 if row.converged else 0.45,
                "label": label if index == 0 else None,
            }
            if marker == "o" and not row.converged:
                kwargs.update(facecolors="none", edgecolors=color)
                kwargs.pop("color")
            ax.scatter(
                eigenvalue.imag / (2.0 * math.pi),
                eigenvalue.real,
                **kwargs,
            )

    ax.axhline(0.0, color="0.25", lw=1.0, ls="--")
    ax.axvline(0.0, color="0.75", lw=0.8)
    ax.set_xlabel(r"Frequency $\mathrm{Im}(\lambda)/(2\pi)$")
    ax.set_ylabel(r"Growth rate $\mathrm{Re}(\lambda)$")
    ax.set_title(f"{case} forward and adjoint spectra")
    ax.grid(True, alpha=0.25)
    ax.legend()
    output = output_dir / "forward_adjoint_spectrum.png"
    fig.savefig(output, dpi=220)
    plt.close(fig)
    print(f"[{case}] wrote {output}")


def read_floquet_spectrum(path: Path) -> list[FloquetRow]:
    rows: list[FloquetRow] = []
    with require_file(path).open(newline="", encoding="utf-8") as stream:
        for row in csv.DictReader(stream):
            try:
                rows.append(
                    FloquetRow(
                        index=int(row["index"]),
                        multiplier=complex(
                            float(row["multiplier_real"]), float(row["multiplier_imag"])
                        ),
                        exponent=complex(
                            float(row["exponent_real"]), float(row["exponent_imag"])
                        ),
                        residual=float(row["residual"]),
                    )
                )
            except (KeyError, TypeError, ValueError) as error:
                raise ValueError(f"Malformed Floquet-spectrum row in {path}: {row}") from error
    if not rows:
        raise ValueError(f"No Floquet spectrum found in {path}")
    return rows


def make_floquet_spectrum(
    case: str,
    case_results: Path,
    output_dir: Path,
    dimensionality: str,
) -> None:
    if dimensionality == "2D":
        direct_dir = case_results / "LNS_floquete"
        adjoint_dir = case_results / "adjoint_floquete"
        adjoint_file = "adjoint_floquet_spectrum.csv"
        output_name = "direct_adjoint_floquet_spectrum.png"
    elif dimensionality == "3D":
        direct_dir = case_results / "LNS_floquete_3D"
        adjoint_dir = case_results / "adjoint_floquete_3D"
        adjoint_file = "adjoint_floquet_3D_spectrum.csv"
        output_name = "direct_adjoint_floquet_3D_spectrum.png"
    else:
        raise ValueError(f"Unknown Floquet dimensionality: {dimensionality}")

    spectra = (
        (
            f"Direct {dimensionality} LNS",
            read_floquet_spectrum(direct_dir / "floquet_spectrum.csv"),
            "tab:blue",
            "o",
        ),
        (
            f"Adjoint {dimensionality} LNS",
            read_floquet_spectrum(adjoint_dir / adjoint_file),
            "tab:red",
            "x",
        ),
    )
    fig, (multiplier_ax, exponent_ax) = plt.subplots(
        1, 2, figsize=(13, 5.8), constrained_layout=True
    )

    theta = np.linspace(0.0, 2.0 * math.pi, 600)
    multiplier_ax.plot(np.cos(theta), np.sin(theta), "--", color="0.45", lw=1.0)
    for label, rows, color, marker in spectra:
        multipliers = np.asarray([row.multiplier for row in rows])
        exponents = np.asarray([row.exponent for row in rows])
        multiplier_ax.scatter(
            multipliers.real,
            multipliers.imag,
            s=68,
            marker=marker,
            color=color,
            label=label,
            zorder=3,
        )
        exponent_ax.scatter(
            exponents.imag / (2.0 * math.pi),
            exponents.real,
            s=68,
            marker=marker,
            color=color,
            label=label,
            zorder=3,
        )
        for row in rows:
            multiplier_ax.annotate(
                str(row.index),
                (row.multiplier.real, row.multiplier.imag),
                xytext=(5, 5),
                textcoords="offset points",
                fontsize=8,
                color=color,
            )

    multiplier_ax.axhline(0.0, color="0.75", lw=0.8)
    multiplier_ax.axvline(0.0, color="0.75", lw=0.8)
    multiplier_ax.set_aspect("equal", adjustable="box")
    multiplier_ax.set_xlabel(r"$\mathrm{Re}(\mu)$")
    multiplier_ax.set_ylabel(r"$\mathrm{Im}(\mu)$")
    multiplier_ax.set_title("Floquet multipliers (unit circle dashed)")
    multiplier_ax.grid(True, alpha=0.25)
    multiplier_ax.legend()

    exponent_ax.axhline(0.0, color="0.25", lw=1.0, ls="--")
    exponent_ax.axvline(0.0, color="0.75", lw=0.8)
    exponent_ax.set_xlabel(r"Frequency $\mathrm{Im}(\lambda)/(2\pi)$")
    exponent_ax.set_ylabel(r"Growth rate $\mathrm{Re}(\lambda)$")
    exponent_ax.set_title("Floquet exponents")
    exponent_ax.grid(True, alpha=0.25)
    exponent_ax.legend()
    fig.suptitle(
        f"{case} direct and adjoint {dimensionality} Floquet spectra",
        fontsize=15,
    )
    output = output_dir / output_name
    fig.savefig(output, dpi=220)
    plt.close(fig)
    print(f"[{case}] wrote {output}")


def make_available_floquet_spectra(
    case: str, case_results: Path, output_dir: Path
) -> None:
    configurations = (
        (
            "2D",
            case_results / "LNS_floquete" / "floquet_spectrum.csv",
            case_results / "adjoint_floquete" / "adjoint_floquet_spectrum.csv",
        ),
        (
            "3D",
            case_results / "LNS_floquete_3D" / "floquet_spectrum.csv",
            case_results
            / "adjoint_floquete_3D"
            / "adjoint_floquet_3D_spectrum.csv",
        ),
    )
    found = False
    for dimensionality, direct_file, adjoint_file in configurations:
        if not direct_file.exists() and not adjoint_file.exists():
            continue
        require_file(direct_file)
        require_file(adjoint_file)
        make_floquet_spectrum(case, case_results, output_dir, dimensionality)
        found = True
    if not found:
        raise FileNotFoundError(
            f"No direct/adjoint Floquet spectrum results found for {case}"
        )


def read_vorticity_sequence(
    files: list[Path],
) -> tuple[SliceMesh, list[np.ndarray], list[float]]:
    """Read velocity snapshots and compute omega_z with SEM derivatives."""
    if not files:
        raise ValueError("No field files were supplied for vorticity")

    comm = MPI.COMM_SELF
    mesh = Mesh(comm)
    field = FieldRegistry(comm)
    geometry: SliceMesh | None = None
    coefficient: Coef | None = None
    connectivity: MeshConnectivity | None = None
    plane = 0
    vorticities: list[np.ndarray] = []
    times: list[float] = []

    for file_index, path in enumerate(files):
        pynekread(
            str(path),
            comm,
            data_dtype=np.float64,
            msh=mesh if file_index == 0 else None,
            fld=field,
            overwrite_fld=True,
        )
        if file_index == 0:
            geometry, plane = build_slice_mesh(mesh)
            coefficient = Coef(mesh, comm, get_area=False)
            # Average only duplicate interface nodes after differentiation;
            # element-interior spectral values are left unchanged.
            connectivity = MeshConnectivity(comm, msh=mesh)
        if "u" not in field.registry or "v" not in field.registry:
            raise KeyError(f"Velocity components u and v are required in {path}")
        assert coefficient is not None and connectivity is not None
        dvdx = coefficient.dudxyz(
            field.registry["v"],
            coefficient.drdx,
            coefficient.dsdx,
            coefficient.dtdx,
        )
        dudy = coefficient.dudxyz(
            field.registry["u"],
            coefficient.drdy,
            coefficient.dsdy,
            coefficient.dtdy,
        )
        vorticity = connectivity.dssum(
            field=dvdx - dudy, msh=mesh, average="multiplicity"
        )
        vorticities.append(
            np.asarray(vorticity[:, plane, :, :], dtype=np.float64)
            .reshape(-1)
            .copy()
        )
        times.append(float(field.t))

    assert geometry is not None
    return geometry, vorticities, times


def make_shedding_frames(
    case: str,
    case_results: Path,
    output_dir: Path,
    xlim: tuple[float, float],
    ylim: tuple[float, float],
    vorticity_limit: float | None,
) -> None:
    movie_dir = (
        case_results
        / "periodic_orbit"
        / "viz"
        / f"periodic_orbit_{case}_movie"
    )
    files = sorted(movie_dir.glob("*.f[0-9]*"), key=field_sort_key)
    if not files:
        raise FileNotFoundError(f"No periodic-orbit movie fields found in {movie_dir}")
    geometry, vorticities, times = read_vorticity_sequence(files)

    if vorticity_limit is None:
        per_frame_limits = [float(np.nanpercentile(np.abs(item), 99.5)) for item in vorticities]
        vmax = max(per_frame_limits)
    else:
        vmax = float(vorticity_limit)
    if not np.isfinite(vmax) or vmax <= 0.0:
        raise ValueError(f"{case} vorticity has zero or invalid amplitude")
    norm = TwoSlopeNorm(vmin=-vmax, vcenter=0.0, vmax=vmax)

    frames_dir = prepare_frames_directory(output_dir / "shedding_vorticity_frames")
    for frame_index, (vorticity, time) in enumerate(zip(vorticities, times)):
        fig, ax = plt.subplots(figsize=(12, 5.6), constrained_layout=True)
        image = draw_slice(
            ax,
            geometry,
            vorticity,
            cmap="RdBu_r",
            norm=norm,
            xlim=xlim,
            ylim=ylim,
        )
        ax.set_title(
            rf"{case} cylinder shedding: $\omega_z=\partial_x v-\partial_y u$, "
            rf"frame {frame_index:03d}, $t={time:.6g}$"
        )
        fig.colorbar(image, ax=ax, pad=0.015, label=r"Spanwise vorticity $\omega_z$")
        frame_path = frames_dir / f"vorticity_{frame_index:04d}.png"
        fig.savefig(frame_path, dpi=150)
        plt.close(fig)
    print(f"[{case}] wrote {len(vorticities)} vorticity frames to {frames_dir}")

def make_mode_frames(
    case: str,
    case_results: Path,
    output_dir: Path,
    kind: str,
    component: str,
    xlim: tuple[float, float],
    ylim: tuple[float, float],
    frames: int,
) -> None:
    if kind == "direct":
        mode_dir = case_results / "LNS" / "viz" / "mode_direct"
        stem = "mode_direct0"
        label = "Forward LNS"
    elif kind == "adjoint":
        mode_dir = case_results / "adjoint" / "viz" / "mode_adjoint"
        stem = "mode_adjoint0"
        label = "Adjoint LNS"
    else:
        raise ValueError(f"Unknown mode kind: {kind}")

    mode_files = [
        require_file(mode_dir / f"{stem}.f00000"),
        require_file(mode_dir / f"{stem}.f00001"),
    ]
    geometry, values, _ = read_slice_sequence(mode_files, component)
    real_part, imaginary_part = values
    amplitude = np.hypot(real_part, imaginary_part)
    vmax = float(np.nanpercentile(amplitude, 99.5))
    if not np.isfinite(vmax) or vmax <= 0.0:
        raise ValueError(f"{case} {kind} mode has zero or invalid amplitude")
    norm = TwoSlopeNorm(vmin=-vmax, vcenter=0.0, vmax=vmax)
    phases = np.linspace(0.0, 2.0 * math.pi, frames, endpoint=False)

    frames_dir = prepare_frames_directory(output_dir / f"leading_{kind}_mode_frames")
    for frame_index, phase in enumerate(phases):
        fig, ax = plt.subplots(figsize=(12, 5.6), constrained_layout=True)
        mode = real_part * math.cos(phase) - imaginary_part * math.sin(phase)
        image = draw_slice(
            ax,
            geometry,
            mode,
            cmap="RdBu_r",
            norm=norm,
            xlim=xlim,
            ylim=ylim,
        )
        ax.set_title(
            f"{case} {label} leading mode: {component}, "
            rf"$\phi={phase / math.pi:.2f}\pi$, z={geometry.z:.3g}"
        )
        fig.colorbar(image, ax=ax, pad=0.015, label=f"{component} perturbation")
        output = frames_dir / f"leading_{kind}_mode_{frame_index:04d}.png"
        fig.savefig(output, dpi=150)
        plt.close(fig)
    print(f"[{case}] wrote {frames} {kind}-mode frames to {frames_dir}")


def read_singular_spectrum(path: Path) -> list[tuple[int, float, float]]:
    rows: list[tuple[int, float, float]] = []
    for line in require_file(path).read_text(encoding="utf-8", errors="replace").splitlines():
        tokens = line.split()
        if len(tokens) != 3:
            continue
        try:
            rows.append((int(tokens[0]), float(tokens[1]), float(tokens[2])))
        except ValueError:
            continue
    if not rows:
        raise ValueError(f"No singular values found in {path}")
    return rows


def make_transient_growth_plots(
    case: str,
    case_results: Path,
    output_dir: Path,
    component: str,
    xlim: tuple[float, float],
    ylim: tuple[float, float],
) -> None:
    transient_dir = case_results / "transient_growth"
    if not transient_dir.is_dir():
        print(f"[{case}] no transient-growth results; skipping singular vectors")
        return

    spectrum = read_singular_spectrum(transient_dir / "singular_spectrum.dat")
    input_files = sorted((transient_dir / "viz" / "prt").glob("*.f[0-9]*"), key=field_sort_key)
    output_files = sorted((transient_dir / "viz" / "rsp").glob("*.f[0-9]*"), key=field_sort_key)
    if not input_files or not output_files:
        raise FileNotFoundError(f"Transient-growth prt/rsp fields are missing in {transient_dir}")

    input_geometry, input_vectors, _ = read_slice_sequence(input_files, component)
    output_geometry, output_vectors, _ = read_slice_sequence(output_files, component)
    if not (
        np.allclose(input_geometry.x, output_geometry.x)
        and np.allclose(input_geometry.y, output_geometry.y)
    ):
        raise ValueError(f"{case}: transient-growth input and output meshes differ")

    count = min(len(spectrum), len(input_vectors), len(output_vectors))
    horizon = read_end_time(transient_dir)

    indices = [row[0] for row in spectrum[:count]]
    singular_values = [row[1] for row in spectrum[:count]]
    fig, ax = plt.subplots(figsize=(7.5, 5.2), constrained_layout=True)
    bars = ax.bar(indices, singular_values, color="tab:blue", alpha=0.82)
    ax.set_xticks(indices)
    ax.set_xlabel("Singular triplet index")
    ax.set_ylabel(r"Singular value $\sigma_i$")
    ax.set_title(f"{case} transient growth at T={horizon:g}")
    ax.grid(True, axis="y", alpha=0.25)
    for bar, sigma in zip(bars, singular_values):
        ax.text(
            bar.get_x() + 0.5 * bar.get_width(),
            bar.get_height(),
            rf"$G={sigma**2:.3f}$",
            ha="center",
            va="bottom",
            fontsize=9,
        )
    spectrum_output = output_dir / "transient_growth_singular_values.png"
    fig.savefig(spectrum_output, dpi=220)
    plt.close(fig)
    print(f"[{case}] wrote {spectrum_output}")

    for vector_index in range(count):
        rank, sigma, residual = spectrum[vector_index]
        optimal_input = input_vectors[vector_index]
        # The SVD convention is A V_i = sigma_i U_i.  rsp stores normalized
        # U_i, so multiplying by sigma gives the actual output generated by a
        # unit-energy input V_i.
        amplified_output = sigma * output_vectors[vector_index]
        vmax = float(
            np.nanpercentile(
                np.abs(np.concatenate((optimal_input, amplified_output))),
                99.5,
            )
        )
        if not np.isfinite(vmax) or vmax <= 0.0:
            raise ValueError(f"{case} singular vector {rank} has zero or invalid amplitude")
        norm = TwoSlopeNorm(vmin=-vmax, vcenter=0.0, vmax=vmax)

        fig, axes = plt.subplots(1, 2, figsize=(14, 5.3), constrained_layout=True)
        input_image = draw_slice(
            axes[0],
            input_geometry,
            optimal_input,
            cmap="RdBu_r",
            norm=norm,
            xlim=xlim,
            ylim=ylim,
        )
        draw_slice(
            axes[1],
            output_geometry,
            amplified_output,
            cmap="RdBu_r",
            norm=norm,
            xlim=xlim,
            ylim=ylim,
        )
        qualifier = "optimal" if rank == 1 else f"rank-{rank}"
        axes[0].set_title(rf"{qualifier.capitalize()} input $V_{{{rank}}}$ at $t=0$")
        axes[1].set_title(
            rf"Amplified output $A V_{{{rank}}}=\sigma_{{{rank}}}U_{{{rank}}}$ at $t={horizon:g}$"
        )
        fig.colorbar(
            input_image,
            ax=axes,
            pad=0.015,
            label=f"{component} perturbation (common scale)",
        )
        fig.suptitle(
            rf"{case} transient-growth triplet {rank}: "
            rf"$\sigma={sigma:.6g}$, $G=\sigma^2={sigma**2:.6g}$, residual={residual:.2e}"
        )
        output = output_dir / f"transient_growth_singular_vector_{rank:02d}.png"
        fig.savefig(output, dpi=200)
        plt.close(fig)
        print(f"[{case}] wrote {output}")


def main() -> None:
    args = parse_args()
    logging.getLogger().setLevel(logging.WARNING)
    results = args.results.resolve()
    output_root = args.output.resolve()
    xlim = tuple(args.xlim)
    ylim = tuple(args.ylim)

    for case in args.reynolds:
        case_results = results / case
        if not case_results.is_dir():
            raise FileNotFoundError(f"Gauntlet results are missing for {case}: {case_results}")
        output_dir = output_root / case
        output_dir.mkdir(parents=True, exist_ok=True)

        if case in ("Re180", "Re200"):
            if not args.skip_baseflow:
                make_periodic_orbit_convergence(
                    case,
                    case_results,
                    output_dir,
                    args.component,
                    xlim,
                    ylim,
                )
            if not args.skip_shedding:
                make_shedding_frames(
                    case,
                    case_results,
                    output_dir,
                    xlim,
                    ylim,
                    args.vorticity_limit,
                )
            if not args.skip_spectra:
                make_available_floquet_spectra(case, case_results, output_dir)
            continue

        if not args.skip_baseflow:
            make_baseflow_convergence(
                case,
                case_results,
                output_dir,
                args.component,
                xlim,
                ylim,
            )
        if not args.skip_spectra:
            make_spectrum(case, case_results, output_dir)
        if not args.skip_modes:
            make_mode_frames(
                case,
                case_results,
                output_dir,
                "direct",
                args.component,
                xlim,
                ylim,
                args.mode_frames,
            )
            make_mode_frames(
                case,
                case_results,
                output_dir,
                "adjoint",
                args.component,
                xlim,
                ylim,
                args.mode_frames,
            )
        if not args.skip_transient_growth:
            make_transient_growth_plots(
                case,
                case_results,
                output_dir,
                args.component,
                xlim,
                ylim,
            )


if __name__ == "__main__":
    main()
