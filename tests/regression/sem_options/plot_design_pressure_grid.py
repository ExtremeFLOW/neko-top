#!/usr/bin/env python3
"""Plot pressure slices from SEM-option design files in a grid."""

from __future__ import annotations

import argparse
import re
from pathlib import Path
from typing import List, Optional, Tuple

import matplotlib.pyplot as plt
import matplotlib.tri as tri
import numpy as np


def parse_iterations(text: str) -> List[int]:
    """Parse a comma-separated iteration list."""
    values: List[int] = []
    for chunk in text.split(","):
        item = chunk.strip()
        if not item:
            continue
        values.append(int(item))
    if not values:
        raise ValueError("No iterations provided")
    return values


def detect_index_width(case_dir: Path, prefix: str) -> int:
    """Detect the field index width from existing files."""
    for path in sorted(case_dir.glob(f"{prefix}.f*")):
        match = re.search(r"\.f(\d+)$", path.name)
        if match:
            return len(match.group(1))
    return 5


def gather_flat(comm, values: np.ndarray) -> Optional[np.ndarray]:
    """Gather flattened data from all ranks onto rank 0."""
    gathered = comm.gather(values, root=0)
    if comm.Get_rank() == 0:
        return np.concatenate(gathered)
    return None


def build_geometry(
    x_all: np.ndarray,
    y_all: np.ndarray,
) -> Tuple[
    Optional[tri.Triangulation],
    Optional[np.ndarray],
    Optional[np.ndarray],
    Optional[str],
]:
    """Build robust xy geometry for repeated scalar projection."""
    coords = np.column_stack((x_all, y_all))
    unique_coords, inverse = np.unique(coords, axis=0, return_inverse=True)

    unique_count = unique_coords.shape[0]
    if unique_count < 3:
        return None, None, None, f"only {unique_count} unique points"

    # Guard against collinear point clouds.
    centered = unique_coords - unique_coords.mean(axis=0, keepdims=True)
    if np.linalg.matrix_rank(centered) < 2:
        return None, None, None, "points are collinear in xy"

    try:
        triang = tri.Triangulation(unique_coords[:, 0], unique_coords[:, 1])
    except Exception as exc:  # pragma: no cover
        return None, None, None, str(exc)

    counts = np.bincount(inverse)
    return triang, inverse, counts, None


def project_to_geometry(
    values: np.ndarray,
    inverse: np.ndarray,
    counts: np.ndarray,
) -> np.ndarray:
    """Average scalar values at duplicate xy coordinates."""
    sums = np.bincount(inverse, weights=values, minlength=counts.size)
    return sums / np.maximum(counts, 1)


def get_k_index(mesh, _k_index: Optional[int]) -> int:
    """Always use the center z-slice for xy plotting."""
    return int(mesh.lz // 2)


def load_case_geometry(
    comm,
    mesh_path: Path,
    k_index: Optional[int],
) -> Tuple[
    Optional[tri.Triangulation],
    Optional[np.ndarray],
    Optional[np.ndarray],
    Optional[int],
    Optional[str],
]:
    """Load mesh geometry once for all iterations in one case row."""
    from pysemtools.io.ppymech.neksuite import preadnek
    from pysemtools.datatypes.msh import Mesh

    mesh_data = preadnek(str(mesh_path), comm)
    mesh = Mesh(comm, data=mesh_data)
    k = get_k_index(mesh, k_index)

    x_local = mesh.x[:, k, :, :].reshape(-1)
    y_local = mesh.y[:, k, :, :].reshape(-1)
    x_all = gather_flat(comm, x_local)
    y_all = gather_flat(comm, y_local)

    if comm.Get_rank() == 0:
        triang, inverse, counts, err = build_geometry(x_all, y_all)
        return triang, inverse, counts, k, err

    return None, None, None, k, None


def read_pressure_xy_slice(
    comm,
    file_path: Path,
    k_index: int,
    inverse: Optional[np.ndarray],
    counts: Optional[np.ndarray],
) -> Optional[np.ndarray]:
    """Read pressure slice and map it onto precomputed xy geometry."""
    from pysemtools.io.ppymech.neksuite import preadnek
    from pysemtools.datatypes.field import Field

    data = preadnek(str(file_path), comm)
    fields = Field(comm, data=data)
    if len(fields.fields["pres"]) < 1:
        raise ValueError(f"No pressure field found in {file_path}")

    p_local = fields.fields["pres"][0][:, k_index, :, :].reshape(-1)
    p_all = gather_flat(comm, p_local)

    if comm.Get_rank() == 0:
        if inverse is None or counts is None:
            raise ValueError("Missing geometry mapping arrays on root")
        return project_to_geometry(p_all, inverse, counts)

    return None


def make_case_label(case_dir: Path) -> str:
    """Create a short row label from the folder name."""
    label = case_dir.name
    label = label.replace("sem_map_option_", "option ")
    label = label.replace("sem_option_", "option ")
    return label


def parse_args() -> argparse.Namespace:
    """Parse command-line arguments."""
    script_dir = Path(__file__).resolve().parent
    default_results = script_dir / "results"
    default_output = default_results / "design_pressure_xy_grid.png"

    parser = argparse.ArgumentParser(
        description="Plot pressure xy slices for multiple SEM-option cases."
    )
    parser.add_argument(
        "--results-dir",
        type=Path,
        default=default_results,
        help="Directory containing sem_map_option_* folders.",
    )
    parser.add_argument(
        "--case-pattern",
        type=str,
        default="sem_map_option_*",
        help="Glob pattern for case folders.",
    )
    parser.add_argument(
        "--iterations",
        type=str,
        default="1,2,3,5,10,20",
        help="Comma-separated iteration numbers.",
    )
    parser.add_argument(
        "--prefix",
        type=str,
        default="design0",
        help="Field-file prefix (default: design0).",
    )
    parser.add_argument(
        "--k-index",
        type=int,
        default=None,
        help="Deprecated: ignored. Script always uses center z-slice.",
    )
    parser.add_argument(
        "--vmin",
        type=float,
        default=0.0,
        help="Colorbar minimum (default: 0.0).",
    )
    parser.add_argument(
        "--vmax",
        type=float,
        default=1.0,
        help="Colorbar maximum (default: 1.0).",
    )
    parser.add_argument(
        "--levels",
        type=int,
        default=64,
        help="Number of contour levels.",
    )
    parser.add_argument(
        "--out",
        type=Path,
        default=default_output,
        help="Output image path.",
    )
    parser.add_argument(
        "--dpi",
        type=int,
        default=150,
        help="Output figure dpi.",
    )
    parser.add_argument(
        "--show",
        action="store_true",
        help="Show figure interactively on rank 0.",
    )
    return parser.parse_args()


def main() -> None:
    """Entrypoint."""
    args = parse_args()
    iterations = parse_iterations(args.iterations)

    try:
        from mpi4py import MPI
    except ModuleNotFoundError as exc:
        raise SystemExit("mpi4py is required to run this script") from exc

    comm = MPI.COMM_WORLD
    rank = comm.Get_rank()

    case_dirs = sorted(
        path
        for path in args.results_dir.glob(args.case_pattern)
        if path.is_dir()
    )
    if not case_dirs:
        raise SystemExit(
            f"No case folders found in {args.results_dir} "
            f"with pattern {args.case_pattern}"
        )

    n_rows = len(case_dirs)
    n_cols = len(iterations)
    fig = None
    axes = None

    if rank == 0:
        fig_w = max(3.0 * n_cols, 8.0)
        fig_h = max(2.2 * n_rows, 4.0)
        fig, axes = plt.subplots(
            n_rows,
            n_cols,
            figsize=(fig_w, fig_h),
            squeeze=False,
            constrained_layout=True,
        )

    levels = np.linspace(args.vmin, args.vmax, max(args.levels, 2))

    for row_idx, case_dir in enumerate(case_dirs):
        index_width = detect_index_width(case_dir, args.prefix)
        row_label = make_case_label(case_dir)
        mesh_file = case_dir / f"{args.prefix}.f{0:0{index_width}d}"

        row_triang = None
        row_inverse = None
        row_counts = None
        row_k = None
        row_geometry_error = None

        if not mesh_file.exists():
            row_geometry_error = "missing mesh file"
        else:
            try:
                row_geometry = load_case_geometry(comm, mesh_file, args.k_index)
                (
                    row_triang,
                    row_inverse,
                    row_counts,
                    row_k,
                    row_geometry_error,
                ) = row_geometry
            except Exception as exc:  # pragma: no cover
                row_geometry_error = str(exc)

        for col_idx, iteration in enumerate(iterations):
            file_name = f"{args.prefix}.f{iteration:0{index_width}d}"
            file_path = case_dir / file_name

            pressure = None
            err_text = None

            if row_geometry_error is not None:
                err_text = f"mesh error: {row_geometry_error}"

            if err_text is None and not file_path.exists():
                err_text = "missing"
            if err_text is None:
                try:
                    pressure = read_pressure_xy_slice(
                        comm,
                        file_path,
                        row_k,
                        row_inverse,
                        row_counts,
                    )
                except Exception as exc:  # pragma: no cover
                    err_text = str(exc)

            if rank != 0:
                continue

            ax = axes[row_idx, col_idx]
            ax.set_xticks([])
            ax.set_yticks([])
            ax.set_aspect("equal", adjustable="box")

            if row_idx == 0:
                ax.set_title(f"iter {iteration}")

            if col_idx == 0:
                ax.text(
                    -0.06,
                    0.5,
                    row_label,
                    transform=ax.transAxes,
                    ha="right",
                    va="center",
                    fontsize=9,
                )

            if row_triang is None or pressure is None:
                ax.set_facecolor("#f2f2f2")
                ax.text(
                    0.5,
                    0.5,
                    err_text or "unavailable",
                    transform=ax.transAxes,
                    ha="center",
                    va="center",
                    fontsize=8,
                )
                continue

            try:
                ax.tricontourf(
                    row_triang,
                    pressure,
                    levels=levels,
                    cmap="viridis",
                    vmin=args.vmin,
                    vmax=args.vmax,
                )
            except Exception as exc:  # pragma: no cover
                ax.set_facecolor("#f2f2f2")
                ax.text(
                    0.5,
                    0.5,
                    f"plot failed:\n{exc}",
                    transform=ax.transAxes,
                    ha="center",
                    va="center",
                    fontsize=8,
                )

    if rank != 0:
        return

    mappable = plt.cm.ScalarMappable(
        norm=plt.Normalize(vmin=args.vmin, vmax=args.vmax),
        cmap="viridis",
    )
    fig.colorbar(
        mappable,
        ax=list(axes.ravel()),
        fraction=0.02,
        pad=0.02,
    )

    args.out.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(args.out, dpi=args.dpi)
    print(f"Wrote {args.out}")

    if args.show:
        plt.show()
    plt.close(fig)


if __name__ == "__main__":
    main()
