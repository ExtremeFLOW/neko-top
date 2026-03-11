#!/usr/bin/env python3
"""Plot SEM option regression optimization histories."""

from __future__ import annotations

import argparse
import csv
import re
from pathlib import Path
from typing import Dict, List, Optional, Tuple

import matplotlib.pyplot as plt


def normalize_header(header: str) -> str:
    """Normalize a header for robust matching."""
    cleaned = re.sub(r"[^a-z0-9]+", " ", header.lower())
    return re.sub(r"\s+", " ", cleaned).strip()


def find_header(
    headers: List[str], required_tokens: List[str]
) -> Optional[str]:
    """Find the first header that contains all required tokens."""
    for header in headers:
        norm = normalize_header(header)
        if all(token in norm for token in required_tokens):
            return header
    return None


def find_header_by_normalized_name(
    headers: List[str], normalized_name: str
) -> Optional[str]:
    """Find a header by exact normalized name."""
    for header in headers:
        if normalize_header(header) == normalized_name:
            return header
    return None


def resolve_data_file(case_dir: Path) -> Optional[Path]:
    """Resolve optimization data file, preferring .txt over .csv."""
    txt_path = case_dir / "optimization_data.txt"
    csv_path = case_dir / "optimization_data.csv"

    if txt_path.exists():
        return txt_path
    if csv_path.exists():
        return csv_path
    return None


def read_case_history(
    data_path: Path,
) -> Tuple[List[float], List[float], List[float]]:
    """Read iteration, total objective and volume columns."""
    with data_path.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.reader(handle)
        try:
            headers = next(reader)
        except StopIteration as exc:
            raise ValueError(f"Empty file: {data_path}") from exc

        iter_col = find_header(headers, ["iter"])
        obj_col = find_header(headers, ["total", "objective"])
        vol_col = find_header_by_normalized_name(
            headers, "volume constraint volume"
        )
        if vol_col is None:
            vol_col = find_header(headers, ["volume", "constraint"])

        if iter_col is None:
            raise ValueError(f"Could not find iteration column in {data_path}")
        if obj_col is None:
            raise ValueError(
                f"Could not find total objective column in {data_path}"
            )
        if vol_col is None:
            raise ValueError(f"Could not find volume column in {data_path}")

        i_iter = headers.index(iter_col)
        i_obj = headers.index(obj_col)
        i_vol = headers.index(vol_col)

        iters: List[float] = []
        objective: List[float] = []
        volume: List[float] = []
        for row in reader:
            if not row or all(cell.strip() == "" for cell in row):
                continue
            iters.append(float(row[i_iter]))
            objective.append(float(row[i_obj]))
            volume.append(float(row[i_vol]))

    if not iters:
        raise ValueError(f"No rows found in {data_path}")

    return iters, objective, volume


def collect_histories(
    results_dir: Path,
) -> Dict[str, Tuple[List[float], List[float], List[float]]]:
    """Collect histories from all case directories."""
    histories: Dict[str, Tuple[List[float], List[float], List[float]]] = {}
    for case_dir in sorted(results_dir.glob("sem_map_option_*")):
        if not case_dir.is_dir():
            continue
        data_file = resolve_data_file(case_dir)
        if data_file is None:
            print(f"Skipping {case_dir.name}: no optimization_data file")
            continue
        try:
            histories[case_dir.name] = read_case_history(data_file)
        except ValueError as exc:
            print(f"Skipping {case_dir.name}: {exc}")
    return histories


def make_plot(
    histories: Dict[str, Tuple[List[float], List[float], List[float]]],
    output_path: Path,
    volume_limit: float,
    show: bool,
) -> None:
    """Create and save the two-panel optimization plot."""
    fig, (ax_left, ax_right) = plt.subplots(1, 2, figsize=(14, 5))

    for label, (iters, objective, volume) in histories.items():
        ax_left.plot(iters, objective, label=label)
        ax_right.plot(iters, volume, label=label)

    ax_left.set_title("Total objective function")
    ax_left.set_xlabel("Iteration")
    ax_left.set_ylabel("Objective")
    ax_left.grid(True, alpha=0.3)
    ax_left.legend(fontsize=8)

    ax_right.set_title("Volume constraint.volume")
    ax_right.axhline(
        volume_limit,
        color="black",
        linestyle="--",
        linewidth=1.0,
        label=f"limit = {volume_limit}",
    )
    ax_right.set_xlabel("Iteration")
    ax_right.set_ylabel("Volume")
    ax_right.grid(True, alpha=0.3)
    ax_right.legend(fontsize=8)

    fig.tight_layout()
    output_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(output_path, dpi=150)
    print(f"Wrote {output_path}")

    if show:
        plt.show()


def parse_args() -> argparse.Namespace:
    """Parse command-line arguments."""
    default_dir = Path(__file__).resolve().parent / "results"
    default_output = default_dir / "optimization_history_sem_options.png"

    parser = argparse.ArgumentParser(
        description="Plot objective and volume histories for SEM options."
    )
    parser.add_argument(
        "--results-dir",
        type=Path,
        default=default_dir,
        help="Directory containing sem_map_option_* case folders.",
    )
    parser.add_argument(
        "--out",
        type=Path,
        default=default_output,
        help="Output image path.",
    )
    parser.add_argument(
        "--volume-limit",
        type=float,
        default=0.2,
        help="Horizontal limit line for volume plot.",
    )
    parser.add_argument(
        "--show",
        action="store_true",
        help="Show the plot interactively.",
    )
    return parser.parse_args()


def main() -> None:
    """Entrypoint."""
    args = parse_args()
    histories = collect_histories(args.results_dir)
    if not histories:
        raise SystemExit(
            f"No valid case histories found under {args.results_dir}"
        )

    make_plot(histories, args.out, args.volume_limit, args.show)


if __name__ == "__main__":
    main()
