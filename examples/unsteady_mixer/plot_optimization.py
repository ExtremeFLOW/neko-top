#!/usr/bin/env python3
import argparse
import csv
import json
import re
from pathlib import Path

import matplotlib.pyplot as plt


def _normalize_header(header):
    return re.sub(r"\s+", " ", header.strip().lower())


def _find_header(headers, required_tokens):
    for header in headers:
        norm = _normalize_header(header)
        if all(token in norm for token in required_tokens):
            return header
    return None


def _load_case_target(case_path):
    text = case_path.read_text()
    # Strip trailing commas before closing braces/brackets.
    cleaned = re.sub(r",\s*([}\]])", r"\1", text)
    data = json.loads(cleaned)
    objectives = data.get("optimization", {}).get("objectives", [])
    for obj in objectives:
        if obj.get("type") == "target_dissipation":
            return obj.get("target")
    raise ValueError("No target_dissipation objective found in case file.")


def _find_case_file(directory):
    case_files = list(directory.glob("*.case"))
    if len(case_files) != 1:
        raise FileNotFoundError(
            f"Expected exactly one .case file in {directory}, found {len(case_files)}."
        )
    return case_files[0]


def _resolve_data_path(directory, data_path):
    if data_path.is_absolute():
        return data_path

    resolved = directory / data_path
    if resolved.exists():
        return resolved

    if data_path.name == "optimization_data.csv":
        txt_path = directory / "optimization_data.txt"
        if txt_path.exists():
            return txt_path

    return resolved


def _load_csv_data(csv_path):
    with csv_path.open(newline="") as f:
        reader = csv.reader(f)
        try:
            headers = next(reader)
        except StopIteration:
            raise ValueError(f"{csv_path} is empty.")

        idx_iter = _find_header(headers, ["iter"])
        idx_total = _find_header(headers, ["total", "objective"])
        idx_ratio = _find_header(headers, ["target", "dissipation", "ratio"])
        scaled_headers = [
            header for header in headers if "scaled" in _normalize_header(header)
        ]

        if idx_iter is None:
            raise ValueError("Could not find iteration column in CSV.")
        if idx_total is None:
            raise ValueError("Could not find total objective column in CSV.")
        if idx_ratio is None:
            raise ValueError("Could not find target dissipation ratio column in CSV.")

        idx_iter = headers.index(idx_iter)
        idx_total = headers.index(idx_total)
        idx_ratio = headers.index(idx_ratio)
        scaled_indices = [headers.index(header) for header in scaled_headers]

        iters = []
        total = []
        ratio = []
        scaled_series = [[] for _ in scaled_indices]

        for row in reader:
            if not row or all(cell.strip() == "" for cell in row):
                continue
            iters.append(float(row[idx_iter]))
            total.append(float(row[idx_total]))
            ratio.append(float(row[idx_ratio]))
            for i, idx_scaled in enumerate(scaled_indices):
                scaled_series[i].append(float(row[idx_scaled]))

    if not total:
        raise ValueError("No data rows found in CSV.")
    return iters, total, ratio, scaled_headers, scaled_series


def _normalization_from_total_at_iteration(iters, total, target_iter=0.0, tol=1.0e-12):
    for i, it in enumerate(iters):
        if abs(it - target_iter) <= tol:
            return total[i]
    return total[0]


def main():
    parser = argparse.ArgumentParser(
        description="Plot optimization data for unsteady_mixer."
    )
    parser.add_argument(
        "--data",
        "--csv",
        dest="data",
        type=Path,
        default=Path("optimization_data.csv"),
        help="Path to optimization data file (CSV or TXT).",
    )
    parser.add_argument(
        "--case",
        type=Path,
        default=None,
        help="Path to .case file (defaults to the only .case in the folder).",
    )
    parser.add_argument(
        "--out",
        type=Path,
        default=Path("optimization_plot.png"),
        help="Output plot filename.",
    )
    parser.add_argument(
        "--show",
        action="store_true",
        help="Display the plot interactively.",
    )

    args = parser.parse_args()
    directory = Path.cwd()

    data_path = _resolve_data_path(directory, args.data)

    case_path = args.case
    if case_path is None:
        case_path = _find_case_file(directory)
    elif not case_path.is_absolute():
        case_path = directory / case_path

    iters, total, ratio, scaled_headers, scaled_series = _load_csv_data(data_path)
    norm_total = _normalization_from_total_at_iteration(iters, total, target_iter=0.0)
    if norm_total == 0.0:
        raise ValueError("Total objective at iteration 0 is zero; cannot normalize.")

    target_value = _load_case_target(case_path)

    fig, (ax_left, ax_right) = plt.subplots(1, 2, figsize=(13, 5))

    ax_left.plot(
        iters,
        [v / norm_total for v in total],
        label="Total objective function",
    )
    for label, series in zip(scaled_headers, scaled_series):
        ax_left.plot(iters, [v / norm_total for v in series], label=label)
    ax_left.set_xlabel("Iteration")
    ax_left.set_ylabel("Value / Total objective(iter=0)")
    ax_left.grid(True, alpha=0.3)
    ax_left.legend()

    ax_right.plot(iters, ratio, label="Target Dissipation.ratio")
    ax_right.axhline(
        target_value, color="black", linestyle="--", linewidth=1.0,
        label=f"Target = {target_value}"
    )
    ax_right.set_xlabel("Iteration")
    ax_right.set_ylabel("Ratio")
    ax_right.grid(True, alpha=0.3)
    ax_right.legend()

    fig.tight_layout()

    fig.savefig(args.out, dpi=150)
    if args.show:
        plt.show()


if __name__ == "__main__":
    main()
