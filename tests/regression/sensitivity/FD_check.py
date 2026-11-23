#!/usr/bin/env python3
"""Compare sensitivity `FD_check_*.csv` output against `reference_data/`.

Two modes, sharing one implementation of the comparison so the tolerances
cannot drift apart:

* **Sweep mode** (no ``--case``) -- glob every ``FD_check*.csv`` in the
  working directory, plot each against its reference and compare the ones
  that have a reference. This is what the manual ``run.sh`` drives. A case
  without a reference is plotted and skipped, because ``cases/`` holds more
  case files than are registered as ctest cases.

* **Gate mode** (``--case NAME``) -- examine exactly one CSV, the one the
  named case just wrote, and treat a **missing reference as a failure**.
  This is what ``run_sensitivity_regression_check.sh`` runs after its
  ``mpirun``, so that every registered ctest case is really gated on the
  recorded solution. The scoping matters: the registered cases share a
  single working directory (they hold ``RESOURCE_LOCK
  "sensitivity_regression_workspace"``), so CSVs from earlier cases are
  still lying around and an unscoped glob would make every case re-check
  every other case's output.

Tolerances (both relative, per row):

* ``VALUE_TOL`` -- ``F`` and ``dFdx``. Measured against the references at
  the time of writing, ``F`` agrees to <= 1.4e-14 and ``dFdx`` to
  <= 5.6e-7, the latter on ``viscous_dissipation`` whose reference was
  generated at a different rank count from the one the ctest lane uses.
  5e-6 sits near the geometric middle of that 5.6e-7 offset and the 7.9e-5
  shift a flipped ``case.numerics.dealias`` / ``optimization.design.dealias``
  produces: ~9x of headroom above the known offset, ~16x of margin below
  the smallest regression it must catch.

* ``ERROR_TOL`` -- the ``error`` column, which the comparison used to
  ignore entirely. It is a difference of differences (the relative gap
  between the finite-difference estimate and the analytic sensitivity) and
  so is far noisier than ``F``: the same sweep shows offsets up to 1.7e-3
  against references generated under a different configuration, where
  ``F`` was bit-identical. Repeat runs of one case under one configuration
  are bit-identical in every column, so the residual risk is a
  configuration change, not jitter. 1e-2 clears the worst observed offset
  by ~6x while still catching a flipped ``dealias_forcing`` (-3.79e-5 ->
  -6.04e-5, 59%) or ``dealias_sensitivity`` (-> -4.23e-5, 12%) by 12x or
  better -- and those two flags move **only** this column, leaving ``F``
  and ``dFdx`` bit-identical, so without this check they are invisible to
  both this comparison and the driver's own assertion.
"""

import argparse
import csv
import glob
import math
import os
import sys

#: Relative tolerance for the `F` and `dFdx` columns.
VALUE_TOL = 5e-6

#: Relative tolerance for the (much noisier) `error` column.
ERROR_TOL = 1e-2

#: Columns compared, and the tolerance each is held to.
TOLERANCES = (("F", VALUE_TOL), ("dFdx", VALUE_TOL), ("error", ERROR_TOL))

#: Columns the driver writes, in order.
COLUMNS = ["perturbation", "F", "dFdx", "error"]

#: Header line the driver writes at the start of each block.
HEADER = ",".join(COLUMNS)


def read_csv(path):
    """Read a driver-written CSV into {column: [values]}.

    Only the most recent block is kept. The driver logs through Neko's
    `csv_file_t`, whose `overwrite` flag defaults to false, so it opens
    with `position="append"`. Re-running a case in a directory that
    already holds its CSV therefore *appends* a second header-plus-rows
    block rather than replacing the file, and parsing the whole file would
    feed the repeated header rows in as data. Take the last block, which
    is the run that just happened.

    Deliberately stdlib-only: this runs inside a ctest gate, and the gate
    should not be able to fail (or, worse, be skipped) because a plotting
    or dataframe package is missing on the machine running it.
    """
    with open(path, newline="") as handle:
        rows = [row for row in csv.reader(handle) if row and any(row)]

    starts = [i for i, row in enumerate(rows) if row == COLUMNS]
    if not starts:
        raise SystemExit(f"{path}: no '{HEADER}' header found.")

    block = rows[starts[-1] + 1:]
    try:
        return {name: [float(row[i]) for row in block]
                for i, name in enumerate(COLUMNS)}
    except (ValueError, IndexError) as error:
        raise SystemExit(f"{path}: malformed data row ({error}).")


def compare(csv_file, data, reference_data):
    """Compare one CSV against its reference; return a list of failures."""
    failures = []

    n_rows = len(data["perturbation"])
    n_rows_ref = len(reference_data["perturbation"])
    if n_rows != n_rows_ref:
        return [
            f"{csv_file}: row count {n_rows} does not match reference "
            f"{n_rows_ref} -- the case or its perturbation list changed."
        ]

    for column, tol in TOLERANCES:
        current = data[column]
        reference = reference_data[column]

        max_rel_error = max(
            abs((cur - ref) / (ref + 1e-30))
            for cur, ref in zip(current, reference))
        if max_rel_error > tol:
            failures.append(
                f"{csv_file}: {column} does not match reference "
                f"(max relative error: {max_rel_error:.6e}, "
                f"tolerance: {tol:.1e})\n"
                f"  current:   {current}\n"
                f"  reference: {reference}"
            )

    return failures


def _abs(values):
    """Element-wise absolute value of a plain list of floats."""
    return [abs(value) for value in values]


def reference_path(csv_file, reference_dir):
    """`FD_check_foo.csv` -> `<reference_dir>/ref_FD_check_foo.csv`."""
    ref_name = os.path.basename(csv_file).replace("FD_check", "ref_FD_check",
                                                  1)
    return os.path.join(reference_dir, ref_name)


def plot(csv_files, frames, references):
    """Plot every CSV against its reference; only used in sweep mode.

    matplotlib is imported here rather than at module scope so the ctest
    gate (which always passes --no-plot) never needs it installed.
    """
    import matplotlib.pyplot as plt

    n_files = len(csv_files)
    ncols = math.ceil(math.sqrt(n_files))
    nrows = math.ceil(n_files / ncols)

    fig, axes = plt.subplots(nrows,
                             ncols,
                             figsize=(5 * ncols, 4 * nrows),
                             squeeze=False)
    axes = axes.ravel()

    for i, csv_file in enumerate(csv_files):
        data = frames[csv_file]
        axes[i].loglog(_abs(data["perturbation"]),
                       _abs(data["error"]),
                       marker="o",
                       linestyle="-",
                       label="Current")

        data_ref = references.get(csv_file)
        if data_ref is not None:
            axes[i].loglog(_abs(data_ref["perturbation"]),
                           _abs(data_ref["error"]),
                           linestyle="--",
                           label="Reference")

        name = csv_file.split("FD_check_")[-1].replace(".csv", "")
        axes[i].set_title(name)
        axes[i].set_ylabel("Error")
        axes[i].grid(True, which="both", linestyle="--", linewidth=0.5)
        axes[i].legend()

    for ax in axes[n_files:]:
        ax.axis("off")

    plt.tight_layout()
    if not os.path.exists("plots"):
        os.makedirs("plots")
    plt.savefig("plots/FD_comparison.png", dpi=200)
    plt.show()


def parse_args(argv):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--case",
        metavar="NAME",
        help="Check only FD_check_<NAME>.csv, and fail if it has no "
        "reference. Registered ctest cases share one working directory, "
        "so the check must be scoped to the case that just ran.")
    parser.add_argument(
        "--reference-dir",
        default="reference_data",
        help="Directory holding the ref_FD_check_*.csv files "
        "(default: %(default)s, relative to the working directory).")
    parser.add_argument("--no-plot",
                        action="store_true",
                        help="Skip plotting entirely, so no matplotlib "
                        "import, no plots/ directory and no window. Used "
                        "by the ctest lane.")
    return parser.parse_args(argv)


def main(argv=None):
    args = parse_args(argv)

    if args.case:
        csv_files = [f"FD_check_{args.case}.csv"]
        if not os.path.exists(csv_files[0]):
            raise SystemExit(
                f"{csv_files[0]} not found in {os.getcwd()} -- the driver "
                "did not write the expected output for this case.")
    else:
        csv_files = sorted(glob.glob("FD_check*.csv"))
        if not csv_files:
            raise SystemExit("No files matching 'FD_check*.csv' found.")

    frames = {}
    references = {}
    failures = []

    for csv_file in csv_files:
        frames[csv_file] = read_csv(csv_file)

        ref_file = reference_path(csv_file, args.reference_dir)
        if not os.path.exists(ref_file):
            if args.case:
                # A silently-skipped comparison is how this lane ended up
                # with reference files that gated nothing at all; in gate
                # mode an absent reference is a failure, not a pass.
                failures.append(
                    f"{csv_file}: no reference data at {ref_file}. Every "
                    "registered case must be gated on a recorded solution; "
                    "generate one (see reference_data/) rather than "
                    "leaving the case unchecked.")
            continue

        references[csv_file] = read_csv(ref_file)
        failures.extend(compare(csv_file, frames[csv_file],
                                references[csv_file]))

    if not args.no_plot:
        plot(csv_files, frames, references)

    if failures:
        for failure in failures:
            print(f"Error: {failure}")
        raise SystemExit("Discrepancies found compared to reference data.")

    for csv_file in csv_files:
        if csv_file in references:
            ref_file = reference_path(csv_file, args.reference_dir)
            print(f"{csv_file} matches {ref_file}.")


if __name__ == "__main__":
    main(sys.argv[1:])
