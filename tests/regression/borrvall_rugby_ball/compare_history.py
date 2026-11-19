"""Compare an optimization history against a reference history.

Every column of optimization_data.csv is compared, row by row: the
objective and its terms, the volume constraint and volume, the KKT
residuals, the design change and q. Integer-like columns (iteration,
weights, scaling factor, q) must match exactly. A real-valued entry
passes when

    |current - reference| <= rtol * max(|reference|, floor * column_max)

where column_max is the largest magnitude in that reference column. The
floor stops entries that are round-off zeros in the reference, such as
the volume constraint of a design started exactly on the limit, from
requiring an exact match.

The reference history is given with the reference option and the current
one as the positional argument.
"""

import argparse
import csv
import itertools
import math
import os
import sys

# Measured on regression_0p8.case (GNU, Release, CPU): repeated runs on
# two ranks are bit-identical, and two against four ranks differ by at most
# 1.1e-9 in any column (3.9e-12 in the objective). Lowering f_max by 4%
# gives a relative difference above 5e-2 somewhere in every real-valued
# column. The tolerance sits about
# 100 times above the rank-count spread and far below any real change.
DEFAULT_RTOL = 1.0e-7
DEFAULT_FLOOR = 1.0e-3

EXACT_COLUMNS = ("iter", "scaling factor", "q")

# Number of failing values listed in full.
MAX_LISTED = 40
INDENT = " " * 2


def read_history(path):
    """Return the header of a history file and its columns by name."""
    with open(path, newline="") as handle:
        rows = [[cell.strip() for cell in row]
                for row in csv.reader(handle) if row]
    header, *data = rows
    columns = {name: [float(row[col]) for row in data]
               for col, name in enumerate(header)}
    return header, columns


def is_exact(name):
    """Return whether a column must match the reference exactly."""
    return name in EXACT_COLUMNS or name.endswith(".weight")


def parse_args():
    """Parse the command line."""
    parser = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    parser.add_argument("--reference", required=True)
    parser.add_argument("--rtol", type=float, default=DEFAULT_RTOL)
    parser.add_argument("--floor", type=float, default=DEFAULT_FLOOR)
    parser.add_argument("current")
    return parser.parse_args()


def check_shapes(ref_header, ref, cur_header, cur):
    """Exit when the histories cannot be compared, or the reference is
    too poor for the comparison to be able to fail."""
    if cur_header != ref_header:
        sys.exit(f"FAIL: columns differ.\n{INDENT}current: {' ' * 2}"
                 f"{cur_header}\n{INDENT}reference: {ref_header}")
    n_ref, n_cur = len(ref["iter"]), len(cur["iter"])
    if n_cur != n_ref:
        sys.exit(f"FAIL: {n_cur} iterations written, reference has {n_ref}.")

    # Guard against a comparison that cannot fail for want of content: the
    # reference must contain the q switch and a design that moved.
    if len(set(ref["q"])) < 2:
        sys.exit("FAIL: the reference never switches q.")
    if max(ref["max_design_change"]) <= 0.0:
        sys.exit("FAIL: the design never changes in the reference.")


def scale_and_tol(name, rval, col_max, args):
    """Return the scale of a reference value and the absolute tolerance."""
    if is_exact(name):
        return max(abs(rval), 1.0), 0.0
    scale = max(abs(rval), args.floor * col_max)
    return scale, args.rtol * scale


def compare_column(name, ref, cur, args):
    """Compare one column, row by row.

    Return the largest relative difference, the iteration where it occurs
    and the description of every value outside the tolerance.
    """
    col_max = max(abs(val) for val in ref[name])
    worst, worst_iter = 0.0, int(ref["iter"][0])
    failures = []
    for row, (rval, cval) in enumerate(zip(ref[name], cur[name])):
        where = f"{name} at iter {int(ref['iter'][row])}"
        if not math.isfinite(cval):
            failures.append(f"{where}: non-finite value {cval}")
            continue
        scale, tol = scale_and_tol(name, rval, col_max, args)
        diff = abs(cval - rval)
        rel = diff / scale if scale > 0.0 else diff
        if rel > worst:
            worst, worst_iter = rel, int(ref["iter"][row])
        if diff > tol:
            failures.append(f"{where}: current {format(cval, '.15e')}, "
                            f"reference {format(rval, '.15e')}, "
                            f"rel.diff {format(rel, '.3e')}")
    print(f"{name.ljust(32)} {format(col_max, '12.4e')} "
          f"{format(worst, '12.3e')} {format(worst_iter, '7d')}")
    return failures


def main():
    """Compare the histories and exit with 1 if any value differs."""
    args = parse_args()
    for path in (args.reference, args.current):
        if not os.path.isfile(path):
            sys.exit(f"FAIL: {path} not found.")
    ref_header, ref = read_history(args.reference)
    cur_header, cur = read_history(args.current)
    check_shapes(ref_header, ref, cur_header, cur)

    print(f"Comparing {len(cur['iter'])} iterations, "
          f"rtol = {format(args.rtol, '.1e')}, "
          f"floor = {format(args.floor, '.1e')}")
    print(f"{'column'.ljust(32)} {'max |ref|'.rjust(12)} "
          f"{'max rel.diff'.rjust(12)} {'at iter'.rjust(7)}")
    failures = []
    for name in ref_header:
        failures.extend(compare_column(name, ref, cur, args))

    if failures:
        print(f"\nFAIL: {len(failures)} value(s) outside the tolerance: "
              .rstrip())
        for failure in itertools.islice(failures, MAX_LISTED):
            print(INDENT + failure)
        sys.exit(1)
    print("\nPASS: history matches the reference.")


if __name__ == "__main__":
    main()
