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

Usage: compare_history.py --reference <ref.csv> [--rtol R] <current.csv>
"""

import argparse
import csv
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


def read_history(path):
    with open(path, newline="") as f:
        rows = [[c.strip() for c in row] for row in csv.reader(f) if row]
    header, data = rows[0], rows[1:]
    columns = {name: [float(r[i]) for r in data]
               for i, name in enumerate(header)}
    return header, columns


def is_exact(name):
    return name in EXACT_COLUMNS or name.endswith(".weight")


def main():
    parser = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    parser.add_argument("--reference", required=True)
    parser.add_argument("--rtol", type=float, default=DEFAULT_RTOL)
    parser.add_argument("--floor", type=float, default=DEFAULT_FLOOR)
    parser.add_argument("current")
    args = parser.parse_args()

    for path in (args.reference, args.current):
        if not os.path.isfile(path):
            sys.exit(f"FAIL: {path} not found.")
    ref_header, ref = read_history(args.reference)
    cur_header, cur = read_history(args.current)
    failures = []

    if cur_header != ref_header:
        sys.exit(f"FAIL: columns differ.\n  current:   {cur_header}\n"
                 f"  reference: {ref_header}")
    n_ref, n_cur = len(ref["iter"]), len(cur["iter"])
    if n_cur != n_ref:
        sys.exit(f"FAIL: {n_cur} iterations written, reference has {n_ref}.")

    # Guard against a comparison that cannot fail for want of content: the
    # reference must contain the q switch and a design that moved.
    if len(set(ref["q"])) < 2:
        sys.exit("FAIL: the reference never switches q.")
    if max(ref["max_design_change"]) <= 0.0:
        sys.exit("FAIL: the design never changes in the reference.")

    print(f"Comparing {n_cur} iterations, rtol = {args.rtol:.1e}, "
          f"floor = {args.floor:.1e}")
    print(f"{'column':<32} {'max |ref|':>12} {'max rel.diff':>12} "
          f"{'at iter':>7}")
    for name in ref_header:
        r, c = ref[name], cur[name]
        col_max = max(abs(v) for v in r)
        worst, worst_iter = 0.0, int(ref["iter"][0])
        for i in range(n_ref):
            if not math.isfinite(c[i]):
                failures.append(f"{name} at iter {int(ref['iter'][i])}: "
                                f"non-finite value {c[i]}")
                continue
            if is_exact(name):
                tol = 0.0
                scale = max(abs(r[i]), 1.0)
            else:
                scale = max(abs(r[i]), args.floor * col_max)
                tol = args.rtol * scale
            diff = abs(c[i] - r[i])
            rel = diff / scale if scale > 0.0 else diff
            if rel > worst:
                worst, worst_iter = rel, int(ref["iter"][i])
            if diff > tol:
                failures.append(
                    f"{name} at iter {int(ref['iter'][i])}: current "
                    f"{c[i]:.15e}, reference {r[i]:.15e}, "
                    f"rel.diff {rel:.3e}")
        print(f"{name:<32} {col_max:12.4e} {worst:12.3e} {worst_iter:7d}")

    if failures:
        print(f"\nFAIL: {len(failures)} value(s) outside the tolerance:")
        for f in failures[:40]:
            print("  " + f)
        sys.exit(1)
    print("\nPASS: history matches the reference.")


if __name__ == "__main__":
    main()
