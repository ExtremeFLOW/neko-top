"""Check a completed rugby-ball optimization against the paper.

Borrvall & Petersson (2003), Sec. 4.4, Table III. The run passes when

- it completed the expected number of MMA iterations,
- the final objective is within --objective-rtol of the paper's value,
- the final design satisfies the volume constraint, and
- the paper's stopping rule, a relative change in the objective below
  --stop-rtol in two consecutive iterations, is met within the final-q
  phase of the continuation before the last iteration. Neko-TOP runs a
  fixed number of iterations instead, so this checks that the run has
  converged in the paper's sense by the time it stops.

Usage: check_paper.py --objective J --iterations N --final-q Q
           [--objective-rtol R] [--stop-rtol S] <optimization_data.csv>
"""

import argparse
import csv
import sys


def read_history(path):
    with open(path, newline="") as f:
        rows = [[c.strip() for c in row] for row in csv.reader(f) if row]
    header, data = rows[0], rows[1:]
    return {name: [float(r[i]) for r in data]
            for i, name in enumerate(header)}


def main():
    parser = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    parser.add_argument("--objective", type=float, required=True)
    parser.add_argument("--objective-rtol", type=float, default=0.02)
    parser.add_argument("--iterations", type=int, required=True)
    parser.add_argument("--final-q", type=float, required=True)
    parser.add_argument("--stop-rtol", type=float, default=1.0e-5)
    parser.add_argument("history")
    args = parser.parse_args()

    h = read_history(args.history)
    it = [int(v) for v in h["iter"]]
    J = h["Total objective function"]
    g = h["Volume constraint"]
    q = h["q"]
    ok = True

    def report(passed, message):
        nonlocal ok
        ok = ok and passed
        print(("PASS: " if passed else "FAIL: ") + message)

    report(it == list(range(args.iterations + 1)),
           f"iterations 0..{it[-1]} written ({len(it)} rows), expected "
           f"0..{args.iterations}")

    rel = (J[-1] - args.objective) / args.objective
    report(abs(rel) <= args.objective_rtol,
           f"final objective {J[-1]:.10g} vs paper {args.objective}: "
           f"{100 * rel:+.3f}% (limit +-{100 * args.objective_rtol:g}%)")

    report(g[-1] <= 0.0, f"final volume constraint {g[-1]:.3e} <= 0")

    # The paper stops at iteration i when |J_k - J_{k-1}| / J_k is below
    # the tolerance for k = i - 1 and k = i. Only count iterations where
    # both differences lie in the final-q phase, so that the jump at the
    # switch cannot take part.
    stop = None
    for i in range(2, len(J)):
        if not all(q[k] == args.final_q for k in (i - 2, i - 1, i)):
            continue
        d1 = abs(J[i - 1] - J[i - 2]) / abs(J[i - 1])
        d2 = abs(J[i] - J[i - 1]) / abs(J[i])
        if d1 < args.stop_rtol and d2 < args.stop_rtol:
            stop = i
            break
    if stop is None:
        report(False, f"stopping rule (relative change < {args.stop_rtol:g}"
               f" twice, q = {args.final_q:g}) never met")
    else:
        report(stop < it[-1],
               f"stopping rule (relative change < {args.stop_rtol:g} "
               f"twice, q = {args.final_q:g}) first met at iteration "
               f"{it[stop]} of {it[-1]}, J = {J[stop]:.10g}")

    if not ok:
        sys.exit(1)
    print("Validation against the paper passed.")


if __name__ == "__main__":
    main()
