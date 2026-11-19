"""Check a completed rugby-ball optimization against the paper.

Borrvall & Petersson (2003), Sec. 4.4, Table III. The run passes when

- it completed the expected number of MMA iterations,
- the final objective is within the given relative tolerance of the
  paper's value,
- the final design satisfies the volume constraint, and
- the paper's stopping rule, a relative change in the objective below
  the stopping tolerance in two consecutive iterations, is met within the
  final-q phase of the continuation before the last iteration. Neko-TOP
  runs a fixed number of iterations instead, so this checks that the run
  has converged in the paper's sense by the time it stops.

The history is the optimization_data.csv written by the run.
"""

import argparse
import csv
import sys

# Options as (name, type, default). Those without a default are required.
OPTIONS = (("--objective", float, None),
           ("--objective-rtol", float, 0.02),
           ("--iterations", int, None),
           ("--final-q", float, None),
           ("--stop-rtol", float, 1.0e-5))


def read_history(path):
    """Return the columns of a history file, keyed by their header."""
    with open(path, newline="") as handle:
        rows = [[cell.strip() for cell in row]
                for row in csv.reader(handle) if row]
    header, *data = rows
    return {name: [float(row[col]) for row in data]
            for col, name in enumerate(header)}


def parse_args():
    """Parse the command line."""
    parser = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    for name, kind, default in OPTIONS:
        parser.add_argument(name, type=kind, default=default,
                            required=default is None)
    parser.add_argument("history")
    return parser.parse_args()


def check_iterations(iters, args):
    """Check that every iteration up to the expected one was written."""
    passed = iters == list(range(args.iterations + 1))
    return passed, (f"iterations 0..{iters[-1]} written ({len(iters)} "
                    f"rows), expected 0..{args.iterations}")


def check_objective(obj, args):
    """Check the final objective against the paper's value."""
    rel = (obj[-1] - args.objective) / args.objective
    passed = abs(rel) <= args.objective_rtol
    return passed, (f"final objective {format(obj[-1], '.10g')} vs paper "
                    f"{args.objective}: {format(100 * rel, '+.3f')}% "
                    f"(limit +-{format(100 * args.objective_rtol, 'g')}%)")


def check_volume(vol):
    """Check that the final design satisfies the volume constraint."""
    return vol[-1] <= 0.0, (f"final volume constraint "
                            f"{format(vol[-1], '.3e')} <= 0")


def first_stop(obj, qval, args):
    """Return the first iteration at which the paper would stop, or None.

    The paper stops at iteration i when |J_k - J_{k-1}| / J_k is below the
    tolerance both for k equal to i - 1 and for k equal to i. Only
    iterations where both differences lie in the final-q phase count, so
    that the jump at the switch cannot take part.
    """
    for idx in range(2, len(obj)):
        if not all(qval[k] == args.final_q for k in (idx - 2, idx - 1, idx)):
            continue
        prev = abs(obj[idx - 1] - obj[idx - 2]) / abs(obj[idx - 1])
        last = abs(obj[idx] - obj[idx - 1]) / abs(obj[idx])
        if prev < args.stop_rtol and last < args.stop_rtol:
            return idx
    return None


def check_stop(iters, obj, qval, args):
    """Check that the paper's stopping rule is met before the last step."""
    rule = (f"stopping rule (relative change < "
            f"{format(args.stop_rtol, 'g')} twice, q"
            f" = {format(args.final_q, 'g')})")
    stop = first_stop(obj, qval, args)
    if stop is None:
        return False, f"{rule} never met"
    return stop < iters[-1], (f"{rule} first met at iteration "
                              f"{iters[stop]} of {iters[-1]}, J"
                              f" = {format(obj[stop], '.10g')}")


def main():
    """Run every check and exit with 1 if any of them fails."""
    args = parse_args()
    hist = read_history(args.history)
    iters = [int(val) for val in hist["iter"]]
    obj = hist["Total objective function"]

    failed = False
    for check in (lambda: check_iterations(iters, args),
                  lambda: check_objective(obj, args),
                  lambda: check_volume(hist["Volume constraint"]),
                  lambda: check_stop(iters, obj, hist["q"], args)):
        passed, message = check()
        failed = failed or not passed
        print(("PASS: " if passed else "FAIL: ") + message)

    if failed:
        sys.exit(1)
    print("Validation against the paper passed.")


if __name__ == "__main__":
    main()
