#!/usr/bin/env python3
"""Compare two runs (<tag>.csv + <tag>.x): iterates, constraints, x and the KKT
columns. Prints per-iteration differences and a one-line summary."""
import csv
import sys

import numpy as np

N = 221184


def load(tag):
    rows = list(csv.DictReader(open(tag + '.csv')))
    X = np.fromfile(tag + '.x', dtype=np.float64).reshape(-1, N)
    return rows, X


def rel(a, b):
    return abs(a - b) / max(abs(b), 1e-300)


def compare(a, b, label, quiet=False):
    ra, Xa = load(a)
    rb, Xb = load(b)
    nit = min(len(Xa), len(Xb))
    worst = dict(f0=0.0, g=0.0, x=0.0, kmax=0.0, knorm=0.0)
    if not quiet:
        print(f'== {label}: {a} vs {b}')
        print(' it  f0[a]           rel df0   max|dg|   max|dx|   '
              'kktmax[a]  kktmax[b]  rel dkmax  rel dknorm')
    for it in range(1, nit + 1):
        fa, fb = float(ra[it]['f0']), float(rb[it]['f0'])
        dg = max(abs(float(ra[it][f'g{i}']) - float(rb[it][f'g{i}']))
                 for i in range(1, 12))
        dx = float(np.max(np.abs(Xa[it - 1] - Xb[it - 1])))
        ka, kb = float(ra[it]['kktmax']), float(rb[it]['kktmax'])
        na, nb = float(ra[it]['kktnorm']), float(rb[it]['kktnorm'])
        row = dict(f0=rel(fa, fb), g=dg, x=dx, kmax=rel(ka, kb),
                   knorm=rel(na, nb))
        for k in worst:
            worst[k] = max(worst[k], row[k])
        if not quiet:
            print(f'{it:3d}  {fa:.12f}  {row["f0"]:8.2e}  {dg:8.2e}  '
                  f'{dx:8.2e}  {ka:9.3e}  {kb:9.3e}  {row["kmax"]:9.2e}  '
                  f'{row["knorm"]:9.2e}')
    print(f'SUMMARY {label}: max rel df0 {worst["f0"]:.2e}, max |dg| '
          f'{worst["g"]:.2e}, max |dx| {worst["x"]:.2e}, max rel dKKTmax '
          f'{worst["kmax"]:.2e}, max rel dKKTnorm {worst["knorm"]:.2e} '
          f'({nit} it)')


if __name__ == '__main__':
    quiet = '-q' in sys.argv
    args = [s for s in sys.argv[1:] if s != '-q']
    compare(args[0], args[1], args[2] if len(args) > 2 else '', quiet)
