#!/usr/bin/env python3
"""Compare two runs written by beam_neko / beam_ref: <tag>.csv and <tag>.x"""
import csv, sys
import numpy as np

N = 221184


def load(tag):
    rows = list(csv.DictReader(open(tag + '.csv')))
    X = np.fromfile(tag + '.x', dtype=np.float64).reshape(-1, N)
    return rows, X


def compare(a, b, label):
    ra, Xa = load(a)
    rb, Xb = load(b)
    print(f'== {label}: {a} vs {b}')
    print(' it   f0[a]            f0[b]            |df0|/f0  max|g_a-g_b|  max|x_a-x_b|')
    for it in range(1, min(len(Xa), len(Xb)) + 1):
        fa, fb = float(ra[it]['f0']), float(rb[it]['f0'])
        dg = max(abs(float(ra[it][f'g{i}']) - float(rb[it][f'g{i}'])) for i in range(1, 12))
        dx = np.max(np.abs(Xa[it - 1] - Xb[it - 1]))
        print(f'{it:3d}  {fa:.12f}  {fb:.12f}  {abs(fa - fb) / fb:8.2e}  {dg:11.2e}  {dx:11.2e}')


if __name__ == '__main__':
    compare(sys.argv[1], sys.argv[2], sys.argv[3] if len(sys.argv) > 3 else '')
