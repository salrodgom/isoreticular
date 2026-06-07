#!/usr/bin/env python3
"""
Tabulate U_eff/A = 1/delta^2 for the isoreticular family G_1..G_5,
Figure S3 (Landau-tricriticality indicator) of the manuscript.

Reads nothing; writes:
    data/figS3_landau_ueff.dat   columns: k-1, U_eff/A [A^-2]

Values taken from Table tab:parameters (delta column) of the manuscript,
using the soft-mode p_c fits (Cowley-Levanyuk extrapolation):
        delta [A]
    G_1   2.210
    G_2   1.236
    G_3   1.014
    G_4   0.744
    G_5   0.464

Run from the repo root:
    python3 scripts/prep_figS3_landau_ueff.py
"""
import os
import numpy as np

DELTA = {1: 2.210, 2: 1.236, 3: 1.014, 4: 0.744, 5: 0.464}


def main():
    here = os.path.dirname(os.path.abspath(__file__))
    out_dir = os.path.normpath(os.path.join(here, '..', 'data'))
    os.makedirs(out_dir, exist_ok=True)

    rows = [(k - 1, 1.0 / DELTA[k] ** 2) for k in sorted(DELTA)]
    arr = np.asarray(rows)
    header = ("u_eff/alpha = 1/delta^2 for the isoreticular RHO family.\n"
              "Columns: k-1  u_eff/alpha [A^-2]")
    np.savetxt(os.path.join(out_dir, 'figS3_landau_ueff.dat'),
               arr, header=header, fmt='%.6f')
    print("Wrote", os.path.join(out_dir, 'figS3_landau_ueff.dat'))
    for k, (km1, uoa) in zip(sorted(DELTA), rows):
        print(f"  G_{k} (k-1={km1:.0f}):  u_eff/alpha = {uoa:.4f} A^-2")


if __name__ == '__main__':
    main()
