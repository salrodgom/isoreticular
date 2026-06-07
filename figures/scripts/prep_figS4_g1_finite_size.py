#!/usr/bin/env python3
"""
Tabulate the G_1 finite-size benchmark of Figure S4 of the manuscript:
D8R distortion Delta vs hydrostatic pressure p, computed in two cell sizes:
the 1x1x1 cubic unit cell (48 T-atoms, red, dense pressure scan) and the
2x2x2 supercell (384 T-atoms, blue, coarser pressure scan).

The data are those previously hard-coded in the pgfplots/TikZ block of
main.tex; they are now consolidated here so manuscript_figures/ is the
reproducible source of truth.

Reads nothing; writes:
    data/figS4_g1_finite_size_1x1x1_stable.dat    columns: p [GPa], Delta [A]
        (1x1x1 points with phonon-stable parent branch, omega_1 >= 0)
    data/figS4_g1_finite_size_1x1x1_unstable.dat  columns: p [GPa], Delta [A]
        (1x1x1 points where omega_1 < 0; metastable cubic basin)
    data/figS4_g1_finite_size_2x2x2.dat           columns: p [GPa], Delta [A]
        (2x2x2 supercell equilibrium points)
"""
import os
import numpy as np


# 1x1x1 dense pressure scan with omega_1 >= 0 (filled red dots).
G1_1x1x1_STABLE = [
    (0.000, 0.000), (0.100, 0.000), (0.200, 0.000), (0.300, 0.000),
    (0.400, 0.000), (0.500, 0.017), (0.600, 0.000), (0.700, 0.000),
    (0.800, 0.000), (0.805, 0.000), (0.805, 0.056), (0.810, 0.000),
    (0.810, 0.084), (0.815, 0.000), (0.815, 0.221), (0.820, 0.000),
    (0.820, 0.220), (0.825, 0.000), (0.825, 0.219), (0.830, 0.000),
    (0.830, 0.247), (0.835, 0.000), (0.835, 0.248), (0.840, 0.000),
    (0.840, 0.246), (0.845, 0.000), (0.845, 0.246), (0.850, 0.000),
    (0.850, 0.263), (0.855, 0.000), (0.855, 0.264), (0.860, 0.000),
    (0.860, 0.301), (0.865, 0.000), (0.865, 0.259), (0.870, 0.000),
    (0.870, 0.298), (0.8725, 0.001), (0.8725, 0.298), (0.875, 0.000),
    (0.875, 0.296), (0.8775, 0.180), (0.878, 0.001), (0.878, 0.295),
    (0.879, 0.001), (0.879, 0.296), (0.880, 0.296), (0.881, 0.001),
    (0.881, 0.295), (0.882, 0.001), (0.882, 0.357), (0.8825, 0.001),
    (0.8825, 0.357), (0.885, 0.002), (0.885, 0.296), (0.8875, 0.002),
    (0.8875, 0.356), (0.890, 0.002), (0.890, 0.356), (0.895, 0.001),
    (0.895, 0.358), (0.900, 0.396), (0.905, 0.000), (0.905, 0.395),
    (0.910, 0.000), (0.910, 0.551), (0.920, 0.000), (0.920, 0.529),
    (0.930, 0.001), (0.930, 0.547), (0.940, 0.002), (0.940, 0.545),
    (0.950, 0.327), (0.960, 0.437), (0.970, 0.523), (0.980, 0.598),
    (0.990, 0.662), (1.000, 0.672), (1.100, 1.111), (1.200, 1.360),
    (1.300, 1.553), (1.400, 1.710), (1.500, 1.843), (1.600, 1.960),
    (1.700, 2.065), (1.800, 2.159), (1.900, 2.244), (2.000, 2.321),
]

# 1x1x1 points where omega_1 < 0 (metastable Im-3m, open red circles).
G1_1x1x1_UNSTABLE = [
    (0.8775, 0.000), (0.880, 0.000), (0.950, 0.008), (0.960, 0.012),
    (0.970, 0.013), (0.980, 0.015), (0.990, 0.017),
]

# 2x2x2 supercell equilibrium points (open blue squares).
G1_2x2x2 = [
    (0.000, 0.000), (0.100, 0.000), (0.200, 0.000), (0.300, 0.000),
    (0.400, 0.000), (0.500, 0.000), (0.600, 0.000), (0.700, 0.000),
    (0.800, 0.012), (0.850, 0.000), (0.860, 0.000), (0.875, 0.000),
    (0.900, 0.000), (0.910, 0.000), (0.925, 0.000), (0.930, 0.000),
    (0.940, 0.000), (0.950, 0.757), (0.968, 0.757), (0.972, 0.757),
    (0.990, 0.757), (1.000, 0.757), (1.050, 0.878), (1.100, 1.111),
    (1.150, 1.246), (1.158, 1.247), (1.200, 1.360), (1.250, 1.463),
    (1.300, 1.553), (1.350, 1.635), (1.400, 1.711), (1.450, 1.779),
    (1.500, 1.844), (1.600, 1.961), (1.700, 2.065), (1.800, 2.159),
    (1.900, 2.251), (2.000, 2.321),
]


def main():
    here = os.path.dirname(os.path.abspath(__file__))
    out_dir = os.path.normpath(os.path.join(here, '..', 'data'))
    os.makedirs(out_dir, exist_ok=True)

    cases = [
        ('figS4_g1_finite_size_1x1x1_stable.dat',
         G1_1x1x1_STABLE,
         '1x1x1 cubic unit cell, phonon-stable parent branch (omega_1 >= 0).'),
        ('figS4_g1_finite_size_1x1x1_unstable.dat',
         G1_1x1x1_UNSTABLE,
         '1x1x1 cubic unit cell, metastable Im-3m branch (omega_1 < 0).'),
        ('figS4_g1_finite_size_2x2x2.dat',
         G1_2x2x2,
         '2x2x2 supercell equilibrium points.'),
    ]
    for fname, rows, header in cases:
        np.savetxt(os.path.join(out_dir, fname), np.asarray(rows),
                   header=f'{header}\nColumns: p [GPa]   Delta [A]',
                   fmt='%.6f')
        print(f'Wrote {fname}  ({len(rows)} rows)')


if __name__ == '__main__':
    main()
