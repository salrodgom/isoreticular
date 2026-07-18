#!/usr/bin/env python3
"""
Automated extraction of free-energy barriers Delta G* from PLUMED reweighted
FES in the isoRHO/ directory (the production set used to generate fig:FES of
the main text).

Each G_k subdirectory contains four 2D FES files (fes_1.dat .. fes_4.dat),
one per hydrostatic pressure shown in the four columns of fig:FES. The
mapping fes_<N>.dat -> P_GPa is fixed and read from the gp_multiplot script
of each member; it is hard-coded in PRESSURES below.

Free-energy file format (one line per (delta, cell) grid point):
    #! FIELDS delta cell file.free
    delta cell free[kJ/mol]
where delta is the D8R distortion (column 1, Angstrom, range [0, 2] or [0, 1.5])
and cell is the cubic cell parameter (column 2, Angstrom). This script reads
the columns in the order they appear in the file (FIXES the column-swap bug
of the original extract_barriers.py).

Barrier definition:
    1) Gaussian-smooth F(delta, cell) with sigma = 2 px.
    2) Subtract the global minimum so F >= 0.
    3) Identify the cubic basin (delta < delta_split) and the broken basin
       (delta > delta_split). Default delta_split = 0.3 Angstrom (separates
       Im-3m centric noise from the I-43m acentric region).
    4) Project F along the cell direction by taking min_cell F(delta, cell)
       at every delta. This is the minimum-energy path connecting the two
       basins through the saddle.
    5) Saddle = max of the projected F between delta_cub and delta_brk.
    6) Delta G* = F_saddle - min(F_cub, F_brk).

Output per file:
    P [GPa], p/p_c, N_T, delta_cub, delta_sad, delta_brk,
    barr [kJ/mol], barr [kBT/Si], barr [1e-2 kBT/Si]

The pressure of fig:FES used in tab:SI-fes-barriers is the one at which
the barrier is most clearly resolved among the four. A LaTeX summary block
at the end of the run echoes that row.
"""

import os
import re
import sys
import glob
import numpy as np

KBT_300K = 298.15 * 0.0083144626  # kJ/mol at the FES temperature (298.15 K)

# Number of T-atoms in the simulation cell used for the FES of each member.
# G_1 uses the 2x2x2 supercell (384 T); G_2-G_5 use their unit cells.
N_SI = {'G_1': 384, 'G_2': 240, 'G_3': 672, 'G_4': 1440, 'G_5': 2640}

# Canonical soft-mode critical pressures (tab:parameters of the manuscript).
P_C = {'G_1': 0.9418, 'G_2': 0.5338, 'G_3': 0.3627, 'G_4': 0.2241, 'G_5': 0.1282}

# Production directories with the 2D reweighted FES.
ISO_ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__),
                                        '..', 'isoRHO'))
DIRS = {f'G_{k}': os.path.join(ISO_ROOT, f'G{k}') for k in range(1, 6)}

# Mapping fes_<N>.dat -> pressure [GPa], extracted from the gp_multiplot
# 'splot' lines of each member.
PRESSURES = {
    'G_1': {'fes_1.dat': 1.000, 'fes_2.dat': 1.200, 'fes_3.dat': 1.225, 'fes_4.dat': 1.300},
    'G_2': {'fes_1.dat': 0.500, 'fes_2.dat': 0.750, 'fes_3.dat': 0.820, 'fes_4.dat': 0.900},
    'G_3': {'fes_1.dat': 0.100, 'fes_2.dat': 0.500, 'fes_3.dat': 0.700, 'fes_4.dat': 1.000},
    'G_4': {'fes_1.dat': 0.200, 'fes_2.dat': 0.500, 'fes_3.dat': 0.600, 'fes_4.dat': 0.800},
    'G_5': {'fes_1.dat': 0.300, 'fes_2.dat': 0.400, 'fes_3.dat': 0.500, 'fes_4.dat': 0.700},
}

DELTA_SPLIT = 0.30   # Angstrom; basins separator between Im-3m and I-43m


def load_fes_2d(fname):
    """Load a 2D FES file with columns (delta, cell, F) into arrays.

    Returns
    -------
    udelta : 1D ndarray
        Sorted unique delta values [Angstrom].
    ucell : 1D ndarray
        Sorted unique cell parameter values [Angstrom].
    F : 2D ndarray of shape (len(udelta), len(ucell))
        Free energy [kJ/mol] on the (delta, cell) grid.
    """
    deltas, cells, fs = [], [], []
    with open(fname) as fh:
        for line in fh:
            line = line.strip()
            if not line or line.startswith('#'):
                continue
            parts = line.split()
            try:
                d = float(parts[0])  # column 1 = delta (FIELDS delta cell file.free)
                c = float(parts[1])  # column 2 = cell
                e = float(parts[2])  # column 3 = free energy
            except (ValueError, IndexError):
                continue
            deltas.append(d)
            cells.append(c)
            fs.append(e)
    deltas = np.asarray(deltas)
    cells = np.asarray(cells)
    fs = np.asarray(fs)
    udelta = np.unique(deltas)
    ucell = np.unique(cells)
    F = np.full((len(udelta), len(ucell)), np.nan)
    # Map (delta, cell) -> index via sorted unique arrays. This is much
    # faster than a per-point mask lookup for the 201x201 grid.
    didx = {d: i for i, d in enumerate(udelta)}
    cidx = {c: j for j, c in enumerate(ucell)}
    for d, c, e in zip(deltas, cells, fs):
        F[didx[d], cidx[c]] = e
    return udelta, ucell, F


def gauss_smooth_2d(F, sigma=2.0):
    """Separable 2D Gaussian smoothing with edge replication."""
    Fp = np.where(np.isfinite(F), F, np.nanmax(F[np.isfinite(F)]))
    k = int(np.ceil(3 * sigma))
    x = np.arange(-k, k + 1)
    g = np.exp(-x ** 2 / (2 * sigma ** 2))
    g /= g.sum()
    F1 = np.zeros_like(Fp)
    for i in range(Fp.shape[0]):
        row = np.pad(Fp[i], k, mode='edge')
        F1[i] = np.convolve(row, g, mode='valid')
    F2 = np.zeros_like(F1)
    for j in range(F1.shape[1]):
        col = np.pad(F1[:, j], k, mode='edge')
        F2[:, j] = np.convolve(col, g, mode='valid')
    return F2


def find_basins_and_saddle(udelta, ucell, F, delta_split=DELTA_SPLIT):
    """Find cubic and broken basins and the saddle on the minimum-cell path.

    Returns None if either basin is empty, or if the broken basin sits at
    the upper grid boundary (suggests an unconverged FES on that side).
    """
    F = F - np.nanmin(F)
    mask_cub = udelta < delta_split
    mask_brk = udelta > delta_split
    if mask_cub.sum() < 2 or mask_brk.sum() < 2:
        return None
    # cubic basin: argmin of F over (delta < split, all cell)
    sub_cub = F[mask_cub]
    di_cub_rel, dj_cub = np.unravel_index(np.nanargmin(sub_cub), sub_cub.shape)
    di_cub = np.where(mask_cub)[0][di_cub_rel]
    cubic_F = F[di_cub, dj_cub]
    cubic_delta = udelta[di_cub]
    cubic_cell = ucell[dj_cub]
    # broken basin: argmin of F over (delta > split, all cell)
    sub_brk = F[mask_brk]
    di_brk_rel, dj_brk = np.unravel_index(np.nanargmin(sub_brk), sub_brk.shape)
    di_brk = np.where(mask_brk)[0][di_brk_rel]
    broken_F = F[di_brk, dj_brk]
    broken_delta = udelta[di_brk]
    broken_cell = ucell[dj_brk]
    if di_brk == len(udelta) - 1:
        return None  # broken basin at grid boundary, not a real minimum
    # Saddle: walk delta from cubic to broken, take min over cell at each delta,
    # then take the max of that 1D profile.
    lo, hi = min(di_cub, di_brk), max(di_cub, di_brk)
    proj = np.nanmin(F[lo:hi + 1], axis=1)
    saddle_local = int(np.nanargmax(proj))
    saddle_F = proj[saddle_local]
    saddle_delta = udelta[lo + saddle_local]
    return {
        'cubic_delta': cubic_delta, 'cubic_cell': cubic_cell, 'cubic_F': cubic_F,
        'broken_delta': broken_delta, 'broken_cell': broken_cell, 'broken_F': broken_F,
        'saddle_delta': saddle_delta, 'saddle_F': saddle_F,
        'barrier_fwd': saddle_F - cubic_F,
        'barrier_bwd': saddle_F - broken_F,
        'barrier_low': saddle_F - min(cubic_F, broken_F),
    }


# === main ============================================================
results = []
hdr = (f"{'G_k':4s} {'file':12s} {'P[GPa]':>7s} {'p/p_c':>6s} {'N_T':>6s} "
       f"{'2bas':4s} {'D_cub':>6s} {'D_sad':>6s} {'D_brk':>6s} "
       f"{'DG*[kJ/mol]':>12s} {'DG*[kBT/Si]':>12s} {'DG*[1e-2kBT/Si]':>15s}")
print(hdr)
print('-' * len(hdr))

for member, dirpath in DIRS.items():
    if not os.path.isdir(dirpath):
        print(f"{member:4s} (directory not found: {dirpath})")
        continue
    n_t = N_SI[member]
    pc = P_C[member]
    pmap = PRESSURES[member]
    fes_files = sorted(glob.glob(os.path.join(dirpath, 'fes_[1-9].dat')))
    if not fes_files:
        print(f"{member:4s} (no fes_<N>.dat files in {dirpath})")
        continue
    for fpath in fes_files:
        fname = os.path.basename(fpath)
        if fname not in pmap:
            continue
        p_gpa = pmap[fname]
        try:
            udelta, ucell, F = load_fes_2d(fpath)
            if F.size == 0:
                continue
            Fs = gauss_smooth_2d(F, sigma=2.0)
            r = find_basins_and_saddle(udelta, ucell, Fs)
        except Exception as exc:
            print(f"{member:4s} {fname:12s} ERROR: {exc}")
            continue
        if r is None:
            print(f"{member:4s} {fname:12s} {p_gpa:7.3f} {p_gpa/pc:6.2f} {n_t:>6d} "
                  f"{'no':4s}")
            continue
        barr = r['barrier_low']
        barr_si = barr / (KBT_300K * n_t)         # per-Si in units of kBT
        barr_intens = barr_si * 100               # in manuscript units of 1e-2 kBT/Si (matches fig:FES colour scale)
        print(f"{member:4s} {fname:12s} {p_gpa:7.3f} {p_gpa/pc:6.2f} {n_t:>6d} "
              f"{'yes':4s} {r['cubic_delta']:6.2f} {r['saddle_delta']:6.2f} "
              f"{r['broken_delta']:6.2f} {barr:12.2f} {barr_si:12.4f} {barr_intens:15.5f}")
        results.append(dict(
            member=member, fname=fname, p_gpa=p_gpa, pc=pc, n_t=n_t,
            barr=barr, barr_si=barr_si, barr_intens=barr_intens,
            cubic_delta=r['cubic_delta'], broken_delta=r['broken_delta'],
            saddle_delta=r['saddle_delta'],
        ))

# === LaTeX summary block: most resolved barrier per G_k ==============
print("\n\n=== Recommended row per G_k for tab:SI-fes-barriers ===")
print("Selection: pressure with the largest extracted intensive barrier among the four panels of fig:FES.")
by_member = {}
for r in results:
    m = r['member']
    if m not in by_member or r['barr_intens'] > by_member[m]['barr_intens']:
        by_member[m] = r

print("\n% Auto-generated by scripts_FES/extract_barriers_isoRHO.py")
print("% Columns: G_k, N_T, p_c [GPa] (soft-mode), p_FES [GPa], DG* [1e-2 kBT/Si], DG* [kBT/cell]")
for m in ['G_1', 'G_2', 'G_3', 'G_4', 'G_5']:
    r = by_member.get(m)
    if r is None:
        print(f"${m}$ & {N_SI[m]} & ${P_C[m]:.4f}$ & --- & --- & --- \\\\  % no usable FES")
    else:
        per_cell = r['barr_si'] * r['n_t']
        print(f"${m}$ & {r['n_t']} & ${r['pc']:.4f}$ & ${r['p_gpa']:.3f}$ "
              f"& ${r['barr_intens']:.1f}$ & ${per_cell:.0f}$ \\\\")
