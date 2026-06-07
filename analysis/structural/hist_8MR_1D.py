#!/usr/bin/env python3
"""
1D probability density of individual 8MR distortions δ_8, at fixed volume
windows, for each G_k. The 2D map can hide whether a distribution is truly
bimodal or just skewed; this script overlays clean 1D histograms.

Output: hist_8MR_1D.png with one panel per G_k, several volume slices each.
"""
import os, sys, glob
import numpy as np

ROOT = '/sessions/wizardly-eager-fermat/mnt/initial_structures_RHO_isoreticular'

SYSTEMS = [
    ('G_1', 'dir_RHO_isoreticular_G1_222_SG_P1/MultiBaric', 384, 0.940, 21.8),
    ('G_2', 'dir_RHO_isoreticular_G2_SG_P1/Multibaric',     270, 0.525, 14.3),
    ('G_3', 'dir_RHO_isoreticular_G3_SG_P1/Multibaric',     672, 0.353, 39.4),
    ('G_4', 'dir_RHO_isoreticular_G4_SG_P1/Multibaric',    1440, 0.199, 86.0),
]

def parse_colvars(fname):
    fields = None; rows = []
    with open(fname) as f:
        for line in f:
            line = line.strip()
            if not line: continue
            if line.startswith('#! FIELDS'):
                fields = line.split()[2:]
                continue
            if line.startswith('#'): continue
            try:
                rows.append([float(x) for x in line.split()])
            except (ValueError, IndexError):
                continue
    if fields is None or not rows: return None, None
    arr = np.array(rows)
    if arr.shape[1] != len(fields):
        n = min(arr.shape[1], len(fields))
        arr = arr[:, :n]; fields = fields[:n]
    return fields, arr

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib import cm

fig, axes = plt.subplots(2, 2, figsize=(11, 8), sharex=True)
axes = axes.flatten()

# Volume windows (as fractions of v_ref) covering the transition region
WINDOWS = [(0.85, 0.93), (0.93, 0.97), (0.97, 1.00), (1.00, 1.03)]
COLORS = ['#1b4f72', '#2874a6', '#e67e22', '#cb4335']
LABELS = ['v < 0.93 (high P, broken)',
          '0.93 < v < 0.97',
          '0.97 < v < 1.00 (near p_c)',
          'v > 1.00 (low P, cubic)']

for ax, (name, dirpath, nsi, pc, vref) in zip(axes, SYSTEMS):
    full = os.path.join(ROOT, dirpath)
    files = sorted(glob.glob(os.path.join(full, 'ALLCOLVARS.*')))
    d_list, vol_list = [], []
    for fn in files:
        flds, arr = parse_colvars(fn)
        if flds is None: continue
        try:
            ds = [flds.index(f'd{i}') for i in range(1, 7)]
            iv = flds.index('vol')
        except ValueError:
            continue
        d_list.append(arr[:, ds] * 10)  # nm -> Å
        vol_list.append(arr[:, iv])
    d_all = np.concatenate(d_list, axis=0)
    vol = np.concatenate(vol_list)
    v_over_v0 = vol / vref

    for (lo, hi), c, lab in zip(WINDOWS, COLORS, LABELS):
        sel = (v_over_v0 >= lo) & (v_over_v0 < hi)
        if sel.sum() < 30: continue
        vals = d_all[sel].ravel()  # 6 × N_sel per-8MR values, in Å
        # 1D histogram with fine bins
        h, edges = np.histogram(vals, bins=np.linspace(0, 2.6, 100), density=True)
        bin_centres = 0.5 * (edges[1:] + edges[:-1])
        ax.plot(bin_centres, h, color=c, lw=1.6,
                label=f"{lab}  (N={sel.sum()})")

    ax.set_title(f"{name}   p_c = {pc:.3f} GPa,  N_Si = {nsi}")
    ax.set_xlabel(r'per-8MR distortion $\delta_8$ [Å]')
    ax.set_ylabel('probability density')
    ax.legend(fontsize=7, loc='upper right')
    ax.grid(alpha=0.3)

fig.suptitle('1D distributions of individual 8MR distortions per cell-volume window')
fig.tight_layout()
out = '/sessions/wizardly-eager-fermat/mnt/outputs/hist_8MR_1D.png'
fig.savefig(out, dpi=140)
print(f"saved {out}")

# Quantitative test: gap between primary and secondary peaks in the near-p_c window
print("\n=== Peak-finding diagnostic (near p_c window, 0.97 < v/v_ref < 1.00) ===")
print(f"{'system':6s} {'#peaks':>7s} {'peak1 [Å]':>10s} {'peak2 [Å]':>10s} {'gap [Å]':>9s} {'min/max ratio':>14s}")
for name, dirpath, nsi, pc, vref in SYSTEMS:
    full = os.path.join(ROOT, dirpath)
    files = sorted(glob.glob(os.path.join(full, 'ALLCOLVARS.*')))
    d_list, vol_list = [], []
    for fn in files:
        flds, arr = parse_colvars(fn)
        if flds is None: continue
        try:
            ds = [flds.index(f'd{i}') for i in range(1, 7)]
            iv = flds.index('vol')
        except ValueError: continue
        d_list.append(arr[:, ds] * 10); vol_list.append(arr[:, iv])
    d_all = np.concatenate(d_list, axis=0)
    vol = np.concatenate(vol_list); v_over_v0 = vol / vref
    sel = (v_over_v0 >= 0.97) & (v_over_v0 < 1.00)
    if sel.sum() < 100:
        print(f"{name:6s}  (too few points in window: {sel.sum()})"); continue
    vals = d_all[sel].ravel()
    h, edges = np.histogram(vals, bins=np.linspace(0, 2.6, 80), density=True)
    bc = 0.5 * (edges[1:] + edges[:-1])
    # Smooth h with a small kernel
    k = 3; ks = np.ones(k)/k
    hs = np.convolve(h, ks, mode='same')
    # Find local maxima
    peaks = []
    for i in range(2, len(hs)-2):
        if hs[i] > hs[i-1] and hs[i] > hs[i+1] and hs[i] > 0.05 * hs.max():
            peaks.append((bc[i], hs[i]))
    peaks.sort(key=lambda x: -x[1])
    top2 = peaks[:2]
    if len(top2) == 2:
        p1, p2 = sorted(top2, key=lambda x: x[0])
        gap = p2[0] - p1[0]
        # Find the minimum between the two peaks
        i1 = int(np.argmin(np.abs(bc - p1[0])))
        i2 = int(np.argmin(np.abs(bc - p2[0])))
        valley = hs[min(i1, i2):max(i1, i2)+1].min()
        ratio = valley / min(p1[1], p2[1])
        print(f"{name:6s}  {len(peaks):>7d} {p1[0]:>10.3f} {p2[0]:>10.3f} {gap:>9.3f} {ratio:>14.3f}")
    else:
        n = len(peaks)
        p1 = peaks[0] if peaks else (None, None)
        p1_str = f"{p1[0]:.3f}" if p1[0] else "-"
        print(f"{name:6s}  {n:>7d} {p1_str:>10s} {'-':>10s} {'-':>9s} {'-':>14s}  (only one peak)")
