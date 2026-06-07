#!/usr/bin/env python3
"""
Refined extraction of free-energy barriers ΔG* from PLUMED reweighted FES.

Auto-detects fes_<P_bar>.dat files (e.g. fes_12250.dat = 1.225 GPa) in each
Multibaric directory and reports per-pressure barriers in kJ/mol, kBT/Si and
in the manuscript units (100 kBT/Si).

Smoothing: 2D Gaussian (σ=2 px) suppresses kernel-density noise.
Barrier definition: max F along the lowest-cost path (min over cell direction
at each Δ) between the cubic basin (Δ < 0.3 Å min) and the broken basin
(Δ > 0.3 Å min); ΔG* = F_saddle − min(F_cub, F_brk).

Output:
  - stdout: a table with one row per (G_k, pressure)
  - LaTeX-ready summary block at the end, ready to paste into the SI table
"""
import os, re, sys, glob
import numpy as np

KBT_300K = 300 * 0.0083144626  # kJ/mol at 300 K

N_SI = {'G_1': 384, 'G_2': 270, 'G_3': 672, 'G_4': 1440, 'G_5': 2640}
P_C  = {'G_1': 0.940, 'G_2': 0.525, 'G_3': 0.353, 'G_4': 0.199, 'G_5': 0.160}

DIRS = {
    'G_1': 'dir_RHO_isoreticular_G1_222_SG_P1/MultiBaric',
    'G_2': 'dir_RHO_isoreticular_G2_SG_P1/Multibaric',
    'G_3': 'dir_RHO_isoreticular_G3_SG_P1/Multibaric',
    'G_4': 'dir_RHO_isoreticular_G4_SG_P1/Multibaric',
    'G_5': 'dir_RHO_isoreticular_G5_SG_P1/Multibaric',
}

# Resolve root directory relative to this script
ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))

def load_fes_2d(fname):
    cells, deltas, fs = [], [], []
    with open(fname) as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith('#'):
                continue
            parts = line.split()
            try:
                c, d, e = float(parts[0]), float(parts[1]), float(parts[2])
                cells.append(c); deltas.append(d); fs.append(e)
            except (ValueError, IndexError):
                continue
    cells = np.array(cells); deltas = np.array(deltas); fs = np.array(fs)
    udel = np.unique(deltas); ucell = np.unique(cells)
    F = np.full((len(udel), len(ucell)), np.nan)
    for i, d in enumerate(udel):
        for j, c in enumerate(ucell):
            mask = (deltas == d) & (cells == c)
            if mask.any():
                F[i, j] = fs[mask][0]
    return udel, ucell, F

def gauss_smooth(F, sigma=2):
    Fp = np.where(np.isfinite(F), F, np.nanmax(F[np.isfinite(F)]))
    k = int(np.ceil(3*sigma))
    x = np.arange(-k, k+1)
    g = np.exp(-x**2/(2*sigma**2)); g /= g.sum()
    F1 = np.zeros_like(Fp)
    for i in range(Fp.shape[0]):
        row = np.pad(Fp[i], k, mode='edge'); F1[i] = np.convolve(row, g, mode='valid')
    F2 = np.zeros_like(F1)
    for j in range(F1.shape[1]):
        col = np.pad(F1[:, j], k, mode='edge'); F2[:, j] = np.convolve(col, g, mode='valid')
    return F2

def find_basins_and_saddle(udel, ucell, F, delta_split=0.3):
    F = F - np.nanmin(F)
    mask_cub = udel < delta_split
    mask_brk = udel > delta_split
    if mask_cub.sum() < 2 or mask_brk.sum() < 2:
        return None
    cub_local = np.unravel_index(np.nanargmin(F[mask_cub]), F[mask_cub].shape)
    di_cub = np.where(mask_cub)[0][cub_local[0]]; dj_cub = cub_local[1]
    cubic_F = F[di_cub, dj_cub]
    brk_local = np.unravel_index(np.nanargmin(F[mask_brk]), F[mask_brk].shape)
    di_brk = np.where(mask_brk)[0][brk_local[0]]; dj_brk = brk_local[1]
    broken_F = F[di_brk, dj_brk]
    if di_brk == len(udel) - 1:
        return None  # broken basin at grid boundary, not a real minimum
    lo, hi = min(di_cub, di_brk), max(di_cub, di_brk)
    proj = np.nanmin(F[lo:hi+1], axis=1)
    saddle_local = int(np.nanargmax(proj))
    saddle_F = proj[saddle_local]
    saddle_delta = udel[lo + saddle_local]
    return {
        'cubic_delta': udel[di_cub], 'cubic_F': cubic_F,
        'broken_delta': udel[di_brk], 'broken_F': broken_F,
        'saddle_delta': saddle_delta, 'saddle_F': saddle_F,
        'barrier_fwd': saddle_F - cubic_F,
        'barrier_bwd': saddle_F - broken_F,
        'barrier_low': saddle_F - min(cubic_F, broken_F),
    }

def extract_pressure_from_name(fname):
    """fes_12250.dat -> 1.225 GPa, fes_11400.dat -> 1.14 GPa, etc."""
    m = re.search(r'fes_(\d+)\.dat$', os.path.basename(fname))
    if m:
        return int(m.group(1)) / 10000.0  # bar -> GPa
    return None

# === Main ============================================================
results = []  # list of dicts
header = f"{'system':6s} {'file':22s} {'P[GPa]':>7s} {'p/p_c':>7s} {'N_Si':6s} {'2bas':5s} {'Δ_cub':>6s} {'Δ_sad':>6s} {'Δ_brk':>6s} {'ΔG*[kJ/mol]':>12s} {'ΔG*[kBT/Si]':>12s} {'ΔG*[100kBT/Si]':>16s}"
print(header)
print('-' * len(header))

for sys_name, dirpath in DIRS.items():
    full_dir = os.path.join(ROOT, dirpath)
    if not os.path.isdir(full_dir):
        print(f"{sys_name:6s} (Multibaric directory not found: {dirpath})")
        continue
    fes_files = sorted(glob.glob(os.path.join(full_dir, 'fes_*.dat')))
    if not fes_files:
        print(f"{sys_name:6s} (no fes_*.dat files yet; run run_reweight.sh {sys_name[-1]})")
        continue
    n_si = N_SI[sys_name]; pc = P_C[sys_name]
    for fpath in fes_files:
        fname = os.path.basename(fpath)
        p_gpa = extract_pressure_from_name(fpath)
        try:
            udel, ucell, F = load_fes_2d(fpath)
            if F.size == 0: continue
            Fs = gauss_smooth(F, sigma=2)
            r = find_basins_and_saddle(udel, ucell, Fs)
        except Exception as e:
            print(f"{sys_name:6s} {fname:22s} ERROR: {e}")
            continue
        if r is None:
            p_str = f"{p_gpa:.3f}" if p_gpa else "?"
            pr = f"{p_gpa/pc:.2f}" if p_gpa else "?"
            print(f"{sys_name:6s} {fname:22s} {p_str:>7s} {pr:>7s} {n_si:6d} {'no':5s}")
            continue
        barr = r['barrier_low']
        barr_si = barr / (KBT_300K * n_si)
        barr_intens = barr_si / 100
        p_str = f"{p_gpa:.3f}" if p_gpa else "?"
        pr = f"{p_gpa/pc:.2f}" if p_gpa else "?"
        print(f"{sys_name:6s} {fname:22s} {p_str:>7s} {pr:>7s} {n_si:6d} {'yes':5s} "
              f"{r['cubic_delta']:>6.2f} {r['saddle_delta']:>6.2f} {r['broken_delta']:>6.2f} "
              f"{barr:>12.2f} {barr_si:>12.4f} {barr_intens:>16.5f}")
        results.append({
            'system': sys_name, 'p_gpa': p_gpa, 'pc': pc, 'n_si': n_si,
            'barr_kJmol': barr, 'barr_kBT_per_Si': barr_si, 'barr_100kBT_per_Si': barr_intens,
        })

# === LaTeX summary: one row per G_k at the pressure with max barrier =====
print("\n\n=== LaTeX rows for tab:SI-fes-barriers (pressure with max barrier per G_k) ===")
print("% Auto-generated by scripts_FES/extract_barriers.py")
by_sys = {}
for r in results:
    if r['p_gpa'] is None: continue
    s = r['system']
    if s not in by_sys or r['barr_100kBT_per_Si'] > by_sys[s]['barr_100kBT_per_Si']:
        by_sys[s] = r
for s in ['G_1', 'G_2', 'G_3', 'G_4', 'G_5']:
    r = by_sys.get(s)
    if r is None:
        print(f"${s}$ & {N_SI[s]} & ${P_C[s]:.3f}$ & --- & --- & --- \\\\  % no data")
    else:
        per_cell = r['barr_kBT_per_Si'] * r['n_si']
        print(f"${s}$ & {r['n_si']} & ${r['pc']:.3f}$ & ${r['p_gpa']:.3f}$ "
              f"& ${r['barr_100kBT_per_Si']:.4f}$ & ${per_cell:.0f}$ \\\\")
