#!/usr/bin/env python3
"""
Per-pressure histograms of individual 8MR distortions for each G_k.

For each G_k with available ALLCOLVARS.{0,1,2} files (the full PLUMED2 output
of the OPES Expanded multithermal-multibaric run):

  1. Loads d1..d6 (six individual 8MR distortions in nm), volume and bias.
  2. Bins all snapshots by cell volume (proxy for pressure: lower V <-> higher P).
  3. For each volume bin aggregates the 6 d_i values across snapshots.
  4. Builds a 2D histogram (volume bin x d_i) and a per-bin 1D histogram.

The purpose is to test the hypothesis that G_3 shows a bimodal distribution
of per-8MR distortions in the broken phase (one population in D8Rs, another
in isolated single-8-rings), distinct from the unimodal distribution of the
other G_k.

Output:
  - hist_8MR_per_pressure.png (2D maps for G_1..G_4, side by side)
  - hist_summary.txt (numerical CV, skewness, bimodality coefficient per bin)
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
    fields = None
    rows = []
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

def summary(vals):
    if len(vals) < 5: return dict(n=len(vals), mean=np.nan, std=np.nan, skew=np.nan, bc=np.nan)
    n = len(vals); m = vals.mean(); s = vals.std()
    if s < 1e-12: return dict(n=n, mean=m, std=s, skew=0.0, bc=0.0)
    z = (vals - m) / s
    skew = (z**3).mean()
    kurt = (z**4).mean() - 3
    # Sarle's bimodality coefficient: > 5/9 ~ 0.555 suggests bimodal
    bc = (skew**2 + 1) / (kurt + 3*((n-1)**2)/((n-2)*(n-3)+1e-9))
    return dict(n=n, mean=m, std=s, skew=skew, bc=bc)

# Storage: per system, per volume bin, list of d_i values
all_data = {}

for name, dirpath, nsi, pc, vref in SYSTEMS:
    full = os.path.join(ROOT, dirpath)
    files = sorted(glob.glob(os.path.join(full, 'ALLCOLVARS.*')))
    if not files:
        print(f"{name}: no ALLCOLVARS files in {dirpath}")
        continue
    d_list, vol_list = [], []
    for fn in files:
        flds, arr = parse_colvars(fn)
        if flds is None: continue
        try:
            ds = [flds.index(f'd{i}') for i in range(1, 7)]
            iv = flds.index('vol')
        except ValueError:
            print(f"{name}/{os.path.basename(fn)}: missing d1..d6 or vol")
            continue
        d_list.append(arr[:, ds] * 10)  # nm -> Å
        vol_list.append(arr[:, iv])
    d_all = np.concatenate(d_list, axis=0)   # shape (Nsnap, 6)
    vol = np.concatenate(vol_list)
    print(f"{name}: {d_all.shape[0]} snapshots, vol range [{vol.min():.2f}, {vol.max():.2f}] nm^3")
    all_data[name] = dict(d=d_all, vol=vol, nsi=nsi, pc=pc, vref=vref)

# ============================================================
# 2D histograms
# ============================================================
try:
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    have_mpl = True
except ImportError:
    have_mpl = False
    print("matplotlib not available; skipping plots")

if have_mpl and all_data:
    n_sys = len(all_data)
    fig, axes = plt.subplots(1, n_sys, figsize=(4*n_sys+1, 4.5), sharey=False)
    if n_sys == 1: axes = [axes]
    for ax, (name, d) in zip(axes, all_data.items()):
        d_vals = d['d'].ravel()             # all 6 values flattened
        vol_rep = np.repeat(d['vol'], 6)    # repeat vol for each of 6 d_i
        # Normalize volume to fraction of v0 (vref)
        v_over_v0 = vol_rep / d['vref']
        H, xe, ye = np.histogram2d(v_over_v0, d_vals, bins=[50, 60],
                                    range=[[0.85, 1.05], [0, 2.5]])
        # column-normalise (so each volume bin shows P(d|V))
        H = H / np.maximum(H.sum(axis=1, keepdims=True), 1)
        ax.imshow(H.T, origin='lower', aspect='auto',
                  extent=[xe[0], xe[-1], ye[0], ye[-1]],
                  cmap='magma_r', vmin=0, vmax=H.max()*0.7)
        ax.set_xlabel('$V/V_0$  (proxy for $1/p$)')
        ax.set_ylabel('per-8MR distortion $\\delta_8$ [Å]')
        ax.set_title(f"{name}  ($p_c={d['pc']:.3f}$ GPa, $N={d['d'].shape[0]}$)")
    fig.suptitle('Distribution of individual 8MR distortions vs cell volume (multibaric OPES)')
    fig.tight_layout()
    out = '/sessions/wizardly-eager-fermat/mnt/outputs/hist_8MR_per_pressure.png'
    fig.savefig(out, dpi=140)
    print(f"\n2D map written to {out}")
    plt.close(fig)

# ============================================================
# Per-bin summary table
# ============================================================
print("\n=== Per-volume-bin bimodality of the 6 d_i values per snapshot ===")
print(f"{'system':6s} {'V/V0 bin':14s} {'<d> [Å]':>10s} {'std':>8s} {'skew':>7s} {'BC':>6s}  bimodal? (BC>0.555)")
for name, d in all_data.items():
    v_over_v0 = d['vol'] / d['vref']
    # 6 volume bins
    edges = np.linspace(0.85, 1.05, 7)
    for i in range(6):
        sel = (v_over_v0 >= edges[i]) & (v_over_v0 < edges[i+1])
        if sel.sum() < 30: continue
        vals = d['d'][sel].ravel() * 1.0  # already in Å
        s = summary(vals)
        bim = '*' if s['bc'] > 0.555 else ''
        print(f"{name:6s} [{edges[i]:.3f},{edges[i+1]:.3f}] {s['mean']:>10.3f} {s['std']:>8.3f} {s['skew']:>7.3f} {s['bc']:>6.3f}  {bim}")
