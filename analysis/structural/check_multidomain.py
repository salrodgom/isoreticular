#!/usr/bin/env python3
"""
Check for spatial heterogeneity / multi-domain artefacts in the OPES Expanded
trajectory.

For each MD snapshot, the PLUMED COLVARS file contains six 8MR distortions
d1..d6 (one per 8MR; pairs (d1,d2), (d3,d4), (d5,d6) belong to the three
orthogonal D8Rs of the cubic cell). In a single-domain configuration the
three D8R-averaged magnitudes should be statistically equivalent at every
snapshot above p_c (the underlying SLC Hamiltonian is cubic-symmetric).

Diagnostics computed below:
  1) D8R averages per snapshot: D_x = (d1+d2)/2, D_y = (d3+d4)/2, D_z = (d5+d6)/2.
  2) Per-snapshot mean μ and CV = std(D_x,D_y,D_z)/μ.
  3) Histogram of CV restricted to the broken-phase window (μ > threshold).

Interpretation:
  - CV < ~0.2 throughout: the three D8Rs are statistically equivalent;
    single-domain regime, your Δ is a clean spatial+time average.
  - CV bimodal or with a tail > 0.5: spatial heterogeneity (one D8R sits in
    Δ ≈ 0 while another is broken). Could be (i) multi-domain in I-43m, (ii)
    a transient during a free-energy basin crossing, or (iii) a numerical
    artefact of the cell.
  - CV broadly distributed but |⟨D_x⟩−⟨D_y⟩|≪std(D_x), i.e. time-averaged
    means are identical: heterogeneity is dynamic only; ensemble Δ unaffected.

Usage:
  cd <Multibaric directory containing ALLCOLVARS.0/1/2>
  python3 ../../scripts_FES/check_multidomain.py
  # or with explicit threshold (in Å):
  python3 check_multidomain.py --threshold 0.5
"""
import sys
import os
import glob
import argparse
import numpy as np

def parse_colvars(fname):
    """Return dict of column-name -> 1D np.array."""
    fields = None
    rows = []
    with open(fname) as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            if line.startswith('#! FIELDS'):
                fields = line.split()[2:]
                continue
            if line.startswith('#'):
                continue
            parts = line.split()
            try:
                rows.append([float(x) for x in parts])
            except ValueError:
                continue
    if fields is None or not rows:
        return None
    arr = np.array(rows)
    if arr.shape[1] != len(fields):
        # tolerate trailing garbage
        ncols = min(arr.shape[1], len(fields))
        arr = arr[:, :ncols]
        fields = fields[:ncols]
    return {name: arr[:, i] for i, name in enumerate(fields)}, arr

def main():
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('--threshold', '-t', type=float, default=None,
                        help='Δ-Å threshold above which "broken phase" diagnostics are computed (default: 0.5 * max(Δ))')
    parser.add_argument('--colvars', '-c', default='ALLCOLVARS',
                        help='glob pattern for the full COLVARS files containing d1..d6 (default: ALLCOLVARS.*; the standard COLVARS.* files only store the aggregated delta and lack the per-8MR columns needed here)')
    parser.add_argument('--plot', action='store_true',
                        help='produce PNG plots if matplotlib is available')
    parser.add_argument('--system', default=None, help='system label for plot titles')
    args = parser.parse_args()

    files = sorted(glob.glob(args.colvars + '.*'))
    files = [f for f in files if not f.endswith('.bak')]
    if not files:
        sys.exit(f"No files matching {args.colvars}.* in {os.getcwd()}; expected ALLCOLVARS.0/1/2")
    print(f"Reading {len(files)} walker files: {files}")

    all_d_x, all_d_y, all_d_z, all_delta = [], [], [], []
    for fn in files:
        parsed = parse_colvars(fn)
        if parsed is None:
            print(f"  warning: {fn} not parseable")
            continue
        cols, _ = parsed
        # Required fields. PLUMED stores δ_8 of each 8MR in d1..d6 in nm;
        # the manuscript Δ is the mean × 10 (Å). We compute the three D8R
        # averages, leaving the factor of 10 implicit (sign- and ratio-
        # diagnostics are scale-invariant).
        needed = ['d1', 'd2', 'd3', 'd4', 'd5', 'd6']
        if not all(k in cols for k in needed):
            print(f"  warning: {fn} missing d1..d6")
            continue
        d1, d2, d3, d4, d5, d6 = [cols[k] for k in needed]
        Dx = 0.5 * (d1 + d2)
        Dy = 0.5 * (d3 + d4)
        Dz = 0.5 * (d5 + d6)
        delta_global = (d1 + d2 + d3 + d4 + d5 + d6) / 6.0
        all_d_x.append(Dx); all_d_y.append(Dy); all_d_z.append(Dz)
        all_delta.append(delta_global)

    Dx = np.concatenate(all_d_x) * 10  # nm -> Å
    Dy = np.concatenate(all_d_y) * 10
    Dz = np.concatenate(all_d_z) * 10
    Dg = np.concatenate(all_delta) * 10
    N = len(Dx)
    print(f"Total snapshots: {N}")
    print()

    # Per-snapshot mean and CV of the three D8R magnitudes
    Dmat = np.column_stack([Dx, Dy, Dz])
    mu  = Dmat.mean(axis=1)
    sig = Dmat.std(axis=1, ddof=0)
    # Coefficient of variation; guard against μ -> 0 in the symmetric phase
    cv  = np.where(mu > 1e-3, sig / mu, np.nan)

    # Threshold to identify "broken-phase" snapshots
    thr = args.threshold if args.threshold is not None else 0.5 * np.nanmax(Dg)
    broken = mu > thr
    n_broken = broken.sum()
    print(f"Threshold for 'broken-phase' diagnostics: μ > {thr:.3f} Å")
    print(f"Snapshots in broken phase: {n_broken} ({100*n_broken/N:.1f}%)")
    print()

    # Diagnostics
    cv_broken = cv[broken & np.isfinite(cv)]
    print("=== Single- vs multi-domain diagnostics (broken phase only) ===")
    print(f"  ⟨D_x⟩ = {Dx[broken].mean():.3f} Å   std = {Dx[broken].std():.3f}")
    print(f"  ⟨D_y⟩ = {Dy[broken].mean():.3f} Å   std = {Dy[broken].std():.3f}")
    print(f"  ⟨D_z⟩ = {Dz[broken].mean():.3f} Å   std = {Dz[broken].std():.3f}")
    print(f"  max |⟨D_i⟩ - ⟨D_j⟩| / σ_D = "
          f"{max(abs(Dx[broken].mean()-Dy[broken].mean()), abs(Dy[broken].mean()-Dz[broken].mean()), abs(Dx[broken].mean()-Dz[broken].mean())) / Dx[broken].std():.3f}")
    print()
    print(f"  CV = σ(D_x,D_y,D_z)/μ per snapshot:")
    print(f"    mean CV  = {cv_broken.mean():.3f}")
    print(f"    median CV = {np.median(cv_broken):.3f}")
    print(f"    95% perc. = {np.percentile(cv_broken, 95):.3f}")
    print(f"    fraction CV > 0.5 (heterogeneous): {100*(cv_broken > 0.5).mean():.1f}%")
    print()
    if cv_broken.mean() < 0.2:
        print("  Conclusion: SINGLE-DOMAIN regime (CV consistently small).")
        print("  The three D8Rs sit at the same magnitude in every broken snapshot;")
        print("  the spatial average over 3 in Δ is statistically clean.")
    elif (cv_broken > 0.5).mean() > 0.1:
        print("  Conclusion: SUSPECTED MULTI-DOMAIN.")
        print("  More than 10% of broken-phase snapshots have one D8R substantially")
        print("  different from the others. The reported Δ is likely diluted by")
        print("  heterogeneous configurations; consider computing the unsigned")
        print("  per-D8R magnitude average separately or restricting to single-")
        print("  domain windows.")
    else:
        print("  Conclusion: MARGINAL.")
        print("  Some heterogeneity exists but no clear bimodality. The ensemble")
        print("  ⟨D_x⟩=⟨D_y⟩=⟨D_z⟩ should still hold; check the time-averaged")
        print("  means above (they should agree within 1 σ).")

    if args.plot:
        try:
            import matplotlib.pyplot as plt
            fig, ax = plt.subplots(1, 2, figsize=(10, 4))
            ax[0].hist(cv_broken, bins=50, color='steelblue')
            ax[0].axvline(0.5, color='red', ls='--', label='heterogeneity threshold')
            ax[0].set_xlabel('CV = σ(D_x,D_y,D_z) / μ (broken phase)')
            ax[0].set_ylabel('count')
            ax[0].set_title('Spatial heterogeneity of the three D8Rs')
            ax[0].legend()
            ax[1].scatter(mu, sig, s=2, alpha=0.4)
            ax[1].set_xlabel('μ = ⟨D_x,D_y,D_z⟩ [Å]')
            ax[1].set_ylabel('σ = std(D_x,D_y,D_z) [Å]')
            ax[1].set_title('Per-snapshot D8R spread vs magnitude')
            ax[1].plot([0, max(mu)], [0, 0.5*max(mu)], 'r--', lw=0.5, label='CV=0.5')
            ax[1].legend()
            sysname = args.system or os.path.basename(os.getcwd())
            fig.suptitle(f'Multi-domain diagnostic — {sysname}')
            fig.tight_layout()
            out = f"multidomain_diag_{sysname}.png"
            fig.savefig(out, dpi=140)
            print(f"\n  Plot saved to {out}")
        except ImportError:
            print("\n  (matplotlib not available; rerun without --plot for text-only output)")

if __name__ == '__main__':
    main()
