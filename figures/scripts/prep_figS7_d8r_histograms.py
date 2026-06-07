#!/usr/bin/env python3
"""
Extract per-D8R distortion histograms from the fixed-pressure NPT
trajectories for G_1..G_5. For each (G_k, pressure) pair the LAMMPS dump
all.lammpstrj is processed frame by frame; on frame 0 the 8MR topology
of the framework is identified once (BFS + SSSR + geometric ring filter)
and the resulting list of Si-cycles + bridging-oxygen map is reused for
the remaining frames.

For each frame and each D8R the Parise-Prince elliptic distortion
    delta_8 = (1/2) max(|d_1 - d_3|, |d_2 - d_4|)
is computed (d_i are the four oxygen-oxygen diagonals of the ring) and
the value (in Angstrom) is appended to the aggregated sample.

Target pressures per G_k follow the user's instruction (0, 0.5 p_c, p_c,
1.5 p_c). For each target, the nearest pressure directory available in
initial_structures_RHO_isoreticular/ is chosen.

Writes:
    data/figS7_d8r_G{k}_p{label}.dat
        single column: per-ring Delta [A]
    data/figS7_d8r_meta.txt
        per-(G_k, target pressure) summary: actual p [GPa], dir name, n_rings,
        n_frames, n_samples, mean, std.

Run from the repo root:
    python3 scripts/prep_figS7_d8r_histograms.py
"""
import os
import sys
import glob
import re
import numpy as np

HERE_PREP = os.path.dirname(os.path.abspath(__file__))
RHO_ROOT  = os.path.normpath(os.path.join(HERE_PREP, '..', '..',
                                          'initial_structures_RHO_isoreticular'))
sys.path.insert(0, os.path.join(RHO_ROOT, 'scripts_FES'))
import ring_distortions as _rd
from ring_distortions import (find_all_minimal_rings,
                              distortion_4MR, distortion_6MR, distortion_8MR)
from collections import defaultdict


def build_graph_fast(pos_si, pos_o, box, cutoff=1.95):
    """Vectorised replacement for ring_distortions.build_graph. The original
    nests two pure-Python loops over (Si, O) pairs with a per-pair numpy
    PBC distance call, which scales O(N_Si x N_O) but with very high
    per-iteration overhead (it timed out on G_5 with 2640 Si x 5280 O).

    This version uses one numpy broadcasted distance computation per Si,
    keeping memory cost to a single (N_O, 3) array but eliminating the
    Python-level inner loop.
    """
    box_arr = np.asarray(box, dtype=float)
    cutoff2 = cutoff ** 2
    o_to_si = defaultdict(list)
    for i, rs in enumerate(pos_si):
        d = pos_o - rs                       # (N_O, 3)
        d -= box_arr * np.round(d / box_arr)
        d2 = np.einsum('ij,ij->i', d, d)
        for k in np.where(d2 < cutoff2)[0]:
            o_to_si[int(k)].append(i)
    adj = defaultdict(set)
    bridge = {}
    for k, sl in o_to_si.items():
        if len(sl) == 2:
            a, b = sl
            adj[a].add(b); adj[b].add(a)
            bridge[(min(a, b), max(a, b))] = k
    return adj, bridge


# Monkey-patch the slow build_graph used by find_all_minimal_rings
# (analyse_structure also calls it, but we bypass analyse_structure here).
_rd.build_graph = build_graph_fast
build_graph = build_graph_fast


SLC_DIR = {
    1: 'dir_RHO_isoreticular_G1_222_SG_P1',
    2: 'dir_RHO_isoreticular_G2_SG_P1',
    3: 'dir_RHO_isoreticular_G3_SG_P1',
    4: 'dir_RHO_isoreticular_G4_SG_P1',
    5: 'dir_RHO_isoreticular_G5_SG_P1',
}
# Critical pressures (GPa) from Table 1 of the manuscript.
PC = {1: 0.940, 2: 0.525, 3: 0.353, 4: 0.199, 5: 0.160}
# Target pressure factors relative to p_c. A factor of None denotes an
# absolute pressure (in GPa) supplied separately via P_ABSOLUTE_TARGETS.
P_FACTORS = (0.0, 0.5, 1.0, 1.5, None)
# Labels for filenames (no '.' or '/').
P_LABELS = ('p0', 'p05pc', 'ppc', 'p15pc', 'p2gpa')
# Absolute pressures (GPa) for the entries with pf = None. Same length as
# P_FACTORS; only the None positions are used.
P_ABSOLUTE_TARGETS = (None, None, None, None, 2.0)

# Common dense pressure list (in bar) used by Figure S8 (D8R-only) and shared
# across all G_k so that the per-panel histograms can be compared directly
# (same pressures => single legend in the first panel).
DENSE_BARS = (0, 1000, 2000, 4000, 6000, 8000, 10000, 15000, 20000)

# Explicit pressure-directory choice per (G_k, pressure factor). For G_k
# with sparse trajectory dumps (notably G_5), the auto "closest-available"
# rule produced collisions; this table is the curated mapping used in the
# manuscript figure. Each entry maps to the pressure in bar of the dir_*_bar
# directory that contains the all.lammpstrj for that point.
EXPLICIT_DIR_BAR = {
    # Last entry (None key, absolute 2 GPa) uses dir_20000.0_bar for every G_k.
    1: {0.0: 0,    0.5: 5000,  1.0: 9400,  1.5: 14000, '2.0': 20000},
    2: {0.0: 0,    0.5: 2000,  1.0: 5250,  1.5: 8000,  '2.0': 20000},
    3: {0.0: 0,    0.5: 2000,  1.0: 3500,  1.5: 5000,  '2.0': 20000},
    4: {0.0: 0,    0.5: 1000,  1.0: 2000,  1.5: 3000,  '2.0': 20000},
    # G_5 has lammpstrj only at sparse 1000-bar increments (plus 350, 750).
    # Best practical mapping (pc=0.160 GPa = 1600 bar):
    5: {0.0: 0,    0.5: 750,   1.0: 1000,  1.5: 3000,  '2.0': 20000},
}


def list_pressure_dirs(g_dir):
    """Return (p_GPa, dirname) for every dir_*_bar entry that has an
    all.lammpstrj trajectory. Dirs without a trajectory (typically the
    minimisation-only fine pressure-scan dirs) are skipped."""
    out = []
    for sub in sorted(glob.glob(os.path.join(g_dir, 'dir_*_bar'))):
        m = re.search(r'dir_([\d.]+)_bar', os.path.basename(sub))
        if not m:
            continue
        if not os.path.exists(os.path.join(sub, 'all.lammpstrj')):
            continue
        p_bar = float(m.group(1))
        out.append((p_bar / 10000.0, sub))
    return sorted(out, key=lambda t: t[0])


def closest_pressure_dir(g_dir, p_target_GPa):
    available = list_pressure_dirs(g_dir)
    if not available:
        return None, None
    closest = min(available, key=lambda t: abs(t[0] - p_target_GPa))
    return closest


def iter_lammpstrj_frames(fname, stride=1):
    """Yield (box, pos_si, pos_o) tuples per frame of a lammpstrj file.
    The trajectory uses scaled (xs, ys, zs) coordinates, atom labels
    starting with 'Si' or 'O' for the framework. Triclinic tilt factors
    in BOX BOUNDS are ignored (assumed << box length).
    """
    with open(fname) as f:
        frame_index = 0
        while True:
            line = f.readline()
            if not line:
                return
            if not line.startswith('ITEM: TIMESTEP'):
                continue
            f.readline()  # timestep value
            assert f.readline().startswith('ITEM: NUMBER OF ATOMS')
            natoms = int(f.readline())
            assert f.readline().startswith('ITEM: BOX BOUNDS')
            xlo, xhi, *_ = (float(x) for x in f.readline().split())
            ylo, yhi, *_ = (float(x) for x in f.readline().split())
            zlo, zhi, *_ = (float(x) for x in f.readline().split())
            assert f.readline().startswith('ITEM: ATOMS')
            ax = xhi - xlo
            ay = yhi - ylo
            az = zhi - zlo
            si, ox = [], []
            for _ in range(natoms):
                parts = f.readline().split()
                if not parts:
                    continue
                sym = parts[0]
                fx, fy, fz = float(parts[1]), float(parts[2]), float(parts[3])
                if sym.lower().startswith('si'):
                    si.append((fx * ax, fy * ay, fz * az))
                elif sym.lower().startswith('o'):
                    ox.append((fx * ax, fy * ay, fz * az))
            if frame_index % stride == 0:
                yield (np.array([ax, ay, az]),
                       np.array(si, dtype=float),
                       np.array(ox, dtype=float))
            frame_index += 1


def classify_d8r_pairs(cycles_8, adj):
    """Identify D8R pairs among the 8MR rings.

    Two 8MR cycles A and B form a D8R if (i) they are vertex-disjoint and
    (ii) every Si in A has exactly one bonded neighbour that lies in B.
    Total: 8 inter-ring Si-Si bonds (= 8 bridging O atoms) connect the two
    rings, the defining geometry of a face-to-face double 8-ring.

    Returns the set of ring indices that participate in at least one D8R
    pair.
    """
    n = len(cycles_8)
    sets_8 = [set(c) for c in cycles_8]
    in_d8r = set()
    for i in range(n):
        A = sets_8[i]
        for j in range(i + 1, n):
            B = sets_8[j]
            if A & B:
                continue
            # Count Si in A that have a bonded neighbour in B.
            matched = sum(1 for a in A if any(b in B for b in adj[a]))
            if matched == 8:
                in_d8r.add(i)
                in_d8r.add(j)
    return in_d8r


def process_one(g_dir, p_target, p_bar_override=None, stride=1, max_frames=None):
    if p_bar_override is not None:
        # Format with .0 because dir names are dir_<N>.0_bar even for integer N
        sub = os.path.join(g_dir, f'dir_{float(p_bar_override):.1f}_bar')
        p_actual = p_bar_override / 10000.0
    else:
        p_actual, sub = closest_pressure_dir(g_dir, p_target)
        if sub is None:
            return None
    traj = os.path.join(sub, 'all.lammpstrj')
    if not os.path.exists(traj):
        return None
    frames = iter_lammpstrj_frames(traj, stride=stride)

    # Frame 0: detect rings.
    try:
        box, pos_si, pos_o = next(frames)
    except StopIteration:
        return None
    adj, bridge = build_graph(pos_si, pos_o, box, cutoff=1.95)
    rings_all = find_all_minimal_rings(adj, pos_si, box, lengths=(4, 6, 8))
    cycles_4 = rings_all[4]
    cycles_6 = rings_all[6]
    cycles_8 = rings_all[8]
    n_rings = len(cycles_8)
    if n_rings == 0:
        return None

    # Identify which rings belong to a D8R pair (face-to-face dimer).
    d8r_set = classify_d8r_pairs(cycles_8, adj)
    n_d8r_rings = len(d8r_set)
    n_4mr = len(cycles_4)
    n_6mr = len(cycles_6)

    deltas_all = []
    deltas_d8r = []
    deltas_iso = []
    lambdas_4 = []
    gammas_6 = []

    def collect(box, pos_o):
        for idx, c in enumerate(cycles_8):
            d = distortion_8MR(c, bridge, pos_o, box)
            if d is not None:
                deltas_all.append(d)
                if idx in d8r_set:
                    deltas_d8r.append(d)
                else:
                    deltas_iso.append(d)
        for c in cycles_4:
            lam = distortion_4MR(c, bridge, pos_o, box)
            if lam is not None:
                lambdas_4.append(lam)
        for c in cycles_6:
            gam = distortion_6MR(c, bridge, pos_o, box)
            if gam is not None:
                gammas_6.append(gam)

    # Frame 0
    collect(box, pos_o)
    n_frames = 1
    # Subsequent frames
    for box, pos_si, pos_o in frames:
        collect(box, pos_o)
        n_frames += 1
        if max_frames is not None and n_frames >= max_frames:
            break

    return dict(
        p_actual=p_actual,
        dir=os.path.basename(sub),
        n_rings=n_rings,
        n_d8r_rings=n_d8r_rings,
        n_iso_rings=n_rings - n_d8r_rings,
        n_4mr=n_4mr,
        n_6mr=n_6mr,
        n_frames=n_frames,
        deltas=np.array(deltas_all),
        deltas_d8r=np.array(deltas_d8r),
        deltas_iso=np.array(deltas_iso),
        lambdas_4=np.array(lambdas_4),
        gammas_6=np.array(gammas_6),
    )


def main():
    out_dir = os.path.normpath(os.path.join(HERE_PREP, '..', 'data'))
    os.makedirs(out_dir, exist_ok=True)
    meta_path = os.path.join(out_dir, 'figS7_d8r_meta.txt')
    # Optional CLI: --k <int>   restricts processing to a single G_k
    #                 --pf <flt>  restricts to a single pressure factor in P_FACTORS
    only_k = None
    only_pf = None
    if '--k' in sys.argv:
        only_k = int(sys.argv[sys.argv.index('--k') + 1])
    if '--pf' in sys.argv:
        only_pf = float(sys.argv[sys.argv.index('--pf') + 1])

    if '--dense' in sys.argv:
        # Dense mode: ignore pf, sweep DENSE_BARS for the requested k (or all).
        only_bar = None
        if '--bar' in sys.argv:
            only_bar = int(sys.argv[sys.argv.index('--bar') + 1])
        for k in SLC_DIR:
            if only_k is not None and k != only_k:
                continue
            g_dir = os.path.join(RHO_ROOT, SLC_DIR[k])
            for p_bar in DENSE_BARS:
                if only_bar is not None and p_bar != only_bar:
                    continue
                fname_d8r = os.path.join(out_dir,
                                         f'figS8d_d8ronly_G{k}_p{p_bar}bar.dat')
                fname_iso = os.path.join(out_dir,
                                         f'figS9d_isolated_G{k}_p{p_bar}bar.dat')
                fname_4mr = os.path.join(out_dir,
                                         f'figS13d_4mr_G{k}_p{p_bar}bar.dat')
                fname_6mr = os.path.join(out_dir,
                                         f'figS13d_6mr_G{k}_p{p_bar}bar.dat')
                # Skip only if all expected outputs already exist.
                all_done = (os.path.exists(fname_d8r)
                            and (k == 1 or os.path.exists(fname_iso))
                            and os.path.exists(fname_4mr)
                            and os.path.exists(fname_6mr))
                if all_done:
                    continue
                p_GPa = p_bar / 10000.0
                # For the largest system (G_5) the 4+6+8 MR loop per frame
                # is heavy enough to bust the workspace bash timeout, so we
                # halve the frame rate. The shape of the histograms is
                # essentially unchanged (~50 frames x 1020 rings still gives
                # 50k+ samples per pressure).
                stride = 2 if k == 5 else 1
                res = process_one(g_dir, p_GPa, p_bar_override=p_bar,
                                  stride=stride)
                if res is None:
                    print(f'  dense G_{k}  p={p_bar:5d} bar  no data')
                    continue
                if not os.path.exists(fname_d8r):
                    np.savetxt(fname_d8r, res['deltas_d8r'],
                               header=(f'G_{k}, p = {p_GPa:.4f} GPa (dir_{p_bar}.0_bar).\n'
                                       f'D8R-only: {res["n_d8r_rings"]} rings x '
                                       f'{res["n_frames"]} frames = '
                                       f'{len(res["deltas_d8r"])} samples.\n'
                                       'Column: per-D8R Delta [A]'),
                               fmt='%.6f')
                if k > 1 and not os.path.exists(fname_iso) and res['n_iso_rings'] > 0:
                    np.savetxt(fname_iso, res['deltas_iso'],
                               header=(f'G_{k}, p = {p_GPa:.4f} GPa (dir_{p_bar}.0_bar).\n'
                                       f'Isolated-8MR (not in any D8R): '
                                       f'{res["n_iso_rings"]} rings x '
                                       f'{res["n_frames"]} frames = '
                                       f'{len(res["deltas_iso"])} samples.\n'
                                       'Column: per-isolated-8MR Delta [A]'),
                               fmt='%.6f')
                if not os.path.exists(fname_4mr) and res['n_4mr'] > 0:
                    np.savetxt(fname_4mr, res['lambdas_4'],
                               header=(f'G_{k}, p = {p_GPa:.4f} GPa (dir_{p_bar}.0_bar).\n'
                                       f'4MR distortion Lambda = (1/2)|d_02 - d_13|.\n'
                                       f'{res["n_4mr"]} rings x {res["n_frames"]} frames'
                                       f' = {len(res["lambdas_4"])} samples.\n'
                                       'Column: per-4MR Lambda [A]'),
                               fmt='%.6f')
                if not os.path.exists(fname_6mr) and res['n_6mr'] > 0:
                    np.savetxt(fname_6mr, res['gammas_6'],
                               header=(f'G_{k}, p = {p_GPa:.4f} GPa (dir_{p_bar}.0_bar).\n'
                                       f'6MR distortion Gamma = (1/2)(max - min) over the three '
                                       'opposite-O diagonals.\n'
                                       f'{res["n_6mr"]} rings x {res["n_frames"]} frames'
                                       f' = {len(res["gammas_6"])} samples.\n'
                                       'Column: per-6MR Gamma [A]'),
                               fmt='%.6f')
                m = float(res['deltas_d8r'].mean())
                s = float(res['deltas_d8r'].std())
                miso = float(res['deltas_iso'].mean()) if len(res['deltas_iso']) else float('nan')
                print(f'  dense G_{k}  p={p_bar:5d} bar = {p_GPa:.2f} GPa  '
                      f'D8R={res["n_d8r_rings"]} iso={res["n_iso_rings"]} '
                      f'<Δ>_d8r={m:.3f}±{s:.3f}  <Δ>_iso={miso:.3f} Å')
        return
    # Append mode if processing a single G_k (so we can build up the meta
    # incrementally across calls); otherwise overwrite.
    meta_mode = 'a' if only_k is not None and os.path.exists(meta_path) else 'w'
    with open(meta_path, meta_mode) as meta:
        if meta_mode == 'w':
            meta.write('# Per-D8R distortion histograms extracted from MD-NPT '
                       'trajectories.\n')
            meta.write('# Columns: k  pf  p_actual[GPa]  dir  n_rings  n_frames  '
                       'n_samples  mean[A]  std[A]\n')
        for k, sub_dir in SLC_DIR.items():
            if only_k is not None and k != only_k:
                continue
            g_dir = os.path.join(RHO_ROOT, sub_dir)
            pc = PC[k]
            for pf, label, pabs in zip(P_FACTORS, P_LABELS, P_ABSOLUTE_TARGETS):
                if only_pf is not None:
                    # accept either matching pf or the absolute-pressure marker
                    if pf is None:
                        if abs(only_pf - 99.0) > 1e-6:
                            continue
                    elif abs(pf - only_pf) > 1e-6:
                        continue
                if pf is None:
                    p_target = pabs
                    key = '2.0'
                else:
                    p_target = pf * pc
                    key = pf
                p_bar_override = EXPLICIT_DIR_BAR.get(k, {}).get(key)
                res = process_one(g_dir, p_target,
                                  p_bar_override=p_bar_override, stride=1)
                if res is None:
                    print(f'  G_{k}  pf={pf:.2f}  no data')
                    continue
                if pf is None:
                    pf_str = f'absolute {pabs:.2f} GPa'
                else:
                    pf_str = f'{pf:.2f} p_c'
                # All 8MRs (Figure S7).
                fname = os.path.join(out_dir,
                                     f'figS7_d8r_G{k}_{label}.dat')
                np.savetxt(fname, res['deltas'],
                           header=(f'G_{k}, target p = {pf_str} = '
                                   f'{p_target:.4f} GPa.  Actual dir = '
                                   f'{res["dir"]} ({res["p_actual"]:.4f} GPa).\n'
                                   f'All 8MRs: {res["n_rings"]} rings x '
                                   f'{res["n_frames"]} frames = '
                                   f'{len(res["deltas"])} samples.\n'
                                   'Column: per-8MR Delta [A]'),
                           fmt='%.6f')
                # D8R-only filter (Figure S8).
                fname2 = os.path.join(out_dir,
                                      f'figS8_d8ronly_G{k}_{label}.dat')
                np.savetxt(fname2, res['deltas_d8r'],
                           header=(f'G_{k}, target p = {pf_str} = '
                                   f'{p_target:.4f} GPa.  Actual dir = '
                                   f'{res["dir"]} ({res["p_actual"]:.4f} GPa).\n'
                                   f'D8R-only: {res["n_d8r_rings"]} rings x '
                                   f'{res["n_frames"]} frames = '
                                   f'{len(res["deltas_d8r"])} samples.\n'
                                   'Column: per-D8R Delta [A]'),
                           fmt='%.6f')
                # Isolated 8MR filter (Figure S9): complement of D8R-only.
                if res['n_iso_rings'] > 0:
                    fname3 = os.path.join(out_dir,
                                          f'figS9_isolated_G{k}_{label}.dat')
                    np.savetxt(fname3, res['deltas_iso'],
                               header=(f'G_{k}, target p = {pf_str} = '
                                       f'{p_target:.4f} GPa.  Actual dir = '
                                       f'{res["dir"]} ({res["p_actual"]:.4f} GPa).\n'
                                       f'Isolated-8MR (not in any D8R): '
                                       f'{res["n_iso_rings"]} rings x '
                                       f'{res["n_frames"]} frames = '
                                       f'{len(res["deltas_iso"])} samples.\n'
                                       'Column: per-isolated-8MR Delta [A]'),
                               fmt='%.6f')
                m = float(res['deltas'].mean())
                s = float(res['deltas'].std())
                m2 = float(res['deltas_d8r'].mean()) if len(res['deltas_d8r']) else float('nan')
                s2 = float(res['deltas_d8r'].std()) if len(res['deltas_d8r']) else float('nan')
                pf_print = f'{pf:.2f}' if pf is not None else 'abs'
                print(f'  G_{k}  pf={pf_print}  p={res["p_actual"]:.3f} GPa  '
                      f'all8MR={res["n_rings"]} d8r-only={res["n_d8r_rings"]} '
                      f'frames={res["n_frames"]}  '
                      f'<Δ>_all={m:.3f}±{s:.3f}  <Δ>_d8r={m2:.3f}±{s2:.3f} Å')
                meta.write(f'{k:d}  {pf_print}  {res["p_actual"]:.4f}  '
                           f'{res["dir"]}  {res["n_rings"]:d}  '
                           f'{res["n_d8r_rings"]:d}  '
                           f'{res["n_frames"]:d}  {len(res["deltas"]):d}  '
                           f'{m:.6f}  {s:.6f}  {m2:.6f}  {s2:.6f}\n')


if __name__ == '__main__':
    main()
