#!/usr/bin/env python3
"""
δ_8 per 8MR analysis on the STATIC (optimised) structures, for ALL 8MRs in
the framework, across the pressure scan.

Pipeline:
  1. For each G_k, read a reference LAMMPS data file at high-symmetry to
     establish Si-O bonds and Si-Si connectivity (the topology does not
     change with pressure; it does not need to be recomputed per pressure).
  2. Enumerate ALL minimal 8-rings in the Si-Si graph (cycles of length 8
     with no chord).
  3. For each cycle, identify the 8 bridging O atoms in cyclic order around
     the ring.
  4. For each pressure, read the last frame of `global_minimum.lammpstrj`
     (the optimised geometry at that pressure).
  5. Compute δ_8 per 8MR using the Parise-Prince definition with the eight
     bridging oxygens.
  6. Build histograms of δ_8 across all 8MRs, one per (G_k, P).

Output:
  - hist_all_8MR_static.png:   2D histograms (P × δ_8) per G_k
  - all_8MR_summary.txt:       per-pressure stats (mean, std, kurtosis, ...)
"""
import os, sys, glob, re
import numpy as np
from collections import defaultdict
from itertools import combinations

ROOT = '/sessions/wizardly-eager-fermat/mnt/initial_structures_RHO_isoreticular'

SYSTEMS_ALL = [
    ('G_1', 'dir_RHO_isoreticular_G1_SG_P1', 48,  0.940),
    ('G_2', 'dir_RHO_isoreticular_G2_SG_P1', 270, 0.525),
    ('G_3', 'dir_RHO_isoreticular_G3_SG_P1', 672, 0.353),
    ('G_4', 'dir_RHO_isoreticular_G4_SG_P1', 1440, 0.199),
    ('G_5', 'dir_RHO_isoreticular_G5_SG_P1', 2640, 0.160),
]
# Limit via env var ONLY_K to validate first
ONLY_K = os.environ.get('ONLY_K')
if ONLY_K:
    SYSTEMS = [s for s in SYSTEMS_ALL if s[0].endswith(ONLY_K)]
else:
    SYSTEMS = SYSTEMS_ALL

# ---------------------------------------------------------------------
# 1. LAMMPS .data parser  (only needs Atoms section: id type x y z)
# ---------------------------------------------------------------------
def parse_lammps_data(fname):
    atoms = []  # list of (id, type, x, y, z)
    box = None
    with open(fname) as f:
        section = None
        for line in f:
            ls = line.strip()
            if not ls or ls.startswith('#'):
                continue
            if 'xlo xhi' in line:
                parts = line.split()
                xlo, xhi = float(parts[0]), float(parts[1])
                box = box or [None, None, None]
                box[0] = (xlo, xhi)
                continue
            if 'ylo yhi' in line:
                parts = line.split()
                ylo, yhi = float(parts[0]), float(parts[1])
                box[1] = (ylo, yhi); continue
            if 'zlo zhi' in line:
                parts = line.split()
                zlo, zhi = float(parts[0]), float(parts[1])
                box[2] = (zlo, zhi); continue
            if ls.startswith('Atoms'):
                section = 'Atoms'; continue
            if ls.startswith('Bonds') or ls.startswith('Velocities') or ls.startswith('Angles'):
                section = ls.split()[0]; continue
            if ls.startswith('Masses') or ls.startswith('Bond Coeffs') or ls.startswith('Angle Coeffs'):
                section = 'Coeffs'; continue
            if section == 'Atoms':
                parts = ls.split()
                if len(parts) >= 7:
                    try:
                        # LAMMPS atom_style full: id molID type charge x y z
                        aid = int(parts[0])
                        atype = int(parts[2])
                        x, y, z = float(parts[4]), float(parts[5]), float(parts[6])
                        atoms.append((aid, atype, x, y, z))
                    except (ValueError, IndexError):
                        continue
    if box is None:
        box = [(0, 0), (0, 0), (0, 0)]
    return atoms, np.array([box[0][1]-box[0][0], box[1][1]-box[1][0], box[2][1]-box[2][0]])

# ---------------------------------------------------------------------
# 2. LAMMPS trajectory parser (last frame)
# ---------------------------------------------------------------------
def read_lammpstrj_last_frame(fname):
    """Returns dict: positions, box, atom_types (per the trajectory)."""
    with open(fname) as f:
        lines = f.readlines()
    # Find all 'ITEM: TIMESTEP' starts
    starts = [i for i, ln in enumerate(lines) if ln.startswith('ITEM: TIMESTEP')]
    if not starts:
        return None
    start = starts[-1]  # last frame
    # Find sub-items
    n_atoms = int(lines[start+3].strip())
    box = []
    box_start = start + 5
    for k in range(3):
        parts = lines[box_start + k].split()
        box.append((float(parts[0]), float(parts[1])))
    box_lengths = np.array([b[1]-b[0] for b in box])
    # Atoms section
    atoms_hdr = lines[start + 8]
    # Format: ITEM: ATOMS element xs ys zs
    cols = atoms_hdr.split()[2:]  # skip "ITEM:" "ATOMS"
    pos_xs = cols.index('xs') if 'xs' in cols else cols.index('x')
    pos_ys = cols.index('ys') if 'ys' in cols else cols.index('y')
    pos_zs = cols.index('zs') if 'zs' in cols else cols.index('z')
    elem_idx = cols.index('element') if 'element' in cols else 0
    is_scaled = 'xs' in cols
    positions = []
    types = []
    for k in range(n_atoms):
        parts = lines[start + 9 + k].split()
        elem = parts[elem_idx]
        x, y, z = float(parts[pos_xs]), float(parts[pos_ys]), float(parts[pos_zs])
        if is_scaled:
            x *= box_lengths[0]; y *= box_lengths[1]; z *= box_lengths[2]
        positions.append([x, y, z])
        types.append(elem)
    return dict(positions=np.array(positions), types=types, box=box_lengths)

# ---------------------------------------------------------------------
# 3. Build Si-O bonds and Si-Si adjacency
# ---------------------------------------------------------------------
def pbc_distance(r1, r2, box):
    """Minimum-image distance, orthorhombic box."""
    d = r1 - r2
    d -= box * np.round(d / box)
    return np.linalg.norm(d)

def build_si_o_bonds(pos_si, pos_o, box, cutoff=1.95):
    """Returns dict si_idx -> list of o_indices bonded."""
    bonds = defaultdict(list)
    for i, rsi in enumerate(pos_si):
        for j, ro in enumerate(pos_o):
            if pbc_distance(rsi, ro, box) < cutoff:
                bonds[i].append(j)
    return bonds

def build_si_si_graph(si_o_bonds, n_si):
    """Two Si are adjacent if they share an O."""
    # First invert: O -> list of Si
    o_to_si = defaultdict(list)
    for si_i, o_list in si_o_bonds.items():
        for o_j in o_list:
            o_to_si[o_j].append(si_i)
    adj = defaultdict(set)
    bridge_o = {}  # (si_i, si_j) -> o_idx (the bridging O)
    for o_j, si_list in o_to_si.items():
        if len(si_list) == 2:
            a, b = si_list
            adj[a].add(b); adj[b].add(a)
            bridge_o[(min(a, b), max(a, b))] = o_j
    return adj, bridge_o

# ---------------------------------------------------------------------
# 4. 8-ring enumeration in Si-Si graph
# ---------------------------------------------------------------------
def find_all_8rings(adj):
    """Returns list of cycles (each = tuple of 8 Si indices in cyclic order)."""
    seen = set()
    cycles = []
    for u in adj:
        # BFS-style: all paths of length 4 from u
        paths4 = []
        # DFS with depth 4 (5 vertices including endpoints)
        stack = [(u, (u,))]
        while stack:
            node, path = stack.pop()
            if len(path) == 5:
                paths4.append(path)
                continue
            for nbr in adj[node]:
                if nbr not in path:
                    stack.append((nbr, path + (nbr,)))
        # Group by endpoint
        by_v = defaultdict(list)
        for p in paths4:
            by_v[p[-1]].append(p)
        # For each pair of paths sharing only endpoints, form 8-cycle
        for v, plist in by_v.items():
            if v == u or v < u:  # restrict to v > u to halve work
                continue
            for i in range(len(plist)):
                for j in range(i+1, len(plist)):
                    p1 = plist[i]; p2 = plist[j]
                    if set(p1[1:-1]).isdisjoint(set(p2[1:-1])):
                        # 8-cycle: u-p1-v-reverse(p2)-u
                        cyc = p1 + p2[::-1][1:-1]
                        if len(set(cyc)) == 8:
                            key = frozenset(cyc)
                            if key not in seen:
                                seen.add(key)
                                cycles.append(cyc)
    return cycles

def is_minimal_ring(cycle, adj):
    """A ring is minimal if no two non-adjacent vertices in the cycle have a chord."""
    n = len(cycle)
    for i in range(n):
        for j in range(i+2, n):
            if i == 0 and j == n-1:  # adjacent in cycle
                continue
            if cycle[j] in adj[cycle[i]]:
                return False
    return True

# ---------------------------------------------------------------------
# 5. δ_8 per ring from the 8 bridging oxygens (Parise-Prince)
# ---------------------------------------------------------------------
def compute_delta8(cycle, bridge_o, pos_o, box):
    """cycle = tuple of 8 Si indices in cyclic order. Returns δ_8 in Å."""
    # Get the 8 bridging O atoms in cyclic order
    o_indices = []
    n = len(cycle)
    for i in range(n):
        a, b = cycle[i], cycle[(i+1) % n]
        key = (min(a, b), max(a, b))
        if key not in bridge_o:
            return None  # incomplete bridge info
        o_indices.append(bridge_o[key])
    # 4 diagonals
    diag = []
    for k in range(4):
        i = k; j = k + 4
        diag.append(pbc_distance(pos_o[o_indices[i]], pos_o[o_indices[j]], box))
    # Group: cross {0,2} and x {1,3}
    d1 = 0.5 * abs(diag[0] - diag[2])
    d2 = 0.5 * abs(diag[1] - diag[3])
    return max(d1, d2)

# ---------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------
def main():
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt

    all_data = {}
    for name, dirpath, n_si, pc in SYSTEMS:
        full = os.path.join(ROOT, dirpath)
        # Reference data file at 0 bar
        ref_dirs = glob.glob(os.path.join(full, 'dir_0.0_bar'))
        if not ref_dirs:
            print(f"{name}: no dir_0.0_bar found"); continue
        ref_data = glob.glob(os.path.join(ref_dirs[0], '*.data'))
        if not ref_data:
            print(f"{name}: no .data file in dir_0.0_bar"); continue
        atoms, box = parse_lammps_data(ref_data[0])
        # Atom types: 1=Si, 2=O core, 3=O shell (per the masses comment in the data file)
        si_atoms = [(a[0], np.array([a[2], a[3], a[4]])) for a in atoms if a[1] == 1]
        o_atoms  = [(a[0], np.array([a[2], a[3], a[4]])) for a in atoms if a[1] == 2]
        pos_si = np.array([s[1] for s in si_atoms])
        pos_o = np.array([o[1] for o in o_atoms])
        print(f"\n{name}: {len(pos_si)} Si, {len(pos_o)} O, box = {box}")

        si_o = build_si_o_bonds(pos_si, pos_o, box, cutoff=1.95)
        adj, bridge_o = build_si_si_graph(si_o, len(pos_si))
        degrees = [len(adj[i]) for i in range(len(pos_si))]
        print(f"  Si-Si graph: avg degree = {np.mean(degrees):.2f}, min = {min(degrees)}, max = {max(degrees)}")
        if min(degrees) < 4:
            print(f"  WARNING: some Si have < 4 neighbours (cutoff problem?). Skipping {name}.")
            continue

        cycles = find_all_8rings(adj)
        # Filter to minimal
        minimal = [c for c in cycles if is_minimal_ring(c, adj)]
        print(f"  Total 8-cycles found: {len(cycles)}; minimal: {len(minimal)}")
        if not minimal:
            print(f"  No minimal 8-rings; skipping"); continue
        all_data[name] = dict(cycles=minimal, bridge_o=bridge_o, box=box,
                              pos_si=pos_si, pos_o=pos_o, atoms=atoms, pc=pc)

        # Now sweep pressure
        scan_pressures = []
        for pdir in sorted(glob.glob(os.path.join(full, 'dir_*_bar')),
                           key=lambda p: float(re.search(r'dir_([\d.]+)_bar', p).group(1))):
            m = re.search(r'dir_([\d.]+)_bar', pdir)
            if not m: continue
            P_bar = float(m.group(1))
            traj_file = os.path.join(pdir, 'global_minimum.lammpstrj')
            if not os.path.exists(traj_file):
                traj_file = os.path.join(pdir, 'opti.lammpstrj')
            if not os.path.exists(traj_file): continue
            frame = read_lammpstrj_last_frame(traj_file)
            if frame is None: continue
            # Identify Si and O from element labels
            elems = frame['types']
            si_pos = [frame['positions'][i] for i, e in enumerate(elems) if e == 'Si2']
            o_pos  = [frame['positions'][i] for i, e in enumerate(elems) if e == 'O2']
            if len(si_pos) != len(pos_si) or len(o_pos) != len(pos_o):
                continue
            si_pos = np.array(si_pos); o_pos = np.array(o_pos); box_p = frame['box']
            # Compute δ_8 for each minimal 8-ring
            deltas = []
            for c in minimal:
                d8 = compute_delta8(c, bridge_o, o_pos, box_p)
                if d8 is not None:
                    deltas.append(d8)
            if deltas:
                scan_pressures.append((P_bar, deltas))
        if scan_pressures:
            all_data[name]['scan'] = scan_pressures
            print(f"  Pressure scan: {len(scan_pressures)} points, "
                  f"{len(scan_pressures[0][1])} 8-rings per point")

    # ---------- Plot histograms ----------
    n_sys = len(all_data)
    if n_sys == 0:
        print("No data to plot"); return
    fig, axes = plt.subplots(1, n_sys, figsize=(4*n_sys + 1, 4.5), sharey=True)
    if n_sys == 1: axes = [axes]
    for ax, (name, d) in zip(axes, all_data.items()):
        if 'scan' not in d: continue
        pressures = [p / 10000 for p, _ in d['scan']]  # bar -> GPa
        all_deltas_2d = []
        all_pressures_rep = []
        for p_gpa, ds in [(p_bar/10000, ds) for p_bar, ds in d['scan']]:
            all_deltas_2d.extend(ds)
            all_pressures_rep.extend([p_gpa] * len(ds))
        H, xe, ye = np.histogram2d(all_pressures_rep, all_deltas_2d,
                                   bins=[len(set(all_pressures_rep)), 50],
                                   range=[[0, max(pressures)+0.1], [0, 2.5]])
        # Column-normalise
        H = H / np.maximum(H.sum(axis=1, keepdims=True), 1)
        ax.imshow(H.T, origin='lower', aspect='auto',
                  extent=[xe[0], xe[-1], ye[0], ye[-1]],
                  cmap='magma_r', vmin=0, vmax=H.max()*0.7)
        ax.axvline(d['pc'], color='cyan', ls='--', lw=1, label=f"$p_c={d['pc']:.3f}$ GPa")
        ax.set_xlabel('Pressure [GPa]')
        ax.set_title(f"{name}: {len(d['cycles'])} 8MRs/UC")
        ax.legend(loc='upper right', fontsize=8)
    axes[0].set_ylabel('per-8MR distortion $\\delta_8$ [Å]')
    fig.suptitle('Distribution of $\\delta_8$ over ALL 8MRs vs pressure (static optimisations)')
    fig.tight_layout()
    out = '/sessions/wizardly-eager-fermat/mnt/outputs/hist_all_8MR_static.png'
    fig.savefig(out, dpi=140)
    print(f"\nSaved: {out}")

if __name__ == '__main__':
    main()
