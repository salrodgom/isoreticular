#!/usr/bin/env python3
"""
Geometric analysis of pore-window distortions in RHO-family zeolites.

Builds a Si-Si graph through O bridges from a LAMMPS data file or CIF,
enumerates minimal rings of length 4, 6 and 8 (no chord; cycle distance
equals graph distance for non-adjacent pairs), and computes the
Parise-Prince elliptic distortion (Γ for 4MR, Λ for 6MR, Δ for 8MR) for
each ring using the cyclically-ordered bridging oxygens:

    δ_4 = (1/2) |r_ik - r_jl|                                   (4MR)
    δ_6 = (1/2) (max - min over the three opposite-O pairs)     (6MR)
    δ_8 = (1/2) max(|r_im - r_ko|, |r_jn - r_pl|)               (8MR)

Output:
  - per-ring tabulation (length, atoms, distortion)
  - histograms of δ_n per ring class
  - aggregate Γ/Λ/Δ per structure (mean over the ring class)
"""
import os, sys, re, glob
import numpy as np
from collections import defaultdict
import argparse

# -----------------------------------------------------------------------
# 1. Structure parsing
# -----------------------------------------------------------------------
def parse_lammps_data(fname):
    """Parse a LAMMPS data file (atom_style full). Returns Si and O positions
    in Å, and the orthorhombic cell vector."""
    atoms = []; box = [None]*3
    section = None
    with open(fname) as f:
        for line in f:
            ls = line.strip()
            if 'xlo xhi' in line:
                p = line.split(); box[0] = float(p[1]) - float(p[0]); continue
            if 'ylo yhi' in line:
                p = line.split(); box[1] = float(p[1]) - float(p[0]); continue
            if 'zlo zhi' in line:
                p = line.split(); box[2] = float(p[1]) - float(p[0]); continue
            if ls.startswith('Atoms'):
                section = 'Atoms'; continue
            if ls.startswith(('Bonds', 'Velocities', 'Angles', 'Masses',
                              'Bond Coeffs', 'Angle Coeffs', 'Pair Coeffs')):
                section = ls.split()[0]; continue
            if section == 'Atoms':
                p = ls.split()
                if len(p) >= 7:
                    try:
                        aid = int(p[0])
                        atype = int(p[2])  # atom_style full: id molID type charge x y z
                        x, y, z = float(p[4]), float(p[5]), float(p[6])
                        atoms.append((aid, atype, x, y, z))
                    except (ValueError, IndexError):
                        continue
    return atoms, np.array(box)

def parse_cif_cubic(fname):
    """Parse a CIF for a cubic cell. Returns Si and O fractional positions
    and the cell parameter a."""
    cell = {}; atoms = []
    in_loop = False; cols = []
    with open(fname) as f:
        for line in f:
            ls = line.strip()
            if ls.startswith('_cell_length_a '): cell['a'] = float(ls.split()[1])
            elif ls.startswith('_cell_length_b '): cell['b'] = float(ls.split()[1])
            elif ls.startswith('_cell_length_c '): cell['c'] = float(ls.split()[1])
            elif ls == 'loop_':
                in_loop = True; cols = []
            elif in_loop and ls.startswith('_'):
                cols.append(ls)
            elif in_loop and ls and not ls.startswith('_'):
                if 'atom_site' not in (cols[0] if cols else ''):
                    in_loop = False; cols = []; continue
                p = ls.split()
                try:
                    ix = next(i for i,c in enumerate(cols) if c.endswith('fract_x'))
                    iy = next(i for i,c in enumerate(cols) if c.endswith('fract_y'))
                    iz = next(i for i,c in enumerate(cols) if c.endswith('fract_z'))
                    il = next(i for i,c in enumerate(cols)
                              if c.endswith('type_symbol') or c.endswith('label'))
                except StopIteration:
                    in_loop = False; cols = []; continue
                if len(p) > max(ix, iy, iz, il):
                    sym = p[il].rstrip('0123456789')
                    atoms.append((sym, float(p[ix]), float(p[iy]), float(p[iz])))
    return cell, atoms

# -----------------------------------------------------------------------
# 2. Build Si-O bonds and Si-Si graph
# -----------------------------------------------------------------------
def pbc_dist(r1, r2, box):
    d = r1 - r2
    d -= box * np.round(d / box)
    return float(np.linalg.norm(d))

def build_graph(pos_si, pos_o, box, cutoff=1.95):
    """Returns:
      adj[i]      = set of Si neighbours of Si i
      bridge[(i,j)] = O index bridging Si i and Si j (i < j)
    """
    si_o = defaultdict(list)
    for i, rs in enumerate(pos_si):
        for k, ro in enumerate(pos_o):
            if pbc_dist(rs, ro, box) < cutoff:
                si_o[i].append(k)
    o_to_si = defaultdict(list)
    for i, oo in si_o.items():
        for k in oo:
            o_to_si[k].append(i)
    adj = defaultdict(set); bridge = {}
    for k, sl in o_to_si.items():
        if len(sl) == 2:
            a, b = sl
            adj[a].add(b); adj[b].add(a)
            bridge[(min(a, b), max(a, b))] = k
    return adj, bridge

# -----------------------------------------------------------------------
# 3. Ring enumeration (minimal, no chord)
# -----------------------------------------------------------------------
def find_rings_of_length(adj, L):
    """Returns list of tuples (Si indices in cyclic order), one per minimal ring."""
    seen = set(); rings = []
    for u in sorted(adj):
        stack = [(u, (u,))]
        while stack:
            node, path = stack.pop()
            if len(path) == L:
                # Try to close
                if u in adj[node]:
                    cyc = path
                    # Canonical form (rotation + reflection invariant)
                    fs = frozenset(cyc)
                    if fs not in seen:
                        seen.add(fs)
                        rings.append(cyc)
                continue
            for nbr in adj[node]:
                if nbr not in path:
                    stack.append((nbr, path + (nbr,)))
    return rings

def bfs_distances(adj, n_nodes):
    """Returns NxN matrix of shortest-path distances between all node pairs."""
    INF = 10**9
    D = np.full((n_nodes, n_nodes), INF, dtype=np.int32)
    for start in range(n_nodes):
        D[start, start] = 0
        from collections import deque
        q = deque([start])
        while q:
            u = q.popleft()
            for v in adj.get(u, ()):
                if D[start, v] == INF:
                    D[start, v] = D[start, u] + 1
                    q.append(v)
    return D

def is_sssr_minimal(cycle, dist):
    """Strict SSSR criterion: for every pair (c_i, c_j) of cycle vertices, the
    graph distance must be at least the cycle distance. If a shortcut exists
    (graph distance < cycle distance), the cycle is composite."""
    n = len(cycle)
    for i in range(n):
        for j in range(i+1, n):
            cd = min(j - i, n - (j - i))
            if dist[cycle[i], cycle[j]] < cd:
                return False
    return True

def unwrap_cycle_coords(cycle, pos, box):
    """Unwrap cycle vertex coordinates using minimum image convention."""
    coords = [pos[cycle[0]].copy()]
    for i in range(1, len(cycle)):
        prev = coords[-1]
        curr = pos[cycle[i]].copy()
        diff = curr - prev
        diff -= box * np.round(diff / box)
        coords.append(prev + diff)
    return np.array(coords)

def is_geometric_ring(cycle, pos_si, box, planarity_tol=0.30, radius_tol=0.20):
    """Test that the cycle vertices form a quasi-planar ring with similar
    centroid-to-vertex distances.

    planarity_tol = ratio σ_3 / σ_1 of singular values of the centered
                    coords (lower = more planar)
    radius_tol    = std / mean of centroid-to-vertex distances
    """
    coords = unwrap_cycle_coords(cycle, pos_si, box)
    centroid = coords.mean(axis=0)
    rs = np.linalg.norm(coords - centroid, axis=1)
    if rs.mean() < 1e-6: return False
    if rs.std() / rs.mean() > radius_tol:
        return False
    centered = coords - centroid
    s = np.linalg.svd(centered, compute_uv=False)
    if s[0] < 1e-9: return False
    if s[2] / s[0] > planarity_tol:
        return False
    return True

def find_all_minimal_rings(adj, pos_si, box, lengths=(4, 6, 8)):
    """Returns dict {L: list of cycles}. Each cycle is a tuple of Si indices in cyclic order.
    Filtered by SSSR minimality (no shortcut) AND geometric ring-likeness
    (planar, equidistant from centroid)."""
    n_nodes = max(adj.keys()) + 1 if adj else 0
    dist = bfs_distances(adj, n_nodes)
    out = {}
    for L in lengths:
        rings = find_rings_of_length(adj, L)
        accepted = []
        for c in rings:
            if not is_sssr_minimal(c, dist): continue
            if not is_geometric_ring(c, pos_si, box): continue
            accepted.append(c)
        out[L] = accepted
    return out

# -----------------------------------------------------------------------
# 4. Per-ring distortion
# -----------------------------------------------------------------------
def ring_bridging_oxygens(cycle, bridge):
    n = len(cycle); o = []
    for k in range(n):
        a, b = cycle[k], cycle[(k+1) % n]
        key = (min(a, b), max(a, b))
        if key not in bridge: return None
        o.append(bridge[key])
    return o

def distortion_4MR(cycle, bridge, pos_o, box):
    """δ_4 = (1/2)|d_02 - d_13|"""
    o = ring_bridging_oxygens(cycle, bridge)
    if o is None: return None
    d02 = pbc_dist(pos_o[o[0]], pos_o[o[2]], box)
    d13 = pbc_dist(pos_o[o[1]], pos_o[o[3]], box)
    return 0.5 * abs(d02 - d13)

def distortion_6MR(cycle, bridge, pos_o, box):
    """δ_6 = (1/2)(max - min over the three opposite-O pairs (0,3), (1,4), (2,5))"""
    o = ring_bridging_oxygens(cycle, bridge)
    if o is None: return None
    d = [pbc_dist(pos_o[o[k]], pos_o[o[(k+3) % 6]], box) for k in range(3)]
    return 0.5 * (max(d) - min(d))

def distortion_8MR(cycle, bridge, pos_o, box):
    """δ_8 = (1/2) max(|d_04 - d_26|, |d_15 - d_37|)"""
    o = ring_bridging_oxygens(cycle, bridge)
    if o is None: return None
    d = [pbc_dist(pos_o[o[k]], pos_o[o[(k+4) % 8]], box) for k in range(4)]
    return 0.5 * max(abs(d[0] - d[2]), abs(d[1] - d[3]))

DIST_FUNC = {4: distortion_4MR, 6: distortion_6MR, 8: distortion_8MR}

# -----------------------------------------------------------------------
# 5. Main: analyse a single structure
# -----------------------------------------------------------------------
def analyse_structure(pos_si, pos_o, box, lengths=(4, 6, 8), verbose=False):
    adj, bridge = build_graph(pos_si, pos_o, box, cutoff=1.95)
    degrees = [len(adj[i]) for i in range(len(pos_si))]
    if verbose:
        print(f"  Graph: avg deg = {np.mean(degrees):.2f}, min = {min(degrees)}, max = {max(degrees)}")
    rings = find_all_minimal_rings(adj, pos_si, box, lengths=lengths)
    out = {}
    for L in lengths:
        deltas = []
        for c in rings[L]:
            d = DIST_FUNC[L](c, bridge, pos_o, box)
            if d is not None:
                deltas.append(d)
        out[L] = dict(n_rings=len(rings[L]), deltas=np.array(deltas))
        if verbose:
            print(f"  {L}MR: {len(rings[L])} rings, mean δ = {np.mean(deltas) if deltas else float('nan'):.4f} Å, "
                  f"max = {max(deltas) if deltas else float('nan'):.4f} Å")
    return out

# -----------------------------------------------------------------------
# 6. CLI
# -----------------------------------------------------------------------
def load_structure_from_cif(fname):
    cell, atoms = parse_cif_cubic(fname)
    a = cell['a']
    si_frac = [(t[1], t[2], t[3]) for t in atoms if t[0].lower().startswith('si')]
    o_frac  = [(t[1], t[2], t[3]) for t in atoms if t[0].lower().startswith('o')]
    pos_si = np.array(si_frac) * a
    pos_o  = np.array(o_frac) * a
    return pos_si, pos_o, np.array([a, a, a])

def load_structure_from_data(fname):
    atoms, box = parse_lammps_data(fname)
    pos_si = np.array([[a[2], a[3], a[4]] for a in atoms if a[1] == 1])
    pos_o  = np.array([[a[2], a[3], a[4]] for a in atoms if a[1] == 2])
    return pos_si, pos_o, box

if __name__ == '__main__':
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument('input', help='LAMMPS .data file or .cif')
    p.add_argument('--lengths', '-L', default='4,6,8',
                   help='Ring lengths to enumerate (default: 4,6,8)')
    p.add_argument('--detail', action='store_true',
                   help='Print per-ring distortion values')
    args = p.parse_args()
    lengths = tuple(int(x) for x in args.lengths.split(','))
    if args.input.endswith('.data'):
        pos_si, pos_o, box = load_structure_from_data(args.input)
    else:
        pos_si, pos_o, box = load_structure_from_cif(args.input)
    print(f"Loaded {len(pos_si)} Si, {len(pos_o)} O from {args.input}")
    print(f"Box = {box}")
    out = analyse_structure(pos_si, pos_o, box, lengths=lengths, verbose=True)
    if args.detail:
        for L in lengths:
            print(f"\n  Per-ring δ_{L} [Å]:")
            for v in sorted(out[L]['deltas']):
                print(f"    {v:.4f}")
