#!/usr/bin/env python3
"""
Topological asymmetry diagnostic for the sampled 8MRs.

For each G_k topology, this script does NOT enumerate all rings (that
requires a heavier topology toolkit). Instead, it samples the LOCAL ENVIRONMENT
of the 8MR atoms used by the order-parameter CV and tests whether the two
8MRs of the same D8R are crystallographically equivalent (symmetric D8R) or
not (asymmetric D8R, i.e., between two different cage types).

Method:
  1. Read the high-symmetry .cif for each G_k (the parent Im-3m structure).
  2. Extract Si positions.
  3. For each Si, count Si neighbours within radius shells [3.0-3.5, 3.5-4.5,
     4.5-5.5, 5.5-7.0] Å. The cumulative coordination number signature
     S(Si) = (n_1, n_2, n_3, n_4) fingerprints the local cage environment
     of each T-atom.
  4. Compute the diversity of S over the unit cell (number of distinct
     fingerprints). G_1 (RHO) should have few; PAU should have many.
  5. As a more direct test: identify D8Rs by clustering Si into rings is
     out of scope here, but the diversity of S already correlates with the
     number of inequivalent cage environments.

Output:
  - For each G_k: number of distinct local environments and their
    populations, indicating the topological complexity that underlies
    the per-8MR bimodality observed in histograms_8MR_per_pressure.py.
"""
import os, sys, glob
import numpy as np

ROOT = '/sessions/wizardly-eager-fermat/mnt/initial_structures_RHO_isoreticular'

CIFS = [
    ('G_1', 'dir_RHO_isoreticular_G1_SG_P1/RHO_isoreticular_G1_SG_P1.cif'),
    ('G_2', 'dir_RHO_isoreticular_G2_SG_P1/RHO_isoreticular_G2_SG_P1_topol.cif'),
    ('G_3', 'dir_RHO_isoreticular_G3_SG_P1/RHO_isoreticular_G3_SG_P1_topol.cif'),
    ('G_4', 'dir_RHO_isoreticular_G4_SG_P1/RHO_isoreticular_G4_SG_P1_topol.cif'),
    ('G_5', 'dir_RHO_isoreticular_G5_SG_P1/RHO_isoreticular_G5_SG_P1_topol.cif'),
]

# Coordination shells (Å). Si-O-Si bridges put first-neighbour Si at ~3.1 Å.
SHELLS = [(3.0, 3.5), (3.5, 4.5), (4.5, 5.5), (5.5, 7.0)]

def parse_cif(fname):
    """Minimal CIF parser: cell parameters and Si fractional positions."""
    cell = {}; atoms = []
    in_loop = False; cols = []
    with open(fname) as f:
        for line in f:
            line = line.rstrip()
            if line.startswith('_cell_length_a'):  cell['a'] = float(line.split()[1])
            elif line.startswith('_cell_length_b'): cell['b'] = float(line.split()[1])
            elif line.startswith('_cell_length_c'): cell['c'] = float(line.split()[1])
            elif line.startswith('_cell_angle_alpha'): cell['al'] = float(line.split()[1])
            elif line.startswith('_cell_angle_beta'):  cell['be'] = float(line.split()[1])
            elif line.startswith('_cell_angle_gamma'): cell['ga'] = float(line.split()[1])
            elif line.strip() == 'loop_':
                in_loop = True; cols = []
            elif in_loop and line.startswith('_'):
                cols.append(line.strip())
            elif in_loop and line.strip() and not line.startswith('_'):
                parts = line.split()
                if len(cols) >= 5 and 'atom_site' in cols[0]:
                    # Find indices of x, y, z, label
                    try:
                        ix = next(i for i, c in enumerate(cols) if 'fract_x' in c)
                        iy = next(i for i, c in enumerate(cols) if 'fract_y' in c)
                        iz = next(i for i, c in enumerate(cols) if 'fract_z' in c)
                        il = next(i for i, c in enumerate(cols) if c.endswith('label')
                                  or c.endswith('type_symbol'))
                    except StopIteration:
                        continue
                    if len(parts) > max(ix, iy, iz, il):
                        sym = parts[il].rstrip('0123456789')
                        atoms.append((sym, float(parts[ix]), float(parts[iy]), float(parts[iz])))
                else:
                    # Not an atom_site loop, stop reading
                    in_loop = False
    return cell, atoms

def min_image_distances(frac_positions, cell):
    """Compute all pairwise distances under PBC, cubic cell assumed."""
    a = cell['a']
    frac = np.array(frac_positions)
    # Cubic only (alpha=beta=gamma=90); generalise if needed
    cart = frac * a
    N = len(cart)
    D = np.zeros((N, N))
    for i in range(N):
        dvec = cart - cart[i]
        # Minimum image
        dvec -= a * np.round(dvec / a)
        D[i] = np.linalg.norm(dvec, axis=1)
    return D

def local_signature(D, shells):
    """For each atom, count number of neighbours in each shell."""
    N = D.shape[0]
    sig = np.zeros((N, len(shells)), dtype=int)
    for k, (rmin, rmax) in enumerate(shells):
        sig[:, k] = ((D > rmin) & (D < rmax)).sum(axis=1)
    return sig

print(f"{'system':6s} {'N_Si':6s} {'cell a':>8s} {'unique env':>11s}  populations (count: signature)")
print('-' * 110)

for name, relpath in CIFS:
    fpath = os.path.join(ROOT, relpath)
    if not os.path.exists(fpath):
        # Try without _topol suffix
        alt = fpath.replace('_topol', '')
        if os.path.exists(alt):
            fpath = alt
        else:
            print(f"{name}: no .cif found at {relpath}")
            continue
    try:
        cell, atoms = parse_cif(fpath)
    except Exception as e:
        print(f"{name}: parse error: {e}")
        continue
    if not atoms or 'a' not in cell:
        print(f"{name}: empty or missing cell")
        continue
    # Filter Si only
    si_atoms = [(a[1], a[2], a[3]) for a in atoms if a[0].startswith('Si') or a[0] == 'T']
    if not si_atoms:
        si_atoms = [(a[1], a[2], a[3]) for a in atoms if not a[0].startswith('O')]
    n_si = len(si_atoms)
    if n_si == 0:
        print(f"{name}: no Si atoms found")
        continue
    D = min_image_distances(si_atoms, cell)
    sig = local_signature(D, SHELLS)
    # Convert each signature row to tuple for hashing
    sig_tuples = [tuple(row) for row in sig]
    unique, counts = np.unique(sig_tuples, axis=0, return_counts=True)
    # Display compact
    pops = sorted(zip(counts, [tuple(u) for u in unique]), reverse=True)
    pops_str = ', '.join(f"{c}: {s}" for c, s in pops[:5])
    if len(pops) > 5:
        pops_str += f", ... ({len(pops)} total)"
    print(f"{name:6s} {n_si:>6d} {cell['a']:>8.3f} {len(unique):>11d}  {pops_str}")
