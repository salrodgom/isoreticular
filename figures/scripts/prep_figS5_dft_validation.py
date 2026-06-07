#!/usr/bin/env python3
"""
Extract the validation-of-SLC datasets for Figure S5 of the manuscript:
G_1 and G_2 cubic frameworks compared across four / two levels of theory.

Reads CIF/CONTCAR files under RHO_MLFF_DFT/ and the SLC scan output, then
writes 6 .dat files:

    data/figS5_g1_dft.dat        G_1 DFT r2SCAN+rVV10
    data/figS5_g1_slc.dat        G_1 SLC (this work)
    data/figS5_g1_nasir.dat      G_1 MACE (Nasir et al.)
    data/figS5_g1_matpes.dat     G_1 MACE-MP foundation (MatPES r2SCAN)
    data/figS5_g2_dft.dat        G_2 DFT r2SCAN+rVV10
    data/figS5_g2_slc.dat        G_2 SLC (this work)

Each file: 3 columns, p [GPa]  V/V_0 [-]  Delta [A]

The Delta extractor and CIF/POSCAR parsers are imported from the existing
RHO_MLFF_DFT/fig3_slc_dft_richard.py to avoid duplication.
"""
import os
import sys
import glob
import re as _re
import numpy as np

HERE_PREP = os.path.dirname(os.path.abspath(__file__))
# Resolve the simulations root by walking up from this script (manuscript_figures/scripts/);
# initial_structures_RHO_isoreticular/ sits next to manuscript_figures/.
RHO_ROOT  = os.path.normpath(os.path.join(HERE_PREP, '..', '..',
                                          'initial_structures_RHO_isoreticular'))
MLFF      = os.path.join(RHO_ROOT, 'RHO_MLFF_DFT')

# Import the loaders from the existing python script.
sys.path.insert(0, MLFF)
sys.path.insert(0, os.path.join(RHO_ROOT, 'scripts_FES'))
# fig3_slc_dft_richard.py imports ring_distortions at module level and
# prints data while loading; we re-use only its helper functions, so we
# need a clean import path.
from ring_distortions import load_structure_from_cif, analyse_structure  # noqa: E402


def cif_to_delta_V(fname):
    pos_si, pos_o, box = load_structure_from_cif(fname)
    V = float(box[0]**3)
    if len(pos_si) == 0 or len(pos_o) == 0:
        return V, None
    out = analyse_structure(pos_si, pos_o, box, lengths=(8,), verbose=False)
    if not out[8]['deltas'].size:
        return V, None
    return V, float(out[8]['deltas'].mean())


def parse_contcar(fname):
    with open(fname) as f:
        lines = f.readlines()
    scale = float(lines[1].split()[0])
    vecs = np.array([[float(x) for x in lines[i].split()] for i in (2, 3, 4)]) * scale
    a_cubic = float(np.mean(np.diag(vecs)))
    syms = lines[5].split()
    counts = [int(x) for x in lines[6].split()]
    mode = lines[7].strip().lower()
    coords = []
    nat = sum(counts)
    for i in range(8, 8 + nat):
        coords.append([float(x) for x in lines[i].split()[:3]])
    coords = np.array(coords)
    if mode.startswith('d'):
        coords = coords @ vecs
    sym_per_atom = []
    for s, n in zip(syms, counts):
        sym_per_atom.extend([s] * n)
    pos_si = np.array([coords[i] for i, s in enumerate(sym_per_atom)
                       if s.lower().startswith('si')])
    pos_o  = np.array([coords[i] for i, s in enumerate(sym_per_atom)
                       if s.lower().startswith('o')])
    return pos_si, pos_o, a_cubic


def contcar_to_delta_V(fname):
    pos_si, pos_o, a = parse_contcar(fname)
    V = float(a**3)
    if len(pos_si) == 0 or len(pos_o) == 0:
        return V, None
    box = np.array([a, a, a])
    out = analyse_structure(pos_si, pos_o, box, lengths=(8,), verbose=False)
    if not out[8]['deltas'].size:
        return V, None
    return V, float(out[8]['deltas'].mean())


def load_slc(path, n_uc=1, p_max_GPa=2.05):
    by_p = {}
    with open(path) as f:
        for line in f:
            parts = line.split()
            if len(parts) < 5:
                continue
            try:
                p_bar = float(parts[1])
                E = float(parts[2])
                delta_nm = float(parts[3])
                V_nm3 = float(parts[4])
            except (ValueError, IndexError):
                continue
            p_GPa = p_bar / 10000.0
            if p_GPa < 0 or p_GPa > p_max_GPa:
                continue
            if p_bar not in by_p or E < by_p[p_bar][0]:
                by_p[p_bar] = (E, delta_nm * 10.0, V_nm3 * 1000.0 / n_uc)
    ps = sorted(by_p.keys())
    P = np.array(ps) / 10000.0
    delta = np.array([by_p[p][1] for p in ps])
    V = np.array([by_p[p][2] for p in ps])
    return P, V, delta


def load_cif_scan(label, pattern):
    files = sorted(glob.glob(os.path.join(MLFF, pattern)),
                   key=lambda f: float(_re.search(r'(\d+\.\d+)GPa', f).group(1)))
    P, V, D = [], [], []
    for fn in files:
        m = _re.search(r'(\d+\.\d+)GPa', fn)
        if not m:
            continue
        p_gpa = float(m.group(1))
        v, delta = cif_to_delta_V(fn)
        if delta is None:
            print(f"  {label}: skipping {os.path.basename(fn)} (no Delta)")
            continue
        P.append(p_gpa); V.append(v); D.append(delta)
    return np.array(P), np.array(V), np.array(D)


def load_dft_g2_from_contcars(root):
    P, V, D = [], [], []
    for sub in sorted(glob.glob(os.path.join(root, '*/'))):
        m = _re.search(r'dir_([\d.]+)_bar', os.path.basename(os.path.normpath(sub)))
        if not m:
            continue
        p_bar = float(m.group(1)); p_gpa = p_bar / 10000.0
        files = glob.glob(os.path.join(sub, 'G2_RHO_dir_*_CONTCAR_*'))
        if not files:
            continue
        v, d = contcar_to_delta_V(files[0])
        if d is None:
            continue
        P.append(p_gpa); V.append(v); D.append(d)
    order = np.argsort(P)
    return np.array(P)[order], np.array(V)[order], np.array(D)[order]


def normalise(p, v):
    v0 = v[np.argmin(np.abs(p))]
    return v / v0


def write_dat(fname, p, v_norm, d, label):
    arr = np.column_stack((p, v_norm, d))
    header = (f'{label}\n'
              'Columns: p [GPa]   V/V_0 [-]   Delta [A]')
    np.savetxt(fname, arr, header=header, fmt='%.6f')
    print(f'  Wrote {os.path.basename(fname)}  ({len(p)} rows, '
          f'p in [{p.min():.2f}, {p.max():.2f}] GPa)')


def main():
    out_dir = os.path.normpath(os.path.join(HERE_PREP, '..', 'data'))
    os.makedirs(out_dir, exist_ok=True)

    print('--- G_1 ---')
    p_slc, V_slc, d_slc = load_slc(os.path.join(RHO_ROOT,
        'dir_RHO_isoreticular_G1_222_SG_P1', 'data_pressure_delta.txt'), n_uc=8)
    write_dat(os.path.join(out_dir, 'figS5_g1_slc.dat'),
              p_slc, normalise(p_slc, V_slc), d_slc, 'G_1 SLC 2x2x2 supercell')

    p_dft, V_dft, d_dft = load_cif_scan('G1 DFT', 'dft/optimized_*.cif')
    write_dat(os.path.join(out_dir, 'figS5_g1_dft.dat'),
              p_dft, normalise(p_dft, V_dft), d_dft,
              'G_1 DFT r2SCAN+rVV10')

    p_rch, V_rch, d_rch = load_cif_scan('G1 Nasir MACE', 'richard/optimized_*.cif')
    mask = np.abs(p_rch - 0.70) > 1e-3
    p_rch, V_rch, d_rch = p_rch[mask], V_rch[mask], d_rch[mask]
    write_dat(os.path.join(out_dir, 'figS5_g1_nasir.dat'),
              p_rch, normalise(p_rch, V_rch), d_rch,
              'G_1 MACE (Nasir et al.); spurious p=0.70 GPa point dropped')

    p_mat, V_mat, d_mat = load_cif_scan('G1 MatPES MACE-MP',
        'mace-matpes-r2scan/optimized_rhoa_*.cif')
    write_dat(os.path.join(out_dir, 'figS5_g1_matpes.dat'),
              p_mat, normalise(p_mat, V_mat), d_mat,
              'G_1 MACE-MP foundation (MatPES r2SCAN)')

    print('--- G_2 ---')
    p_slc2, V_slc2, d_slc2 = load_slc(os.path.join(RHO_ROOT,
        'dir_RHO_isoreticular_G2_SG_P1', 'data_pressure_delta.txt'),
        n_uc=1, p_max_GPa=2.05)
    write_dat(os.path.join(out_dir, 'figS5_g2_slc.dat'),
              p_slc2, normalise(p_slc2, V_slc2), d_slc2, 'G_2 SLC unit cell')

    p_dft2, V_dft2, d_dft2 = load_dft_g2_from_contcars(
        os.path.join(MLFF, 'G2', 'G2_RHO_r2scan_CONTCAR_20260520'))
    write_dat(os.path.join(out_dir, 'figS5_g2_dft.dat'),
              p_dft2, normalise(p_dft2, V_dft2), d_dft2,
              'G_2 DFT r2SCAN+rVV10')


if __name__ == '__main__':
    main()
