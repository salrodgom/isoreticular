#!/usr/bin/env python3
"""
Extract the framework-density vs cohesive-energy data of Figure S6:
all SLC-optimised structures of G_1 ... G_5 over the available pressure scan.

Reads:
    ../../initial_structures_RHO_isoreticular/dir_RHO_isoreticular_G{k}_*/
        data_pressure_delta.txt

Writes 6 .dat files in data/:

    figS6_density_energy_G{k}.dat   for k = 1..5
        columns: p [GPa]   rho [T-atoms / 10^3 A^3]   E_per_Si [eV/Si]
    figS6_density_energy_equilibrium.dat
        columns: rho_eq [T/nm3]   E_eq [eV/Si]   k (integer)
        one row per G_k, the lowest-pressure point of each scan.

The Balestra et al. (CGD 2024) reference parabola is hard-coded in the
gnuplot script; only the anchor (E_quartz from G_1) is written to a
small text file for convenience.

    figS6_density_energy_meta.txt
        single line: E_quartz [eV/Si] anchored at G_1, plus k_f and FD_0
        of the published parabola.

Run:
    python3 scripts/prep_figS6_density_energy.py
"""
import os
import numpy as np

HERE_PREP = os.path.dirname(os.path.abspath(__file__))
RHO_ROOT  = os.path.normpath(os.path.join(HERE_PREP, '..', '..',
                                          'initial_structures_RHO_isoreticular'))

# G_1 uses the 2x2x2 supercell (384 = 48*8); G_2..G_5 use the unit cell.
SLC_DIR = {
    1: 'dir_RHO_isoreticular_G1_222_SG_P1',
    2: 'dir_RHO_isoreticular_G2_SG_P1',
    3: 'dir_RHO_isoreticular_G3_SG_P1',
    4: 'dir_RHO_isoreticular_G4_SG_P1',
    5: 'dir_RHO_isoreticular_G5_SG_P1',
}
N_T = {1: 384, 2: 240, 3: 672, 4: 1440, 5: 2640}

# Balestra et al. CGD 2024 parabola E(FD) = E_quartz + k_f (FD - FD_0)^2.
KF_BAL2024  = 1.36529029359175e-3   # eV per T per (T/10^3 A^3)^2
FD0_BAL2024 = 27.735093             # T-atoms / (10^3 A^3) (parabola vertex)


def load_one(k):
    path = os.path.join(RHO_ROOT, SLC_DIR[k], 'data_pressure_delta.txt')
    if not os.path.exists(path):
        return None
    rows = []
    with open(path) as f:
        for line in f:
            parts = line.split()
            if len(parts) < 5:
                continue
            try:
                p_bar = float(parts[1])
                E_per_Si = float(parts[2])
                V_nm3 = float(parts[4])
            except (ValueError, IndexError):
                continue
            if V_nm3 <= 0:
                continue
            p_GPa = p_bar / 10000.0
            rho = N_T[k] / V_nm3
            rows.append((p_GPa, rho, E_per_Si))
    if not rows:
        return None
    arr = np.array(rows)
    _, idx = np.unique(np.round(arr, 8), axis=0, return_index=True)
    return arr[np.sort(idx)]


def main():
    out_dir = os.path.normpath(os.path.join(HERE_PREP, '..', 'data'))
    os.makedirs(out_dir, exist_ok=True)

    equil_rows = []
    for k in sorted(SLC_DIR):
        arr = load_one(k)
        if arr is None:
            print(f"  G_{k}: data file missing, skipping")
            continue
        fname = f'figS6_density_energy_G{k}.dat'
        np.savetxt(os.path.join(out_dir, fname), arr,
                   header=f'G_{k} pressure scan (SLC).\n'
                          'Columns: p [GPa]   rho [T/10^3 A^3]   E_per_Si [eV/Si]',
                   fmt='%.6f')
        j = int(np.argmin(arr[:, 0]))
        equil_rows.append((arr[j, 1], arr[j, 2], k))
        print(f'  G_{k}: {len(arr)} pts  rho in [{arr[:,1].min():.2f}, '
              f'{arr[:,1].max():.2f}]  p in [{arr[:,0].min():.2f}, '
              f'{arr[:,0].max():.2f}] GPa')

    # Equilibrium points
    np.savetxt(os.path.join(out_dir, 'figS6_density_energy_equilibrium.dat'),
               np.asarray(equil_rows),
               header='Equilibrium (p->0) point per G_k.\n'
                      'Columns: rho_eq [T/10^3 A^3]   E_eq [eV/Si]   k',
               fmt='%.6f')

    # Anchor the Balestra parabola at G_1 and dump metadata
    if equil_rows:
        rho_G1, E_G1, _ = equil_rows[0]
        E_quartz = E_G1 - KF_BAL2024 * (rho_G1 - FD0_BAL2024) ** 2
        meta_path = os.path.join(out_dir, 'figS6_density_energy_meta.txt')
        with open(meta_path, 'w') as f:
            f.write(f'# Balestra et al. CGD 2024 parabola anchored at G_1.\n')
            f.write(f'E_quartz   = {E_quartz:.6f}   # eV/Si\n')
            f.write(f'k_f        = {KF_BAL2024:.10e}   # eV per (T/10^3 A^3)^2\n')
            f.write(f'FD_0       = {FD0_BAL2024:.6f}   # T-atoms / 10^3 A^3\n')
        print(f"  E_quartz (anchor at G_1) = {E_quartz:.4f} eV/Si")


if __name__ == '__main__':
    main()
