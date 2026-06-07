#!/usr/bin/env python3
"""
Prepare tabulated inputs for Figure 3 of the manuscript:

  Figure 3(a)  V/V_0 vs hydrostatic pressure for G_1 to G_5
  Figure 3(b)  D8R distortion Delta vs hydrostatic pressure for G_1 to G_5

Reads the raw PLUMED/GULP scans collected in

    initial_structures_RHO_isoreticular/dir_RHO_isoreticular_G<k>_*/data_pressure_delta.txt

(column conventions are documented in extract.sh of that repository;
relevant columns here are col 2 = pressure in bar, col 4 = D8R distortion
Delta in nm, col 5 = cubic-cell volume V in nm^3, col 9 = signed soft-mode
frequency, used as a stable/unstable classifier of the parent branch).

For each member k = 1..5 the script writes four .dat files in data/:

    fig3_G<k>_volume_stable.dat     p [GPa]   V/V_0 [-]    (omega1 >= 0)
    fig3_G<k>_volume_unstable.dat   p [GPa]   V/V_0 [-]    (omega1 <  0)
    fig3_G<k>_delta_stable.dat      p [GPa]   Delta [A]    (omega1 >= 0)
    fig3_G<k>_delta_unstable.dat    p [GPa]   Delta [A]    (omega1 <  0)

V_0 is the V (in nm^3) reported at the first row of each input scan
(p = 0 bar) for that member. G_1 uses the 2x2x2 supercell (384 T-atoms);
all other members use their crystallographic unit cell (672, 1440 and
2640 T-atoms for G_3, G_4 and G_5; G_2 uses 1x1x1 of PWN, 240 T-atoms).
Normalisation per cell does not matter for v = V/V_0.

The script writes nothing else and produces no plot.
"""
import os

# Resolve the path to the raw GULP/PLUMED scans relative to this script.
# scripts/   -> manuscript_figures/   -> isoreticular_RHO/   -> initial_structures...
_HERE = os.path.dirname(os.path.abspath(__file__))
BASE = os.environ.get(
    'RHO_RAW_DATA',
    os.path.normpath(os.path.join(_HERE, '..', '..',
                                  'initial_structures_RHO_isoreticular')))

INPUTS = {
    1: os.path.join(BASE, 'dir_RHO_isoreticular_G1_222_SG_P1',
                    'data_pressure_delta.txt'),
    2: os.path.join(BASE, 'dir_RHO_isoreticular_G2_SG_P1',
                    'data_pressure_delta.txt'),
    3: os.path.join(BASE, 'dir_RHO_isoreticular_G3_SG_P1',
                    'data_pressure_delta.txt'),
    4: os.path.join(BASE, 'dir_RHO_isoreticular_G4_SG_P1',
                    'data_pressure_delta.txt'),
    5: os.path.join(BASE, 'dir_RHO_isoreticular_G5_SG_P1',
                    'data_pressure_delta.txt'),
}

# G_5 has no col-9 signed soft-mode column populated; treat all its points
# as "stable" for plotting so that the same script handles the whole family
# uniformly.
ALWAYS_STABLE = {5}

OUT_DIR = os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'data')


def parse(path):
    rows = []
    with open(path) as f:
        for line in f:
            parts = line.split()
            if len(parts) < 5:
                continue
            try:
                p_bar = float(parts[1])
                Delta_nm = float(parts[3])
                V_nm3 = float(parts[4])
            except (ValueError, IndexError):
                continue
            omega_signed = 0.0
            if len(parts) >= 9:
                try:
                    omega_signed = float(parts[8])
                except ValueError:
                    omega_signed = 0.0
            rows.append((p_bar / 1e4, Delta_nm * 10.0, V_nm3, omega_signed))
    return rows


def write_dat(path, header, rows):
    with open(path, 'w') as f:
        f.write(header)
        for p, y in rows:
            f.write(f'{p:.6f} {y:.6f}\n')


def main():
    os.makedirs(OUT_DIR, exist_ok=True)
    for k, path in INPUTS.items():
        rows = parse(path)
        if not rows:
            print(f'G_{k}: no rows in {path}')
            continue

        # V_0 = volume at the first p = 0 bar row of the scan
        V0 = next((V for p, _, V, _ in rows if abs(p) < 1e-9), rows[0][2])

        always_stable = (k in ALWAYS_STABLE)

        vol_stable, vol_unstable = [], []
        del_stable, del_unstable = [], []
        for p, Delta_A, V_nm3, w in rows:
            v = V_nm3 / V0
            is_stable = always_stable or (w >= 0.0)
            (vol_stable if is_stable else vol_unstable).append((p, v))
            (del_stable if is_stable else del_unstable).append((p, Delta_A))

        prefix = os.path.join(OUT_DIR, f'fig3_G{k}')
        write_dat(f'{prefix}_volume_stable.dat',
                  '# G_{} volume scan (omega1 >= 0)\n'
                  '# Columns: p [GPa]   V/V_0 [-]\n'.format(k),
                  vol_stable)
        write_dat(f'{prefix}_volume_unstable.dat',
                  '# G_{} volume scan (omega1 < 0; parent metastable)\n'
                  '# Columns: p [GPa]   V/V_0 [-]\n'.format(k),
                  vol_unstable)
        write_dat(f'{prefix}_delta_stable.dat',
                  '# G_{} D8R distortion (omega1 >= 0)\n'
                  '# Columns: p [GPa]   Delta [A]\n'.format(k),
                  del_stable)
        write_dat(f'{prefix}_delta_unstable.dat',
                  '# G_{} D8R distortion (omega1 < 0; parent metastable)\n'
                  '# Columns: p [GPa]   Delta [A]\n'.format(k),
                  del_unstable)

        print(f'G_{k}: V_0 = {V0:.6f} nm^3, '
              f'rows = {len(rows)} '
              f'(stable {len(vol_stable)}, unstable {len(vol_unstable)})')


if __name__ == '__main__':
    main()
