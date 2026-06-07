#!/usr/bin/env python3
"""
Tabulate the *real-parameter* effective Landau wells for G_1 and G_5
(panels c and d of the merged Figure S2 of the manuscript).

Reads nothing; writes:
    data/figS2_landau_real_G1.dat        columns: Delta, f/alpha at p = 0.5 p_c, p_c, 1.5 p_c
    data/figS2_landau_real_G1_minima.dat columns: p_factor, Delta_eq, f_eq
    data/figS2_landau_real_G5.dat        idem for G_5
    data/figS2_landau_real_G5_minima.dat

Functional:
    f_eff(Delta; p) / A = (1/2) (1 - p/p_c) Delta^2 + (1/4) (U_eff/A) Delta^4
with parameters taken from Table tab:SI-landau of the manuscript:
        delta^2 [A^2]   U_eff/A [A^-2]   p_c [GPa]
    G_1     4.884           0.205            0.9418
    G_5     0.215           4.640            0.1282

Plot range: Delta in [0, 1.5 * delta_eq]; three pressures 0.5 p_c, p_c, 1.5 p_c.
"""
import os
import numpy as np

LANDAU = {
    1: dict(label='G_1',  delta2=4.884, u_over_a=0.205, pc=0.9418),
    5: dict(label='G_5',  delta2=0.215, u_over_a=4.640, pc=0.1282),
}
P_FACTORS = (0.5, 1.0, 1.5)


def f_eff_real(D, p, pc, u_over_a):
    return 0.5 * (1.0 - p / pc) * D**2 + 0.25 * u_over_a * D**4


def main():
    here = os.path.dirname(os.path.abspath(__file__))
    out_dir = os.path.normpath(os.path.join(here, '..', 'data'))
    os.makedirs(out_dir, exist_ok=True)

    # G_1 needs Delta up to 2.5 A (its delta = 2.21 A); G_5 lives in
    # [0, 0.8] A. Sample G_1 on [0, 2.5] and G_5 on [0, 0.8] with the same
    # density so each panel is well resolved on its own axis.
    D_RANGE = {1: (0.0, 2.5), 5: (0.0, 0.8)}
    for k, L in LANDAU.items():
        D = np.linspace(D_RANGE[k][0], D_RANGE[k][1], 600)
        delta_eq = float(np.sqrt(L['delta2']))
        cols = [D]
        minima = []
        for pf in P_FACTORS:
            p = pf * L['pc']
            cols.append(f_eff_real(D, p, L['pc'], L['u_over_a']))
            if p > L['pc']:
                D_eq = float(np.sqrt((p / L['pc'] - 1.0) / L['u_over_a']))
                y_eq = float(f_eff_real(D_eq, p, L['pc'], L['u_over_a']))
                minima.append((pf, D_eq, y_eq))
        header = (f"f_eff(Delta;p)/alpha for {L['label']} at three pressures "
                  f"(0.5 p_c, p_c, 1.5 p_c).\n"
                  f"p_c = {L['pc']} GPa, delta = {delta_eq:.4f} A, "
                  f"u_eff/alpha = {L['u_over_a']} A^-2.\n"
                  "Columns: Delta [A]  f/alpha [A^2] at 0.5 p_c  at p_c  at 1.5 p_c")
        np.savetxt(os.path.join(out_dir, f'figS2_landau_real_G{k}.dat'),
                   np.column_stack(cols), header=header, fmt='%.6f')
        np.savetxt(os.path.join(out_dir, f'figS2_landau_real_G{k}_minima.dat'),
                   np.asarray(minima) if minima else np.empty((0, 3)),
                   header='Broken-phase minima. Columns: p/p_c, Delta_eq [A], f_eq/alpha [A^2]',
                   fmt='%.6f')
        print(f"Wrote G_{k} files (delta_eq = {delta_eq:.3f} A)")


if __name__ == '__main__':
    main()
