#!/usr/bin/env python3
"""
Tabulate the dimensionless conceptual Landau wells of Figure S2.

Reads nothing; writes:
    data/figS2_landau_concept_a.dat   columns: Delta, f(p_red=-0.6, -0.2, 0, +0.2, +0.6)
    data/figS2_landau_concept_a_minima.dat   columns: p_red, Delta_eq, f_eq
    data/figS2_landau_concept_b.dat   columns: v, df(Delta=0.0, 0.4, 0.8, 1.2)
    data/figS2_landau_concept_b_minima.dat   columns: Delta, v_min, df(v_min)

Sign convention of the manuscript:

    f_eff(Delta; p) = (1/2) alpha (1 - p/p_c) Delta^2 + (1/4) u Delta^4
    f(Delta, v; p) = f_eff + (1/2) K_0 (v - 1)^2 + g Delta^2 (v - 1)

so the broken phase appears for p > p_c (p_red > 0). Δ >= 0 and v <= 1
are enforced by the plotting ranges, not here.

Run from the repo root:
    python3 scripts/prep_figS2_landau_wells_concept.py
"""
import os
import numpy as np

# Conceptual (dimensionless) coefficients --- arbitrary units.
ALPHA = 1.0
U     = 1.0
K0    = 1.0
G     = 0.25

# Panel (a): five pressures, Delta on [0, 1.7]
P_REDS = (-0.6, -0.2, 0.0, +0.2, +0.6)
DELTA  = np.linspace(0.0, 1.7, 600)

def f_eff(D, pr):
    """Effective Landau in Delta alone, dimensionless."""
    return 0.5 * ALPHA * (-pr) * D**2 + 0.25 * U * D**4

# Panel (b): four Deltas, v on [0.5, 1.0], p_red = +0.3
DELTAS_B = (0.0, 0.4, 0.8, 1.2)
P_RED_B  = +0.3
V_AXIS   = np.linspace(0.5, 1.0, 600)

def f_full(D, v, pr):
    return (f_eff(D, pr)
            + 0.5 * K0 * (v - 1.0)**2
            + G * D**2 * (v - 1.0))

def v_min_of_Delta(D):
    """Argmin of f wrt v at fixed Delta: v_min - 1 = -(g/K_0) Delta^2."""
    return 1.0 - G * D**2 / K0


def main():
    here = os.path.dirname(os.path.abspath(__file__))
    out_dir = os.path.normpath(os.path.join(here, '..', 'data'))
    os.makedirs(out_dir, exist_ok=True)

    # ---------- panel (a) ----------
    cols_a = [DELTA] + [f_eff(DELTA, pr) for pr in P_REDS]
    header_a = ('Conceptual f_eff(Delta; p_red) for the manuscript Landau '
                'functional with alpha = u = 1.\n'
                'Columns: Delta  f_p=-0.6  f_p=-0.2  f_p=0  f_p=+0.2  f_p=+0.6')
    np.savetxt(os.path.join(out_dir, 'figS2_landau_concept_a.dat'),
               np.column_stack(cols_a),
               header=header_a, fmt='%.6f')

    minima_a = []
    for pr in P_REDS:
        if pr > 0:
            D_eq = float(np.sqrt(ALPHA * pr / U))
            minima_a.append((pr, D_eq, float(f_eff(D_eq, pr))))
    np.savetxt(os.path.join(out_dir, 'figS2_landau_concept_a_minima.dat'),
               np.asarray(minima_a),
               header='Broken-phase minima for p_red > 0. Columns: p_red, Delta_eq, f_eq',
               fmt='%.6f')

    # ---------- panel (b) ----------
    refs = [f_full(D, 1.0, P_RED_B) for D in DELTAS_B]
    cols_b = [V_AXIS] + [f_full(D, V_AXIS, P_RED_B) - r
                         for D, r in zip(DELTAS_B, refs)]
    header_b = ('Conceptual f(Delta, v; p_red=+0.3) - f(Delta, 1; p_red=+0.3) '
                'for four fixed Delta values.\n'
                'Columns: v  df_D=0.0  df_D=0.4  df_D=0.8  df_D=1.2')
    np.savetxt(os.path.join(out_dir, 'figS2_landau_concept_b.dat'),
               np.column_stack(cols_b),
               header=header_b, fmt='%.6f')

    minima_b = []
    for D, r in zip(DELTAS_B, refs):
        vm = float(v_min_of_Delta(D))
        if vm >= V_AXIS.min():
            minima_b.append((D, vm, float(f_full(D, vm, P_RED_B)) - r))
    np.savetxt(os.path.join(out_dir, 'figS2_landau_concept_b_minima.dat'),
               np.asarray(minima_b) if minima_b else np.empty((0, 3)),
               header='v_min markers per Delta. Columns: Delta, v_min, df(v_min)',
               fmt='%.6f')

    print(f"Wrote 4 files to {out_dir}/")


if __name__ == '__main__':
    main()
