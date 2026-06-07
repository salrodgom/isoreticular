#!/usr/bin/env python3
"""
Tabulate the dimensionless Landau collapse of the isoreticular RHO family.

After substituting the closed-form parametrisations of Table tab:SI-landau-fits
of the manuscript into the sextic functional, the effective Landau density
reduces to a single-parameter polynomial in the reduced variables
    pi  = p / p_c(k) - 1     (reduced pressure)
    phi = Delta / delta(k)   (reduced order parameter)
    fbar(phi; pi, eta) = 2 f_eff / (alpha delta^2(k))
                       = - pi phi^2 + 1/2 phi^4 + (eta/12) phi^6.
The crossover parameter eta(k) = 4 w alpha / u_eff^2(k) is the only relevant
control parameter once the framework-specific scales p_c(k) and delta(k) have
been absorbed into the axes.

Writes:
    data/figS_landau_collapse_phi_eq.dat
        columns: pi, phi_eq(pi; eta_k) for k=1..5 and the mean-field limit eta=0
    data/figS_landau_collapse_wells.dat
        columns: phi, fbar(phi; pi=0.5, eta_k) for k=1..5

The broken-phase minimum is the closed form
    phi^2_eq(pi; eta) = 2 (sqrt(1 + eta pi) - 1) / eta   (eta > 0)
                     =  pi                                (eta = 0, mean-field)
"""
import os
import numpy as np

# eta(k) from the DIRECT sextic fit of Delta(p) (tab:SI-eta-direct of the SI).
# For G_4 and G_5 the sextic fit is consistent with eta = 0 within precision,
# so we use a tiny value here purely to keep the closed-form phi_eq numerically
# stable at the rendering grid.
ETA = {1: 5.32, 2: 0.12, 3: 0.05, 4: 0.01, 5: 0.01}


def phi2_eq(pi, eta):
    """Broken-phase minimum, vectorised in pi."""
    pi = np.asarray(pi, dtype=float)
    out = np.zeros_like(pi)
    mask = pi > 0
    if eta == 0:
        out[mask] = pi[mask]
    else:
        out[mask] = 2.0 * (np.sqrt(1.0 + eta * pi[mask]) - 1.0) / eta
    return out


def fbar(phi, pi, eta):
    """Dimensionless effective free energy."""
    return -pi * phi**2 + 0.5 * phi**4 + (eta / 12.0) * phi**6


def main():
    here = os.path.dirname(os.path.abspath(__file__))
    out_dir = os.path.normpath(os.path.join(here, '..', 'data'))
    os.makedirs(out_dir, exist_ok=True)

    # ------------------------------------------------------------------
    # (a) Order-parameter trajectory phi_eq(pi) for the five eta values.
    # Logarithmic pi grid from 1e-3 to 1e1.
    # ------------------------------------------------------------------
    pi = np.logspace(-3, 1, 401)
    cols = [pi]
    header_a = (
        "Order-parameter trajectory in the dimensionless collapse.\n"
        "Columns: pi  phi_eq(eta=0)  phi_eq(eta_1)  phi_eq(eta_2)  phi_eq(eta_3) "
        " phi_eq(eta_4)  phi_eq(eta_5)\n"
        "eta_k = 0, 160, 6, 6, 1, 1 (mean-field reference and G_1..G_5)."
    )
    # Mean-field reference (eta = 0)
    cols.append(np.sqrt(phi2_eq(pi, 0.0)))
    for k in sorted(ETA):
        cols.append(np.sqrt(phi2_eq(pi, ETA[k])))
    arr_a = np.column_stack(cols)
    out_a = os.path.join(out_dir, 'figS_landau_collapse_phi_eq.dat')
    np.savetxt(out_a, arr_a, header=header_a, fmt='%.6e')
    print(f"Wrote {out_a}")

    # ------------------------------------------------------------------
    # (b) Wells fbar(phi; pi=0.5, eta) for the five eta values.
    # ------------------------------------------------------------------
    pi_w = 0.5
    phi = np.linspace(0.0, 1.6, 601)
    cols = [phi]
    header_b = (
        f"Dimensionless effective Landau wells at fixed pi = p/p_c - 1 = {pi_w}.\n"
        "Columns: phi  fbar(eta_1)  fbar(eta_2)  fbar(eta_3)  fbar(eta_4)"
        "  fbar(eta_5)  fbar(eta=0, mean-field reference)\n"
        "eta_k = 160, 6, 6, 1, 1, 0 (G_1..G_5 + mean-field reference)."
    )
    for k in sorted(ETA):
        cols.append(fbar(phi, pi_w, ETA[k]))
    cols.append(fbar(phi, pi_w, 0.0))
    arr_b = np.column_stack(cols)
    out_b = os.path.join(out_dir, 'figS_landau_collapse_wells.dat')
    np.savetxt(out_b, arr_b, header=header_b, fmt='%.6e')
    print(f"Wrote {out_b}")

    # ------------------------------------------------------------------
    # Sanity check: write the minima of the wells for the markers in (b).
    # ------------------------------------------------------------------
    rows = []
    for k in sorted(ETA):
        phi_eq = float(np.sqrt(phi2_eq(np.array([pi_w]), ETA[k]))[0])
        rows.append((ETA[k], phi_eq, float(fbar(phi_eq, pi_w, ETA[k]))))
    out_min = os.path.join(out_dir, 'figS_landau_collapse_wells_minima.dat')
    np.savetxt(out_min, np.asarray(rows),
               header=('Broken-phase minima of the dimensionless wells at pi=0.5.\n'
                       'Columns: eta_k, phi_eq, fbar_eq'),
               fmt='%.6e')
    print(f"Wrote {out_min}")


if __name__ == '__main__':
    main()
