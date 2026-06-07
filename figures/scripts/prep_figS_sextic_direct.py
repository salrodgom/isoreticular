#!/usr/bin/env python3
"""
Sextic-Landau DIRECT fit of Delta(p) on the broken phase for G_1-G_5.

Model:
    Delta^2(p) = C [sqrt(1 + eta (p/p_c - 1)) - 1]
where C = U_eff/(2W) [A^2], eta = 4 W A / U_eff^2 [-].
In A = 1 units:
    U_eff/A = 2/(C eta)
    W/A     = 1/(C^2 eta)

Writes:
    data/figS_sextic_direct.dat
    data/figS_sextic_direct_pars.dat
"""
import os
import numpy as np
from scipy.optimize import curve_fit

RAW = '../../initial_structures_RHO_isoreticular'

MEMBERS = [
    dict(label='G_1', file='dir_RHO_isoreticular_G1_222_SG_P1', pc=0.9418),
    dict(label='G_2', file='dir_RHO_isoreticular_G2_SG_P1',      pc=0.5338),
    dict(label='G_3', file='dir_RHO_isoreticular_G3_SG_P1',      pc=0.3627),
    dict(label='G_4', file='dir_RHO_isoreticular_G4_SG_P1',      pc=0.2241),
    dict(label='G_5', file='dir_RHO_isoreticular_G5_SG_P1',      pc=0.1282),
]
PMAX = 2.0


def read_broken_phase(path, pc):
    P, D = [], []
    for line in open(path):
        toks = line.split()
        if len(toks) < 4: continue
        try:
            p = float(toks[1])/1e4; d = float(toks[3])*10.0
        except ValueError:
            continue
        if p > PMAX: continue
        if p > pc + 0.002 and d > 0.1:
            P.append(p); D.append(d)
    return np.asarray(P), np.asarray(D)


def sextic_d2(p, C, eta, pc):
    arg = np.maximum(1.0 + eta*(p/pc - 1.0), 0.0)
    return C*(np.sqrt(arg) - 1.0)


def main():
    here = os.path.dirname(os.path.abspath(__file__))
    out_dir = os.path.normpath(os.path.join(here, '..', 'data'))
    raw_root = os.path.normpath(os.path.join(here, RAW))
    os.makedirs(out_dir, exist_ok=True)

    pars_lines = ["# Direct sextic-Landau fit to Delta(p), broken phase.",
                  "# Columns: G_k_label  N_pts  pc[GPa]  C[A^2]+/-err  eta+/-err  U_eff/A[A^-2]  W/A[A^-4]  R^2"]
    fit_curves = {}
    for M in MEMBERS:
        path = os.path.join(raw_root, M['file'], 'data_pressure_delta.txt')
        p, d = read_broken_phase(path, M['pc'])
        if len(p) < 5:
            print(f"{M['label']}: too few points"); continue
        d2 = d**2
        try:
            popt, pcov = curve_fit(lambda x, C, eta: sextic_d2(x, C, eta, M['pc']),
                                   p, d2, p0=[1.0, 5.0],
                                   bounds=([1e-4, 1e-4], [1e3, 1e3]),
                                   maxfev=40000)
            C, eta = popt
            Cerr, eta_err = np.sqrt(np.diag(pcov))
            ua = 2.0/(C*eta); wa = 1.0/(C**2 * eta)
            ss_res = np.sum((d2 - sextic_d2(p, C, eta, M['pc']))**2)
            ss_tot = np.sum((d2 - d2.mean())**2)
            r2 = 1.0 - ss_res/ss_tot
            pars_lines.append(f"{M['label']:6} {len(p):4d}  {M['pc']:.4f}  "
                              f"{C:8.3f} {Cerr:8.3f}  {eta:9.4f} {eta_err:9.4f}  "
                              f"{ua:8.4f}  {wa:9.5f}  {r2:.4f}")
            print(pars_lines[-1])
            # Save fit curve (200 pts)
            pp = np.linspace(M['pc']+0.005, PMAX, 200)
            fit_curves[M['label']] = (pp, np.sqrt(np.maximum(sextic_d2(pp, C, eta, M['pc']), 0.0)),
                                      p, d)
        except Exception as e:
            print(f"{M['label']}: fit failed: {e}")
    with open(os.path.join(out_dir, 'figS_sextic_direct_pars.dat'), 'w') as f:
        f.write('\n'.join(pars_lines) + '\n')
    # Write per-member (p, Delta) data + fit on the same grid
    for label, (pp, dd, p_obs, d_obs) in fit_curves.items():
        outp = os.path.join(out_dir, f'figS_sextic_direct_{label.replace("_","")}_data.dat')
        with open(outp, 'w') as f:
            f.write(f"# {label} broken-phase observed (p[GPa], Delta[A])\n")
            for pi_, di_ in zip(p_obs, d_obs):
                f.write(f"{pi_:.5f}  {di_:.5f}\n")
        outf = os.path.join(out_dir, f'figS_sextic_direct_{label.replace("_","")}_fit.dat')
        with open(outf, 'w') as f:
            f.write(f"# {label} sextic fit curve (p[GPa], Delta[A])\n")
            for pi_, di_ in zip(pp, dd):
                f.write(f"{pi_:.5f}  {di_:.5f}\n")
    print(f"\nWrote {os.path.join(out_dir, 'figS_sextic_direct_pars.dat')}")

if __name__ == '__main__':
    main()
