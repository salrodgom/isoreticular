#!/usr/bin/env python3
"""
Two-step self-consistent fit of the critical pressure p_c and the
phenomenological parameters per isoreticular member G_k.

Step 1 (canonical p_c from soft mode):
    omega^2(p) = alpha (p_c - p)
    Linear regression over all rows where the parent-symmetry soft mode is
    known (real positive in col 12, or imaginary with magnitude in col 9).
    p_c = -intercept/slope. Single dataset, two parameters (alpha, p_c).

Step 2a (Delta with fixed p_c):
    Delta(p) = delta (p/p_c - 1)^beta  for p > p_c
    Fit on the early-growth window (Delta < cap A) to suppress the
    high-pressure saturation bias of the Heaviside-power form.

Step 2b (kappa_1, kappa_2 with fixed p_c):
    V(p)/V_0 = 1 + kappa_1 p              p < p_c    (linear, intercept 1)
    V(p)/V_0 = b2 + kappa_2 (p - p_c)     p > p_c    (linear, free intercept)
    Two independent linear fits, no continuity at p_c. Column 5 of
    data_pressure_delta.txt is the volume in nm^3 (from PLUMED's VOLUME
    action, via extract.sh), so V/V_0 is the genuine bulk-volume ratio
    and kappa here is the volumetric compressibility -dV/Vdp.

Output: table comparing the new self-consistent values per G_k with the
values currently in tab:parameters and tab:SI-volume of the manuscript.
"""

import numpy as np
import os

PAPER = {
    1: dict(pc=0.940, delta=2.213, beta=0.358, k1=-0.01344, k2=-0.1404),
    2: dict(pc=0.525, delta=1.221, beta=0.452, k1=-0.01340, k2=-0.0810),
    3: dict(pc=0.353, delta=0.964, beta=0.460, k1=-0.01371, k2=-0.0745),
    4: dict(pc=0.199, delta=0.657, beta=0.511, k1=-0.01521, k2=-0.0799),
    5: dict(pc=0.160, delta=0.589, beta=0.504, k1=-0.01694, k2=-0.0570),
}

BASE = os.path.dirname(os.path.abspath(__file__))
FILES = {
    1: os.path.join(BASE, "dir_RHO_isoreticular_G1_222_SG_P1",
                    "data_pressure_delta.txt"),
    2: os.path.join(BASE, "dir_RHO_isoreticular_G2_SG_P1",
                    "data_pressure_delta.txt"),
    3: os.path.join(BASE, "dir_RHO_isoreticular_G3_SG_P1",
                    "data_pressure_delta.txt"),
    4: os.path.join(BASE, "dir_RHO_isoreticular_G4_SG_P1",
                    "data_pressure_delta.txt"),
    5: os.path.join(BASE, "dir_RHO_isoreticular_G5_SG_P1",
                    "data_pressure_delta.txt"),
}

DELTA_FIT_CAP_ANGSTROM = 1.0   # Δ < 1 Å for the Heaviside-power fit window


def load(path, delta_threshold_nm=0.005):
    """Return p_GPa, Delta_A, V_nm3, omega_signed_sq numpy arrays.

    omega_signed_sq is restricted to the PARENT-SYMMETRY branch only: rows
    where Delta_nm < delta_threshold_nm (cubic Im-3m, distorted or not).
    Post-transition points where the I-43m soft mode "reappears" (Delta
    finite, col 12 large) are not part of the linear Cowley-Levanyuk branch
    and are excluded so that they do not bias the omega^2(p) fit.
    """
    rows = []
    with open(path) as f:
        for line in f:
            parts = line.split()
            if len(parts) < 5:
                continue
            try:
                p_bar = float(parts[1])
                D_nm = float(parts[3])
                V_nm3 = float(parts[4])
            except (ValueError, IndexError):
                continue
            w_imag = w_real = np.nan
            if len(parts) >= 9:
                try:    w_imag = float(parts[8])   # col 9
                except: pass
            if len(parts) >= 12:
                try:    w_real = float(parts[11])  # col 12
                except: pass

            # Restrict to parent-symmetry branch (Delta near zero).
            is_parent = D_nm < delta_threshold_nm
            w2 = np.nan
            if is_parent:
                if np.isfinite(w_real) and w_real > 0.5:
                    w2 = w_real**2
                elif np.isfinite(w_imag) and w_imag < -0.5:
                    w2 = -(w_imag**2)
            rows.append((p_bar/10000.0, D_nm*10.0, V_nm3, w2))
    arr = np.array(rows)
    return arr[:, 0], arr[:, 1], arr[:, 2], arr[:, 3]


def fit_pc_softmode(p, w2):
    """Linear least squares: omega^2 = alpha*(p_c - p).

    The Cowley-Levanyuk extrapolation is valid only on the parent-symmetry
    branch where omega^2 decreases monotonically with p (real, going to
    zero at p_c, then continues into negative-omega^2 imaginary territory
    just above p_c). We enforce monotonic decrease in p by sorting and
    truncating at the first reversal. This drops contaminating points
    where the soft mode "reappears" because the system has fallen into
    a different local minimum (a known artefact in GULP optimisations
    near the transition pressure).

    Returns (pc, alpha, sigma_pc, sigma_alpha, N).
    """
    ok = np.isfinite(w2)
    p_, w2_ = p[ok], w2[ok]
    # Sort by pressure, then enforce monotonic-decreasing omega^2
    order = np.argsort(p_)
    p_, w2_ = p_[order], w2_[order]
    if len(p_) >= 2:
        keep = [0]
        last_w2 = w2_[0]
        for i in range(1, len(p_)):
            if w2_[i] <= last_w2 + 1.0:  # allow 1 cm^-2 noise tolerance
                keep.append(i)
                last_w2 = w2_[i]
            else:
                break  # reversal: stop
        p_, w2_ = p_[keep], w2_[keep]
    if len(p_) < 2:
        return np.nan, np.nan, np.nan, np.nan, len(p_)
    # Model: w2 = b - a*p  (b = alpha*pc, a = alpha)
    A = np.column_stack([np.ones_like(p_), -p_])
    coef, residuals, rank, _ = np.linalg.lstsq(A, w2_, rcond=None)
    b, a_slope = coef
    pc = b / a_slope
    # Covariance estimate
    res = w2_ - A @ coef
    dof = max(1, len(p_) - 2)
    sigma2 = (res @ res) / dof
    try:
        cov = sigma2 * np.linalg.inv(A.T @ A)
        # Propagate sigma_pc from sigma_b, sigma_a via pc = b/a_slope
        sb, sa = np.sqrt(np.maximum(np.diag(cov), 0.0))
        sigma_pc = abs(pc) * np.sqrt((sb/b)**2 + (sa/a_slope)**2)
    except np.linalg.LinAlgError:
        sb = sa = sigma_pc = np.nan
    return pc, a_slope, sigma_pc, sa, len(p_)


def fit_delta_window(p, D, pc, cap=None):
    """Nonlinear fit Delta = delta*(p/pc - 1)**beta with pc fixed.
    If cap is given, restrict to D < cap; otherwise use all p > pc points.
    Returns (delta, beta, sigma_delta, sigma_beta, N).
    """
    from scipy.optimize import curve_fit
    ok = np.isfinite(p) & np.isfinite(D) & (p > pc) & (D > 0)
    if cap is not None:
        ok = ok & (D < cap)
    p_, D_ = p[ok], D[ok]
    if len(p_) < 3:
        return np.nan, np.nan, np.nan, np.nan, len(p_)
    f = lambda x, delta, beta: delta * (x/pc - 1.0)**beta
    try:
        popt, pcov = curve_fit(f, p_, D_, p0=[1.0, 0.5],
                               bounds=([0.05, 0.05], [10.0, 1.50]))
        sigmas = np.sqrt(np.maximum(np.diag(pcov), 0.0))
        return popt[0], popt[1], sigmas[0], sigmas[1], len(p_)
    except Exception:
        return np.nan, np.nan, np.nan, np.nan, len(p_)


def fit_delta_meanfield(p, D, pc, cap):
    """Heaviside-power fit with beta = 1/2 fixed (mean-field). Returns
    (delta_mf, sigma_delta_mf, rms_residual, N)."""
    from scipy.optimize import curve_fit
    ok = np.isfinite(p) & np.isfinite(D) & (p > pc) & (D > 0) & (D < cap)
    p_, D_ = p[ok], D[ok]
    if len(p_) < 2:
        return np.nan, np.nan, np.nan, len(p_)
    arg = p_/pc - 1.0
    f = lambda x, delta: delta * np.sqrt(np.maximum(x, 0.0))
    try:
        popt, pcov = curve_fit(f, p_, D_, p0=[1.0])
        delta_mf = popt[0]
        sigma_delta_mf = float(np.sqrt(max(0.0, pcov[0, 0])))
        residual = D_ - delta_mf * np.sqrt(np.maximum(arg, 0.0))
        rms = float(np.sqrt(np.mean(residual**2)))
        return delta_mf, sigma_delta_mf, rms, len(p_)
    except Exception:
        return np.nan, np.nan, np.nan, len(p_)


def fit_kappa(p, V_nm3, pc):
    """Two independent linear fits of v = V/V_0 vs p:
        v = 1 + k1 * p          p < pc
        v = b2 + k2 * (p - pc)  p > pc
    Returns (k1, k2, sigma_k1, sigma_k2, b2, N_below, N_above).
    Here V_nm3 is the volume in nm^3 from PLUMED (column 5 of
    data_pressure_delta.txt). kappa = dv/dp is the genuine volumetric
    compressibility -dV/Vdp.
    """
    ok = np.isfinite(p) & np.isfinite(V_nm3)
    p_, V_ = p[ok], V_nm3[ok]
    idx_p0 = np.argmin(np.abs(p_))
    V0 = V_[idx_p0]
    v = V_ / V0

    below = p_ < pc
    above = p_ > pc
    # Below pc: v - 1 = k1 * p (single slope, intercept fixed at 1)
    k1 = sk1 = np.nan
    if below.sum() >= 2:
        p_b, v_b = p_[below], v[below]
        A = p_b.reshape(-1, 1)
        sol, _, _, _ = np.linalg.lstsq(A, v_b - 1.0, rcond=None)
        k1 = sol[0]
        res = v_b - 1.0 - k1 * p_b
        if len(p_b) > 1:
            sk1 = np.sqrt((res @ res) / max(1, len(p_b) - 1) / (p_b @ p_b))

    # Above pc: v = b2 + k2 * (p - pc) (slope k2, intercept b2 free)
    k2 = sk2 = b2 = np.nan
    if above.sum() >= 2:
        p_a = p_[above] - pc
        v_a = v[above]
        A = np.column_stack([np.ones_like(p_a), p_a])
        sol, _, _, _ = np.linalg.lstsq(A, v_a, rcond=None)
        b2, k2 = sol
        res = v_a - b2 - k2 * p_a
        dof = max(1, len(p_a) - 2)
        sigma2 = (res @ res) / dof
        try:
            cov = sigma2 * np.linalg.inv(A.T @ A)
            sk2 = np.sqrt(max(0.0, cov[1, 1]))
        except np.linalg.LinAlgError:
            pass

    return k1, k2, sk1, sk2, b2, int(below.sum()), int(above.sum())


def main():
    print(f"\n{'='*30}  Self-consistent fit results with soft-mode p_c  {'='*30}")
    print("Step 1: p_c from omega^2(p) = alpha (p_c - p) (Cowley-Levanyuk)")
    print("Step 2 (with p_c fixed): delta, beta from Heaviside-power on all p>pc points")
    print("Step 3 (with p_c fixed): kappa_1 (linear v=1+k1*p below pc),")
    print("                          kappa_2 (linear v=b2+k2*(p-pc) above pc, independent)")
    print("Also reported: delta_MF with beta=1/2 fixed (mean-field) for comparison")
    print()
    print(f"{'k':>2} | "
          f"{'p_c (soft)':>15} | {'p_c paper':>9} | "
          f"{'delta (free β)':>16} | {'δ paper':>8} | "
          f"{'beta (free)':>14} | {'β paper':>8} | "
          f"{'δ_MF (β=½)':>13} | "
          f"{'kappa_1':>16} | {'k1 paper':>9} | "
          f"{'kappa_2':>16} | {'k2 paper':>9} | "
          f"{'Nph':>3} {'NΔ':>3} {'Nv<':>3} {'Nv>':>3}")
    print("-" * 220)

    rows_out = []
    for k in (1, 2, 3, 4, 5):
        p, D, V_nm3, w2 = load(FILES[k])

        pc, alpha, spc, sa, n_ph = fit_pc_softmode(p, w2)
        # Free β fit, full broken-phase window
        delta, beta, sd, sb, n_d = fit_delta_window(p, D, pc, cap=None)
        # Reference: β=1/2 fixed
        delta_mf, sd_mf, rms_mf, n_d_mf = fit_delta_meanfield(
            p, D, pc, DELTA_FIT_CAP_ANGSTROM)
        k1, k2, sk1, sk2, b2, nv1, nv2 = fit_kappa(p, V_nm3, pc)

        r = PAPER[k]
        print(f"{k:>2} | "
              f"{pc:>8.4f}+/-{spc:.4f} | {r['pc']:>9.4f} | "
              f"{delta:>9.3f}+/-{sd:.3f} | {r['delta']:>8.3f} | "
              f"{beta:>7.3f}+/-{sb:.3f} | {r['beta']:>8.3f} | "
              f"{delta_mf:>8.3f}+/-{sd_mf:.3f} | "
              f"{k1:>9.5f}+/-{sk1:.5f} | {r['k1']:>9.5f} | "
              f"{k2:>9.4f}+/-{sk2:.4f} | {r['k2']:>9.4f} | "
              f"{n_ph:>3d} {n_d:>3d} {nv1:>3d} {nv2:>3d}")
        rows_out.append((k, pc, spc, delta, sd, beta, sb, delta_mf, sd_mf,
                         k1, sk1, k2, sk2))
    return rows_out


if __name__ == "__main__":
    main()
