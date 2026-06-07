#!/usr/bin/gnuplot
# Figure 3(a) of the manuscript:
# Volume shrinkage v = V/V_0 vs hydrostatic pressure for G_1 to G_5, with
# the piecewise-linear equation of state of Equation \ref{eq:kappa}:
#     v(p) = 1 + kappa_1 * p             p < p_c   (centric branch, intercept = 1)
#     v(p) = v_c^+ + kappa_2 * (p - p_c) p > p_c   (acentric branch, free intercept)
#
# The volume jump at the transition is Delta v = v_c^+ - (1 + kappa_1 * p_c).
#
# The fit is performed end-to-end inside this script (no external Python
# pre-processing). The data filter is omega1 >= 0, mechanically stable rows
# of the raw GULP/PLUMED pressure scan. The critical pressures p_c are held
# fixed at the values obtained from fig3b_delta.gp (Heaviside-power free fit
# of Delta vs pressure on the same data). A gap of half-width h = 0.06 GPa
# is excluded around each p_c during the linear fits.
#
# Inputs (raw pressure scans, NOT preprocessed):
#   ../../initial_structures_RHO_isoreticular/dir_RHO_isoreticular_G<k>_*/data_pressure_delta.txt

set term postscript eps color enhanced blacktext 'Helvetica,26'
set output 'fig3a_volume.eps'
set encoding iso_8859_1
set locale "en_GB.UTF-8"

set fit quiet
set fit errorvariables

# ---- Path to the raw scans (resolved relative to manuscript_figures/build/)
RAW = '../../initial_structures_RHO_isoreticular'

# ---- p_c per member: HELD FIXED at the soft-mode values from the
#      Cowley-Levanyuk extrapolation omega^2(p) = alpha (p_c - p) on the
#      parent-symmetry branch (see fig3e_softmode.gp). Identical to the
#      values listed in Table 1 of the manuscript.
pc_1 = 0.9418 ; pc_2 = 0.5338 ; pc_3 = 0.3627 ; pc_4 = 0.2241 ; pc_5 = 0.1282
h    = 0.06   # GPa gap excluded around each p_c during the linear fits

# Stability filter: short-format rows (NF<9) are treated as mechanically
# stable (they come from converged structural relaxations); rows with the
# full 21-column format use the sign of col 9 directly.
STABLE_FILT   = "(NF<9 || $9>=0)"
UNSTABLE_FILT = "(NF>=9 && $9<0)"

# ---- Read V_0 per member (first row, col 5) via stats --------------------
stats RAW."/dir_RHO_isoreticular_G1_222_SG_P1/data_pressure_delta.txt" every ::0::0 u 5 nooutput
V1 = STATS_min
stats RAW."/dir_RHO_isoreticular_G2_SG_P1/data_pressure_delta.txt" every ::0::0 u 5 nooutput
V2 = STATS_min
stats RAW."/dir_RHO_isoreticular_G3_SG_P1/data_pressure_delta.txt" every ::0::0 u 5 nooutput
V3 = STATS_min
stats RAW."/dir_RHO_isoreticular_G4_SG_P1/data_pressure_delta.txt" every ::0::0 u 5 nooutput
V4 = STATS_min
stats RAW."/dir_RHO_isoreticular_G5_SG_P1/data_pressure_delta.txt" every ::0::0 u 5 nooutput
V5 = STATS_min

# ---- Piecewise-linear models ---------------------------------------------
v1_lo(p) = 1.0 + k1_1*p
v1_hi(p) = b2_1 + k2_1*(p - pc_1)
v2_lo(p) = 1.0 + k1_2*p
v2_hi(p) = b2_2 + k2_2*(p - pc_2)
v3_lo(p) = 1.0 + k1_3*p
v3_hi(p) = b2_3 + k2_3*(p - pc_3)
v4_lo(p) = 1.0 + k1_4*p
v4_hi(p) = b2_4 + k2_4*(p - pc_4)
v5_lo(p) = 1.0 + k1_5*p
v5_hi(p) = b2_5 + k2_5*(p - pc_5)

# Initial guesses
k1_1=-0.014; k2_1=-0.14; b2_1=0.98
k1_2=-0.014; k2_2=-0.08; b2_2=0.99
k1_3=-0.014; k2_3=-0.08; b2_3=0.99
k1_4=-0.014; k2_4=-0.08; b2_4=0.99
k1_5=-0.014; k2_5=-0.07; b2_5=1.00

# ---- Fits (omega1 >= 0 filter, gap of width 2h excluded around p_c) -----
fit [0:pc_1-h] v1_lo(x) "< awk '{if".STABLE_FILT." print $2,$5}' ".RAW."/dir_RHO_isoreticular_G1_222_SG_P1/data_pressure_delta.txt" u ($1/1e4):($2/V1) via k1_1
fit [pc_1+h:]  v1_hi(x) "< awk '{if".STABLE_FILT." print $2,$5}' ".RAW."/dir_RHO_isoreticular_G1_222_SG_P1/data_pressure_delta.txt" u ($1/1e4):($2/V1) via k2_1, b2_1

fit [0:pc_2-h] v2_lo(x) "< awk '{if".STABLE_FILT." print $2,$5}' ".RAW."/dir_RHO_isoreticular_G2_SG_P1/data_pressure_delta.txt" u ($1/1e4):($2/V2) via k1_2
fit [pc_2+h:]  v2_hi(x) "< awk '{if".STABLE_FILT." print $2,$5}' ".RAW."/dir_RHO_isoreticular_G2_SG_P1/data_pressure_delta.txt" u ($1/1e4):($2/V2) via k2_2, b2_2

fit [0:pc_3-h] v3_lo(x) "< awk '{if".STABLE_FILT." print $2,$5}' ".RAW."/dir_RHO_isoreticular_G3_SG_P1/data_pressure_delta.txt" u ($1/1e4):($2/V3) via k1_3
fit [pc_3+h:]  v3_hi(x) "< awk '{if".STABLE_FILT." print $2,$5}' ".RAW."/dir_RHO_isoreticular_G3_SG_P1/data_pressure_delta.txt" u ($1/1e4):($2/V3) via k2_3, b2_3

fit [0:pc_4-h] v4_lo(x) "< awk '{if".STABLE_FILT." print $2,$5}' ".RAW."/dir_RHO_isoreticular_G4_SG_P1/data_pressure_delta.txt" u ($1/1e4):($2/V4) via k1_4
fit [pc_4+h:]  v4_hi(x) "< awk '{if".STABLE_FILT." print $2,$5}' ".RAW."/dir_RHO_isoreticular_G4_SG_P1/data_pressure_delta.txt" u ($1/1e4):($2/V4) via k2_4, b2_4

fit [0:pc_5-h] v5_lo(x) RAW."/dir_RHO_isoreticular_G5_SG_P1/data_pressure_delta.txt" u ($2/1e4):($5/V5) via k1_5
fit [pc_5+h:]  v5_hi(x) RAW."/dir_RHO_isoreticular_G5_SG_P1/data_pressure_delta.txt" u ($2/1e4):($5/V5) via k2_5, b2_5

# ---- Diagnostic print ----------------------------------------------------
print "============================================================"
print "fig3a: piecewise-linear v(p) fits (omega1>=0; pc from fig3b)"
print "Delta_v = v_c^+ - (1 + k1*pc)"
print sprintf("G_1: pc=%.4f k1=%.5f k2=%.4f b2=%.5f Dv=%.5f", pc_1, k1_1, k2_1, b2_1, b2_1 - (1.0 + k1_1*pc_1))
print sprintf("G_2: pc=%.4f k1=%.5f k2=%.4f b2=%.5f Dv=%.5f", pc_2, k1_2, k2_2, b2_2, b2_2 - (1.0 + k1_2*pc_2))
print sprintf("G_3: pc=%.4f k1=%.5f k2=%.4f b2=%.5f Dv=%.5f", pc_3, k1_3, k2_3, b2_3, b2_3 - (1.0 + k1_3*pc_3))
print sprintf("G_4: pc=%.4f k1=%.5f k2=%.4f b2=%.5f Dv=%.5f", pc_4, k1_4, k2_4, b2_4, b2_4 - (1.0 + k1_4*pc_4))
print sprintf("G_5: pc=%.4f k1=%.5f k2=%.4f b2=%.5f Dv=%.5f", pc_5, k1_5, k2_5, b2_5, b2_5 - (1.0 + k1_5*pc_5))
print "============================================================"

# ---- Plot ----------------------------------------------------------------
set xlabel "Hydrostatic pressure, {/Helvetica-Italic p} [GPa]"
set ylabel "{/Helvetica-Italic V}/{/Helvetica-Italic V}_0 [-]"
set xrange [0 : 2]
set yrange [0.80 : 1.01]
set xtics 0.5
set ytics 0.05
unset grid
set key bottom left reverse Left font ",22"

# Vertical dotted markers at each p_c
set arrow 1 from pc_1, 0.80 to pc_1, 1.01 nohead lc rgb '#cb4335' dt 2 lw 0.8
set arrow 2 from pc_2, 0.80 to pc_2, 1.01 nohead lc rgb '#f1c40f' dt 2 lw 0.8
set arrow 3 from pc_3, 0.80 to pc_3, 1.01 nohead lc rgb '#27ae60' dt 2 lw 0.8
set arrow 4 from pc_4, 0.80 to pc_4, 1.01 nohead lc rgb '#2874a6' dt 2 lw 0.8
set arrow 5 from pc_5, 0.80 to pc_5, 1.01 nohead lc rgb '#7d3c98' dt 2 lw 0.8

set samples 1001

plot \
    '+' u 1:($1 < pc_1 ? v1_lo($1) : v1_hi($1)) w l lw 2.0 lc rgb '#cb4335' notitle,\
    '+' u 1:($1 < pc_2 ? v2_lo($1) : v2_hi($1)) w l lw 2.0 lc rgb '#f1c40f' notitle,\
    '+' u 1:($1 < pc_3 ? v3_lo($1) : v3_hi($1)) w l lw 2.0 lc rgb '#27ae60' notitle,\
    '+' u 1:($1 < pc_4 ? v4_lo($1) : v4_hi($1)) w l lw 2.0 lc rgb '#2874a6' notitle,\
    '+' u 1:($1 < pc_5 ? v5_lo($1) : v5_hi($1)) w l lw 2.0 lc rgb '#7d3c98' notitle,\
    "< awk '{if".STABLE_FILT." print $2,$5}' ".RAW."/dir_RHO_isoreticular_G1_222_SG_P1/data_pressure_delta.txt" u ($1/1e4):($2/V1) w p pt 7 ps 1.0 lc rgb '#cb4335' title '{/Helvetica-Italic G}_1',\
    "< awk '{if".UNSTABLE_FILT." print $2,$5}' ".RAW."/dir_RHO_isoreticular_G1_222_SG_P1/data_pressure_delta.txt" u ($1/1e4):($2/V1) w p pt 6 ps 1.4 lc rgb '#cb4335' notitle,\
    "< awk '{if".STABLE_FILT." print $2,$5}' ".RAW."/dir_RHO_isoreticular_G2_SG_P1/data_pressure_delta.txt" u ($1/1e4):($2/V2) w p pt 7 ps 1.0 lc rgb '#f1c40f' title '{/Helvetica-Italic G}_2',\
    "< awk '{if".UNSTABLE_FILT." print $2,$5}' ".RAW."/dir_RHO_isoreticular_G2_SG_P1/data_pressure_delta.txt" u ($1/1e4):($2/V2) w p pt 6 ps 1.4 lc rgb '#f1c40f' notitle,\
    "< awk '{if".STABLE_FILT." print $2,$5}' ".RAW."/dir_RHO_isoreticular_G3_SG_P1/data_pressure_delta.txt" u ($1/1e4):($2/V3) w p pt 7 ps 1.0 lc rgb '#27ae60' title '{/Helvetica-Italic G}_3',\
    "< awk '{if".UNSTABLE_FILT." print $2,$5}' ".RAW."/dir_RHO_isoreticular_G3_SG_P1/data_pressure_delta.txt" u ($1/1e4):($2/V3) w p pt 6 ps 1.4 lc rgb '#27ae60' notitle,\
    "< awk '{if".STABLE_FILT." print $2,$5}' ".RAW."/dir_RHO_isoreticular_G4_SG_P1/data_pressure_delta.txt" u ($1/1e4):($2/V4) w p pt 7 ps 1.0 lc rgb '#2874a6' title '{/Helvetica-Italic G}_4',\
    "< awk '{if".UNSTABLE_FILT." print $2,$5}' ".RAW."/dir_RHO_isoreticular_G4_SG_P1/data_pressure_delta.txt" u ($1/1e4):($2/V4) w p pt 6 ps 1.4 lc rgb '#2874a6' notitle,\
    "< awk '{if".STABLE_FILT." print $2,$5}' ".RAW."/dir_RHO_isoreticular_G5_SG_P1/data_pressure_delta.txt" u ($1/1e4):($2/V5) w p pt 7 ps 1.0 lc rgb '#7d3c98' title '{/Helvetica-Italic G}_5',\
    "< awk '{if".UNSTABLE_FILT." print $2,$5}' ".RAW."/dir_RHO_isoreticular_G5_SG_P1/data_pressure_delta.txt" u ($1/1e4):($2/V5) w p pt 6 ps 1.4 lc rgb '#7d3c98' notitle
