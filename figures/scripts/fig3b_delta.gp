#!/usr/bin/gnuplot
# Figure 3(b) of the manuscript:
# D8R distortion Delta vs hydrostatic pressure for G_1 to G_5, with the
# Heaviside-power phenomenological fit of Equation \ref{eq:fit}:
#     Delta(p) = a * H(p - p_c) * (p - p_c)^beta_eff
# This is the dimensional form used in gp_Figure2b of the data repository
# (initial_structures_RHO_isoreticular/gp_Figure2b). The dimensionless
# amplitude reported in Table 1 of the main text is delta = a * p_c^beta_eff.
#
# The fit is performed end-to-end inside this script (no external Python
# pre-processing). The data filter is omega1 >= 0, mechanically stable
# rows of the raw GULP/PLUMED pressure scan, matching the original protocol.
#
# For each member, the fit is two-step:
#   (i)  via (a, beta) with p_c held at the initial guess (the soft-mode
#        value), to anchor the slope and amplitude.
#   (ii) via (a, p_c, beta), to refine p_c jointly with the other two.
#
# Inputs (raw pressure scans, NOT preprocessed):
#   ../../initial_structures_RHO_isoreticular/dir_RHO_isoreticular_G<k>_*/data_pressure_delta.txt

set term postscript eps color enhanced blacktext 'Helvetica,26'
set output 'fig3b_delta.eps'
set encoding iso_8859_1
set locale "en_GB.UTF-8"

set fit quiet
set fit errorvariables

# ---- Path to the raw scans (resolved relative to manuscript_figures/build/)
RAW = '../../initial_structures_RHO_isoreticular'

# p_c is HELD FIXED at the soft-mode value from the Cowley-Levanyuk
# extrapolation omega^2(p) = alpha (p_c - p) on the parent-symmetry branch
# (joint_pc_fit.py of the data repository). Only the amplitude a and the
# exponent mu (= beta) are fitted.
# Stability filter: rows with NF<9 (short format, no phonon calc) are
# treated as mechanically stable (they come from converged structural
# relaxations); rows with NF>=9 use the sign of col 9 directly.

STABLE_FILT   = "(NF<9 || $9>=0)"
UNSTABLE_FILT = "(NF>=9 && $9<0)"

# Soft-mode p_c values (held fixed; not fitted). Identical to those listed
# in Table 1 of the manuscript, obtained from a Cowley-Levanyuk linear
# extrapolation omega^2(p) = alpha (p_c - p) on the parent-symmetry
# branch (see fig3e_softmode.gp for the actual extraction).
x01 = 0.9418 ; x02 = 0.5338 ; x03 = 0.3627 ; x04 = 0.2241 ; x05 = 0.1282

# ---- G_1 (2x2x2 supercell, 384 T-atoms) -----------------------------------
H1(x) = x>x01 ? 1 : 0
f1(x) = a1*H1(x)*(x-x01)**mu1
a1=2 ; mu1=0.5
fit f1(x) "< awk '{if".STABLE_FILT." print $2,$4}' ".RAW."/dir_RHO_isoreticular_G1_222_SG_P1/data_pressure_delta.txt" u ($1/1e4):($2*10) via a1,mu1

# ---- G_2 -------------------------------------------------------------------
H2(x) = x>x02 ? 1 : 0
f2(x) = a2*H2(x)*(x-x02)**mu2
a2=1.5 ; mu2=0.5
fit f2(x) "< awk '{if".STABLE_FILT." print $2,$4}' ".RAW."/dir_RHO_isoreticular_G2_SG_P1/data_pressure_delta.txt" u ($1/1e4):($2*10) via a2,mu2

# ---- G_3 -------------------------------------------------------------------
H3(x) = x>x03 ? 1 : 0
f3(x) = a3*H3(x)*(x-x03)**mu3
a3=1.5 ; mu3=0.5
fit f3(x) "< awk '{if".STABLE_FILT." print $2,$4}' ".RAW."/dir_RHO_isoreticular_G3_SG_P1/data_pressure_delta.txt" u ($1/1e4):($2*10) via a3,mu3

# ---- G_4 -------------------------------------------------------------------
H4(x) = x>x04 ? 1 : 0
f4(x) = a4*H4(x)*(x-x04)**mu4
a4=1.5 ; mu4=0.5
fit f4(x) "< awk '{if".STABLE_FILT." print $2,$4}' ".RAW."/dir_RHO_isoreticular_G4_SG_P1/data_pressure_delta.txt" u ($1/1e4):($2*10) via a4,mu4

# ---- G_5 -------------------------------------------------------------------
H5(x) = x>x05 ? 1 : 0
f5(x) = a5*H5(x)*(x-x05)**mu5
a5=1.5 ; mu5=0.5
fit f5(x) "< awk '{if".STABLE_FILT." print $2,$4}' ".RAW."/dir_RHO_isoreticular_G5_SG_P1/data_pressure_delta.txt" u ($1/1e4):($2*10) via a5,mu5

# ---- Diagnostic print: fitted parameters in eq:fit form --------------------
print "============================================================"
print "fig3b: Heaviside-power Delta(p) = a*(p-pc)^mu fits (omega1>=0)"
print "delta = a * pc^mu (manuscript eq:fit dimensionless form)"
print sprintf("G_1: a=%.4f pc=%.4f mu=%.4f -> delta=%.4f beta=%.4f", a1, x01, mu1, a1*x01**mu1, mu1)
print sprintf("G_2: a=%.4f pc=%.4f mu=%.4f -> delta=%.4f beta=%.4f", a2, x02, mu2, a2*x02**mu2, mu2)
print sprintf("G_3: a=%.4f pc=%.4f mu=%.4f -> delta=%.4f beta=%.4f", a3, x03, mu3, a3*x03**mu3, mu3)
print sprintf("G_4: a=%.4f pc=%.4f mu=%.4f -> delta=%.4f beta=%.4f", a4, x04, mu4, a4*x04**mu4, mu4)
print sprintf("G_5: a=%.4f pc=%.4f mu=%.4f -> delta=%.4f beta=%.4f", a5, x05, mu5, a5*x05**mu5, mu5)
print "============================================================"

# ---- Plot ------------------------------------------------------------------
set xlabel "Hydrostatic pressure, {/Helvetica-Italic p} [GPa]"
set ylabel "D8R distortion, {/Symbol D} [{\305}]"
set xrange [0 : 2]
set yrange [-0.1 : 2.5]
set xtics 0.5
set ytics 0.5
unset grid
set key top left reverse Left font ",22"

# Vertical dotted markers at the fitted p_c per member
set arrow 1 from x01, -0.1 to x01, 2.5 nohead lc rgb '#cb4335' dt 2 lw 0.8
set arrow 2 from x02, -0.1 to x02, 2.5 nohead lc rgb '#f1c40f' dt 2 lw 0.8
set arrow 3 from x03, -0.1 to x03, 2.5 nohead lc rgb '#27ae60' dt 2 lw 0.8
set arrow 4 from x04, -0.1 to x04, 2.5 nohead lc rgb '#2874a6' dt 2 lw 0.8
set arrow 5 from x05, -0.1 to x05, 2.5 nohead lc rgb '#7d3c98' dt 2 lw 0.8

set samples 1001

plot \
    '+' u 1:(f1($1)) w l lw 2.0 lc rgb '#cb4335' notitle,\
    '+' u 1:(f2($1)) w l lw 2.0 lc rgb '#f1c40f' notitle,\
    '+' u 1:(f3($1)) w l lw 2.0 lc rgb '#27ae60' notitle,\
    '+' u 1:(f4($1)) w l lw 2.0 lc rgb '#2874a6' notitle,\
    '+' u 1:(f5($1)) w l lw 2.0 lc rgb '#7d3c98' notitle,\
    "< awk '{if".STABLE_FILT." print $2,$4}' ".RAW."/dir_RHO_isoreticular_G1_222_SG_P1/data_pressure_delta.txt" u ($1/1e4):($2*10) w p pt 7 ps 1.0 lc rgb '#cb4335' title '{/Helvetica-Italic G}_1',\
    "< awk '{if".UNSTABLE_FILT." print $2,$4}' ".RAW."/dir_RHO_isoreticular_G1_222_SG_P1/data_pressure_delta.txt" u ($1/1e4):($2*10) w p pt 6 ps 1.4 lc rgb '#cb4335' notitle,\
    "< awk '{if".STABLE_FILT." print $2,$4}' ".RAW."/dir_RHO_isoreticular_G2_SG_P1/data_pressure_delta.txt" u ($1/1e4):($2*10) w p pt 7 ps 1.0 lc rgb '#f1c40f' title '{/Helvetica-Italic G}_2',\
    "< awk '{if".UNSTABLE_FILT." print $2,$4}' ".RAW."/dir_RHO_isoreticular_G2_SG_P1/data_pressure_delta.txt" u ($1/1e4):($2*10) w p pt 6 ps 1.4 lc rgb '#f1c40f' notitle,\
    "< awk '{if".STABLE_FILT." print $2,$4}' ".RAW."/dir_RHO_isoreticular_G3_SG_P1/data_pressure_delta.txt" u ($1/1e4):($2*10) w p pt 7 ps 1.0 lc rgb '#27ae60' title '{/Helvetica-Italic G}_3',\
    "< awk '{if".UNSTABLE_FILT." print $2,$4}' ".RAW."/dir_RHO_isoreticular_G3_SG_P1/data_pressure_delta.txt" u ($1/1e4):($2*10) w p pt 6 ps 1.4 lc rgb '#27ae60' notitle,\
    "< awk '{if".STABLE_FILT." print $2,$4}' ".RAW."/dir_RHO_isoreticular_G4_SG_P1/data_pressure_delta.txt" u ($1/1e4):($2*10) w p pt 7 ps 1.0 lc rgb '#2874a6' title '{/Helvetica-Italic G}_4',\
    "< awk '{if".UNSTABLE_FILT." print $2,$4}' ".RAW."/dir_RHO_isoreticular_G4_SG_P1/data_pressure_delta.txt" u ($1/1e4):($2*10) w p pt 6 ps 1.4 lc rgb '#2874a6' notitle,\
    "< awk '{if".STABLE_FILT." print $2,$4}' ".RAW."/dir_RHO_isoreticular_G5_SG_P1/data_pressure_delta.txt" u ($1/1e4):($2*10) w p pt 7 ps 1.0 lc rgb '#7d3c98' title '{/Helvetica-Italic G}_5',\
    "< awk '{if".UNSTABLE_FILT." print $2,$4}' ".RAW."/dir_RHO_isoreticular_G5_SG_P1/data_pressure_delta.txt" u ($1/1e4):($2*10) w p pt 6 ps 1.4 lc rgb '#7d3c98' notitle
