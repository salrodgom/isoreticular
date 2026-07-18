#!/usr/bin/gnuplot
# Figure 3(e) of the manuscript (new subpanel of fig:Phonons):
# Squared soft-mode frequency omega^2 vs hydrostatic pressure for G_1 to G_5,
# with the Cowley-Levanyuk linear extrapolation
#     omega^2(p) = alpha (p_c - p)
# fitted independently for each member on the parent-symmetry monotonic
# branch. The zero-crossing of each fit identifies the critical pressure
# p_c, which is the canonical value used throughout the manuscript and
# tabulated in Table 1.
#
# Sign convention: omega^2 is positive for real (mechanically stable)
# modes and negative for imaginary (unstable) ones. The convention follows
# the original gp_Figure phonon scripts of the data repository.
#
# Inputs (raw pressure scans, NOT preprocessed):
#   ../../initial_structures_RHO_isoreticular/dir_RHO_isoreticular_G<k>_*/data_pressure_delta.txt
#
# The awk filter selects rows on the parent-symmetry branch (D8R
# distortion D < 0.005 nm = 0.05 angstrom) with the lowest-mode frequency
# populated either as a real value (col 12 > 0.5 cm^-1) or as an
# imaginary value (col 9 < -0.5 cm^-1). The signed omega^2 is then
# +(col 12)^2 or -(col 9)^2, respectively.

set term postscript eps color enhanced blacktext 'Helvetica,26'
set output 'fig3e_softmode.eps'
set encoding iso_8859_1
set locale "en_GB.UTF-8"

set fit quiet
set fit errorvariables

RAW = '../../initial_structures_RHO_isoreticular'

# Awk filters (parameterised by PMAX per member): extract (p [GPa], omega^2
# [cm^-2]) on the parent-symmetry branch up to a hard pressure cutoff.
# Two complementary streams, mirroring the convention used in fig3a_volume.gp
# and fig3b_delta.gp:
#   REAL  (filled circles, pt 7): RFO Hessian converged with all real
#         frequencies. col 12 > 0.5 cm^-1, parent-symmetry (col 4 < 0.015 nm)
#         and col 9 >= 0. omega^2 = +col12^2.
#   IMAG  (open circles,  pt 6): RFO Hessian carries a negative-curvature
#         direction; the structure was relaxed instead with the GULP
#         lower-symmetry option, which yields a parent-branch imaginary
#         soft mode reported in col 9 < -0.5 cm^-1. omega^2 = -col9^2.
# The PMAX cutoff drops post-transition metastable-centric reentries that
# would otherwise pollute both the fit and the plot (see Section S
# G5-phonon of the SI). The fit uses both streams jointly.
W2_REAL(pmax) = sprintf("awk -v pmax=%g 'NF>=12 { p=$2/1e4; if(p>pmax) next; if($12>0.5 && $4<0.015 && $9>=0){print p, $12*$12} }'", pmax)
W2_IMAG(pmax) = sprintf("awk -v pmax=%g 'NF>=12 { p=$2/1e4; if(p>pmax) next; if($9<-0.5){print p, -($9*$9)} }'", pmax)
W2(pmax)      = sprintf("awk -v pmax=%g 'NF>=12 { p=$2/1e4; if(p>pmax) next; if($9<-0.5){w2=-($9*$9)} else if($12>0.5 && $4<0.015 && $9>=0){w2=$12*$12} else next; print p, w2 }'", pmax)

# Pressure cutoffs per member (also used as the fit window). Chosen to
# include all monotonically-decreasing omega^2 points and to drop any
# spurious reappearance of a positive real mode above p_c.
PMAX1 = 1.10
PMAX2 = 0.75
PMAX3 = 0.50
PMAX4 = 0.25
PMAX5 = 0.145

# Cowley-Levanyuk linear law per member
f(x, alpha, pc) = alpha * (pc - x)

a1=200; pc_1=0.94
a2=200; pc_2=0.53
a3=200; pc_3=0.34
a4=200; pc_4=0.23
a5=200; pc_5=0.13

fit f(x, a1, pc_1) "< ".W2(PMAX1)." ".RAW."/dir_RHO_isoreticular_G1_222_SG_P1/data_pressure_delta.txt" u 1:2 via a1, pc_1
fit f(x, a2, pc_2) "< ".W2(PMAX2)." ".RAW."/dir_RHO_isoreticular_G2_SG_P1/data_pressure_delta.txt"      u 1:2 via a2, pc_2
fit f(x, a3, pc_3) "< ".W2(PMAX3)." ".RAW."/dir_RHO_isoreticular_G3_SG_P1/data_pressure_delta.txt"      u 1:2 via a3, pc_3
fit f(x, a4, pc_4) "< ".W2(PMAX4)." ".RAW."/dir_RHO_isoreticular_G4_SG_P1/data_pressure_delta.txt"      u 1:2 via a4, pc_4
fit f(x, a5, pc_5) "< ".W2(PMAX5)." ".RAW."/dir_RHO_isoreticular_G5_SG_P1/data_pressure_delta.txt"      u 1:2 via a5, pc_5

print "============================================================"
print "fig3e: Cowley-Levanyuk omega^2(p) = alpha (pc - p)  fits"
print "============================================================"
print sprintf("G_1: pc=%.4f +/- %.4f GPa  alpha=%.2f cm^-2/GPa", pc_1, pc_1_err, a1)
print sprintf("G_2: pc=%.4f +/- %.4f GPa  alpha=%.2f cm^-2/GPa", pc_2, pc_2_err, a2)
print sprintf("G_3: pc=%.4f +/- %.4f GPa  alpha=%.2f cm^-2/GPa", pc_3, pc_3_err, a3)
print sprintf("G_4: pc=%.4f +/- %.4f GPa  alpha=%.2f cm^-2/GPa", pc_4, pc_4_err, a4)
print sprintf("G_5: pc=%.4f +/- %.4f GPa  alpha=%.2f cm^-2/GPa", pc_5, pc_5_err, a5)

set xlabel "Hydrostatic pressure, {/Helvetica-Italic p} [GPa]"
set ylabel "{/Symbol w}_1^{ 2} [cm^{-2}]"
set xrange [-0.05 : 1.20]
set yrange [-200 : 800]
set xtics 0.2
set ytics 200
unset grid
set key top right reverse Left font ",22"

# Horizontal reference at omega^2 = 0 (the Cowley-Levanyuk zero-crossing)
set arrow 100 from -0.05,0 to 1.20,0 nohead lc rgb '#888888' dt 2 lw 0.6

# Vertical dashed line at each fitted p_c (same convention as fig3a/fig3b)
set arrow 101 from pc_1,-200 to pc_1,800 nohead lc rgb '#cb4335' dt 2 lw 0.8
set arrow 102 from pc_2,-200 to pc_2,800 nohead lc rgb '#f1c40f' dt 2 lw 0.8
set arrow 103 from pc_3,-200 to pc_3,800 nohead lc rgb '#27ae60' dt 2 lw 0.8
set arrow 104 from pc_4,-200 to pc_4,800 nohead lc rgb '#2874a6' dt 2 lw 0.8
set arrow 105 from pc_5,-200 to pc_5,800 nohead lc rgb '#7d3c98' dt 2 lw 0.8

set samples 1001

plot \
    '+' u 1:(f($1, a1, pc_1)) w l lw 2.0 lc rgb '#cb4335' notitle,\
    '+' u 1:(f($1, a2, pc_2)) w l lw 2.0 lc rgb '#f1c40f' notitle,\
    '+' u 1:(f($1, a3, pc_3)) w l lw 2.0 lc rgb '#27ae60' notitle,\
    '+' u 1:(f($1, a4, pc_4)) w l lw 2.0 lc rgb '#2874a6' notitle,\
    '+' u 1:(f($1, a5, pc_5)) w l lw 2.0 lc rgb '#7d3c98' notitle,\
    "< ".W2_REAL(PMAX1)." ".RAW."/dir_RHO_isoreticular_G1_222_SG_P1/data_pressure_delta.txt" u 1:2 w p pt 7 ps 1.2 lc rgb '#cb4335' title '{/Helvetica-Italic G}_1',\
    "< ".W2_IMAG(PMAX1)." ".RAW."/dir_RHO_isoreticular_G1_222_SG_P1/data_pressure_delta.txt" u 1:2 w p pt 6 ps 1.4 lc rgb '#cb4335' notitle,\
    "< ".W2_REAL(PMAX2)." ".RAW."/dir_RHO_isoreticular_G2_SG_P1/data_pressure_delta.txt"      u 1:2 w p pt 7 ps 1.2 lc rgb '#f1c40f' title '{/Helvetica-Italic G}_2',\
    "< ".W2_IMAG(PMAX2)." ".RAW."/dir_RHO_isoreticular_G2_SG_P1/data_pressure_delta.txt"      u 1:2 w p pt 6 ps 1.4 lc rgb '#f1c40f' notitle,\
    "< ".W2_REAL(PMAX3)." ".RAW."/dir_RHO_isoreticular_G3_SG_P1/data_pressure_delta.txt"      u 1:2 w p pt 7 ps 1.4 lc rgb '#27ae60' title '{/Helvetica-Italic G}_3',\
    "< ".W2_IMAG(PMAX3)." ".RAW."/dir_RHO_isoreticular_G3_SG_P1/data_pressure_delta.txt"      u 1:2 w p pt 6 ps 1.4 lc rgb '#27ae60' notitle,\
    "< ".W2_REAL(PMAX4)." ".RAW."/dir_RHO_isoreticular_G4_SG_P1/data_pressure_delta.txt"      u 1:2 w p pt 7 ps 1.4 lc rgb '#2874a6' title '{/Helvetica-Italic G}_4',\
    "< ".W2_IMAG(PMAX4)." ".RAW."/dir_RHO_isoreticular_G4_SG_P1/data_pressure_delta.txt"      u 1:2 w p pt 6 ps 1.4 lc rgb '#2874a6' notitle,\
    "< ".W2_REAL(PMAX5)." ".RAW."/dir_RHO_isoreticular_G5_SG_P1/data_pressure_delta.txt"      u 1:2 w p pt 7 ps 1.4 lc rgb '#7d3c98' title '{/Helvetica-Italic G}_5',\
    "< ".W2_IMAG(PMAX5)." ".RAW."/dir_RHO_isoreticular_G5_SG_P1/data_pressure_delta.txt"      u 1:2 w p pt 6 ps 1.4 lc rgb '#7d3c98' notitle
