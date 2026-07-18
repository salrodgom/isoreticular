#!/usr/bin/gnuplot
# Generates fig5_isopressure.eps (panel (e) of Figure 5, fig:Mechanical).
# Plots the soft-mode-derived critical pressure p_c vs isoreticular order
# k-1 for G_1-G_5 on a log y-axis, together with the unanchored
# exponential fit p_c(k-1) = pi exp(-alpha (k-1)) extrapolated to G_6-G_8.
# The per-Si enthalpy relative to alpha-quartz is overlaid as blue bars
# on a secondary right y-axis.
#
# Input files (consumed via gnuplot 'using' on column $0 = row index):
#   ../data/fig5_isopressure_data.dat        soft-mode p_c for G_1-G_5
#                                             $1=label  $2=p_c  $3=err
#                                             $4=N_T    $5=E_per_Si[eV]
#   ../data/fig5_isopressure_predictions.dat extrapolated p_c for G_6-G_8
#                                             $1=label  $2=p_c  $3=N_T

set term postscript eps color enhanced blacktext 'Helvetica,26'
set output 'fig5_isopressure.eps'
set encoding iso_8859_1
set locale "en_GB.UTF-8"
unset grid
set style fill solid
set boxwidth 0.4
set fit brief errorvariables errorscaling prescale nowrap v5

DATA = '../data/fig5_isopressure_data.dat'
PRED = '../data/fig5_isopressure_predictions.dat'

# Fit p_c(k-1) = pi * exp(-alpha (k-1)) with both prefactor and decay
# rate free, on the five soft-mode-derived p_c of G_1-G_5.
pi_ = 0.94
alpha_ = 0.5
gg(x) = pi_*exp(-alpha_*x)
fit gg(x) DATA u ($0):2 via pi_, alpha_

# Reference per-Si energy of alpha-quartz (eV/Si); the offset that makes
# Delta h_quartz non-negative for the family.
E_quartz = -128.703350991

set xlabel  "{/Helvetica-Italic k} - 1 [-]"
set xrange  [-0.5 : 7.5]
set xtics   0,1,7
set ylabel  "{/Helvetica-Italic p_c} [GPa]"
set logscale y
set yrange  [0.02 : 1.05]
set ytics ("0.02" 0.02, "0.05" 0.05, "0.10" 0.1, "0.20" 0.2, "0.50" 0.5, "1.00" 1) nomirror
set y2tics 0.16,0.01,0.2 nomirror
unset mytics

# y2: zoom around the actual range of Delta h_quartz (~0.167-0.189 eV/Si).
set y2label "{/Symbol D}{/Helvetica-Italic h}_q_u_a_r_t_z / [eV/Si]"
set y2range [0.16 : 0.20]
set boxwidth 0.55 relative

set key top right reverse Left font ",18" samplen 1.5 spacing 1.1 box opaque
unset label

plot DATA u ($0):($5-E_quartz) axes x1y2 w boxes lc rgb 'blue' fs transparent solid .25 title '{/Symbol D}{/Helvetica-Italic h}_q_u_a_r_t_z',\
     gg(x) w l lw 2 lc rgb 'red' notitle,\
     DATA u ($0):2:3 w yerrorbars pt 7 ps 1.5 lc rgb 'red' title 'computed',\
     PRED u ($0+5):2 pt 6 ps 1.7 lc rgb 'red' title 'predicted'

print "================================================================"
print sprintf("Exponential fit p_c(k-1) = pi * exp(-alpha*(k-1)):")
print sprintf("  pi    = %.4f +/- %.4f GPa", pi_, pi__err)
print sprintf("  alpha = %.4f +/- %.4f",      alpha_, alpha__err)
do for [k=1:8] {
    print sprintf("  p_c(G_%d) = %.4f GPa  (k-1=%d)", k, pi_*exp(-alpha_*(k-1)), k-1)
}
print "================================================================"
