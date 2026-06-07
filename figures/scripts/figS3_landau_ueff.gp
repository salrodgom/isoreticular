#!/usr/bin/gnuplot
# Figure S3 - dimensionless ratio U_eff/A = 1/delta^2 vs isoreticular
# order k-1, on log scale. Quantifies the proximity to tricriticality:
#   U_eff/A -> 0    : tricritical (beta = 1/4)
#   U_eff/A large   : mean-field plateau (beta = 1/2)
#
# Annotations placed by hand: dashed reference at U_eff/A = 0.30,
# and two regime labels.
#
# Input:
#   ../data/figS3_landau_ueff.dat   columns: k-1, U_eff/A [A^-2]

set term postscript eps color enhanced blacktext 'Helvetica,26'
set output 'figS3_landau_ueff.eps'
set encoding iso_8859_1
set locale "en_GB.UTF-8"

set xlabel "{/Helvetica-Italic k} - 1 [-]"
set ylabel "{/Helvetica-Italic U}_{eff}/{/Helvetica-Italic A} = 1/{/Symbol d}^2 [{\305}^{-2}]"
set logscale y
set xrange [-0.5 : 4.5]
set yrange [0.10 : 8.0]
set xtics 0,1,4
set ytics ("0.1" 0.1, "0.2" 0.2, "0.5" 0.5, "1" 1, "2" 2, "5" 5)
unset mytics
unset grid
unset key

# Reference dashed line at U_eff/A = 0.30, conventional lower bound for
# the mean-field-compatible regime.
set arrow from -0.5,0.30 to 4.5,0.30 nohead lc rgb 'dark-grey' dt 2 lw 1

# Regime labels
set label "tricritical regime ({/Helvetica-Italic U}_{eff} {/Symbol \256} 0)" \
    at 4.0,0.16 right font ",18" tc rgb '#333333'
set label "mean-field plateau ({/Symbol b} = 1/2)" \
    at 4.0,6.0 right font ",18" tc rgb '#333333'

DATA = '../data/figS3_landau_ueff.dat'

plot DATA u 1:2 w lp lw 2.5 pt 7 ps 2.0 lc rgb '#1b4f72' notitle
