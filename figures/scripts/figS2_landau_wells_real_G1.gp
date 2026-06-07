#!/usr/bin/gnuplot
# Figure S2(c) - effective Landau f_eff/A for G_1 (RHO) with the
# parameters extracted in tab:SI-landau of the manuscript:
#   delta = 2.210 A, U_eff/A = 0.205 A^-2, p_c = 0.9418 GPa
# Plotted at three pressures: 0.5 p_c (cubic), p_c (transition), 1.5 p_c
# (broken). Delta restricted to the physical quadrant Delta >= 0.
#
# Key (legend) shown without box border.
#
# Inputs:
#   ../data/figS2_landau_real_G1.dat
#   ../data/figS2_landau_real_G1_minima.dat

set term postscript eps color enhanced blacktext 'Helvetica,26'
set output 'figS2_landau_wells_real_G1.eps'
set encoding iso_8859_1
set locale "en_GB.UTF-8"

set xlabel "{/Symbol D} [{\305}]"
set ylabel "{/Helvetica-Italic f}_{eff}({/Symbol D}; {/Helvetica-Italic p}) / {/Helvetica-Italic A} [{\305}^2]"
set xrange [0.0 : 2.5]
set yrange [-0.5 : 1.0]
unset grid
set ytics 0.25
set xtics 0.5

set key top left reverse Left font ",20"

set arrow from 0,0 to 2.5,0 nohead lc rgb 'dark-grey' lw 1

DATA = '../data/figS2_landau_real_G1.dat'
MIN  = '../data/figS2_landau_real_G1_minima.dat'

plot \
    DATA u 1:2 w l lw 3.0 lc rgb '#1b4f72' title '{/Helvetica-Italic p} = 0.5 {/Helvetica-Italic p}_c',\
    DATA u 1:3 w l lw 3.0 lc rgb '#5dade2' title '{/Helvetica-Italic p} =      {/Helvetica-Italic p}_c',\
    DATA u 1:4 w l lw 3.0 lc rgb '#cb4335' title '{/Helvetica-Italic p} = 1.5 {/Helvetica-Italic p}_c',\
    MIN  u 2:3 w p pt 7 ps 2.0 lc rgb '#cb4335' notitle
