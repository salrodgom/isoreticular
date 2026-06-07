#!/usr/bin/gnuplot
# Figure S2(a) - conceptual effective Landau f_eff(Delta; p) vs Delta.
# Sign convention of the manuscript:
#   f_eff = (1/2) alpha (1 - p/p_c) Delta^2 + (1/4) u Delta^4
# with alpha = u = 1 (arbitrary units). Broken phase appears for p > p_c.
# Delta restricted to the physical quadrant Delta >= 0.
#
# Key (legend) shown without box border.
#
# Inputs:
#   ../data/figS2_landau_concept_a.dat          curves (5 pressures)
#   ../data/figS2_landau_concept_a_minima.dat   broken-phase minima

set term postscript eps color enhanced blacktext 'Helvetica,26'
set output 'figS2_landau_wells_concept_a.eps'
set encoding iso_8859_1
set locale "en_GB.UTF-8"
set xlabel "{/Symbol D} [a. u.]"
set ylabel "{/Helvetica-Italic f}_{eff}({/Symbol D}; {/Helvetica-Italic p}) [a. u.]"
set xrange [0.0 : 1.5]
set yrange [-0.15 : 0.5]
unset grid
set ytics 0.1
set xtics 0.5

set key top left reverse Left font ",20"
set arrow from 0,0 to 1.5,0 nohead lc rgb 'dark-grey' lw 1

DATA = '../data/figS2_landau_concept_a.dat'
MIN  = '../data/figS2_landau_concept_a_minima.dat'

plot \
    DATA u 1:2 w l lw 3.0 lc rgb '#1b4f72' title '{/Helvetica-Italic p}/{/Helvetica-Italic p}_c - 1 = -0.6',\
    DATA u 1:3 w l lw 3.0 lc rgb '#2874a6' title '{/Helvetica-Italic p}/{/Helvetica-Italic p}_c - 1 = -0.2',\
    DATA u 1:4 w l lw 3.0 lc rgb '#5dade2' title '{/Helvetica-Italic p}/{/Helvetica-Italic p}_c - 1 =  0',\
    DATA u 1:5 w l lw 3.0 lc rgb '#e67e22' title '{/Helvetica-Italic p}/{/Helvetica-Italic p}_c - 1 = +0.2',\
    DATA u 1:6 w l lw 3.0 lc rgb '#cb4335' title '{/Helvetica-Italic p}/{/Helvetica-Italic p}_c - 1 = +0.6',\
    MIN  u (column(1)==0.2 ? $2 : NaN):3 w p pt 7 ps 2.0 lc rgb '#e67e22' notitle,\
    MIN  u (column(1)==0.6 ? $2 : NaN):3 w p pt 7 ps 2.0 lc rgb '#cb4335' notitle
