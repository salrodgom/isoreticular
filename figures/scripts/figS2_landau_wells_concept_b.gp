#!/usr/bin/gnuplot
# Figure S2(b) - strain well of the conceptual Landau functional at fixed Delta.
# Plots f(Delta, v; p) - f(Delta, 1; p) vs v, at p/p_c - 1 = +0.3, for four
# fixed Delta values. v restricted to the physical compression range v <= 1.
#
# Key (legend) shown without box border.
#
# Inputs:
#   ../data/figS2_landau_concept_b.dat          four curves
#   ../data/figS2_landau_concept_b_minima.dat   v_min markers

set term postscript eps color enhanced blacktext 'Helvetica,26'
set output 'figS2_landau_wells_concept_b.eps'
set encoding iso_8859_1
set locale "en_GB.UTF-8"

set xlabel "{/Helvetica-Italic v} = {/Helvetica-Italic V}/{/Helvetica-Italic V}_0 [-]"
set ylabel "{/Helvetica-Italic f}({/Symbol D}, {/Helvetica-Italic v}; {/Helvetica-Italic p}) - {/Helvetica-Italic f}({/Symbol D}, 1; {/Helvetica-Italic p}) [a. u.]"
set xrange [0.5 : 1.0]
set yrange [-0.10 : 0.01]
unset grid

set key bottom right font ",20"
set arrow from 0.5,0 to 1.0,0 nohead lc rgb 'dark-grey' lw 1
set arrow from 1.0,-0.10 to 1.0,0.10 nohead lc rgb '#a6a6a6' dt 3 lw 0.8

DATA = '../data/figS2_landau_concept_b.dat'
MIN  = '../data/figS2_landau_concept_b_minima.dat'

plot \
    DATA u 1:2 w l lw 3.0 lc rgb '#1b4f72' title '{/Symbol D} = 0.0',\
    DATA u 1:3 w l lw 3.0 lc rgb '#5dade2' title '{/Symbol D} = 0.4',\
    DATA u 1:4 w l lw 3.0 lc rgb '#e67e22' title '{/Symbol D} = 0.8',\
    DATA u 1:5 w l lw 3.0 lc rgb '#cb4335' title '{/Symbol D} = 1.2',\
    MIN  u (column(1)==0.0 ? $2 : NaN):3 w p pt 7 ps 2.0 lc rgb '#1b4f72' notitle,\
    MIN  u (column(1)==0.4 ? $2 : NaN):3 w p pt 7 ps 2.0 lc rgb '#5dade2' notitle,\
    MIN  u (column(1)==0.8 ? $2 : NaN):3 w p pt 7 ps 2.0 lc rgb '#e67e22' notitle,\
    MIN  u (column(1)==1.2 ? $2 : NaN):3 w p pt 7 ps 2.0 lc rgb '#cb4335' notitle
