#!/usr/bin/gnuplot
# Figure S5(d) - G_2 (PST-29) D8R distortion Delta vs hydrostatic pressure,
# benchmarked against DFT r2SCAN+rVV10. Same axes as Figure S5(b).
# No legend; the colour/marker code is given in panel (a).
#
# Inputs:
#   ../data/figS5_g2_dft.dat, _slc.dat

set term postscript eps color enhanced blacktext 'Helvetica,26'
set output 'figS5_g2_delta.eps'
set encoding iso_8859_1
set locale "en_GB.UTF-8"

set xlabel "Hydrostatic pressure {/Helvetica-Italic p} [GPa]"
set ylabel "D8R distortion {/Symbol D} [{\305}]"
set xrange [0.0 : 2.05]
set yrange [-0.05 : 3.2]
unset grid
set xtics 0.5
set ytics 0.5
unset key

DFT = '../data/figS5_g2_dft.dat'
SLC = '../data/figS5_g2_slc.dat'

plot \
    SLC u 1:3 w lp lw 2.0 pt 5 ps 1.4 lc rgb '#cb4335',\
    DFT u 1:3 w lp lw 2.0 pt 7 ps 1.8 lc rgb '#1b4f72'
