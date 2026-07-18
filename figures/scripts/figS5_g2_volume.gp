#!/usr/bin/gnuplot
# Figure S5(c) - G_2 (PST-29) volume shrinkage V/V_0 vs hydrostatic pressure,
# benchmarked against DFT r2SCAN+rVV10. Same axes as Figure S5(a).
# No legend; the colour/marker code is given in panel (a).
#
# Inputs:
#   ../data/figS5_g2_dft.dat, _slc.dat

set term postscript eps color enhanced blacktext 'Helvetica,26'
set output 'figS5_g2_volume.eps'
set encoding iso_8859_1
set locale "en_GB.UTF-8"

set xlabel "Hydrostatic pressure, {/Helvetica-Italic p} [GPa]"
set ylabel "{/Helvetica-Italic V}/{/Helvetica-Italic V}_0 [-]"
set xrange [0.0 : 2.05]
set yrange [0.70 : 1.02]
unset grid
set xtics 0.5
set ytics 0.05
unset key

DFT = '../data/figS5_g2_dft.dat'
SLC = '../data/figS5_g2_slc.dat'
NAS = '../data/figS5_g2_nasir.dat'
MAT = '../data/figS5_g2_matpes.dat'

plot \
    SLC u 1:2 w lp lw 2.0 pt 5 ps 1.4 lc rgb '#cb4335',\
    DFT u 1:2 w lp lw 2.0 pt 7 ps 1.8 lc rgb '#1b4f72',\
    NAS u 1:2 w lp lw 2.0 pt 9  ps 1.6 lc rgb '#e67e22',\
    MAT u 1:2 w lp lw 2.0 pt 13 ps 1.6 lc rgb '#117a65'
