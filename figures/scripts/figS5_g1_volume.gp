#!/usr/bin/gnuplot
# Figure S5(a) - G_1 (RHO) volume shrinkage V/V_0 vs hydrostatic pressure,
# benchmarked across four levels of theory.
#
# This is the only panel of Figure S5 that carries a legend; the other three
# panels reuse the same colour/marker code.
#
# Inputs:
#   ../data/figS5_g1_dft.dat        DFT r2SCAN+rVV10
#   ../data/figS5_g1_slc.dat        SLC (this work)
#   ../data/figS5_g1_nasir.dat      MACE (Nasir et al.)
#   ../data/figS5_g1_matpes.dat     MACE-MP foundation (MatPES r2SCAN)

set term postscript eps color enhanced blacktext 'Helvetica,26'
set output 'figS5_g1_volume.eps'
set encoding iso_8859_1
set locale "en_GB.UTF-8"

set xlabel "Hydrostatic pressure, {/Helvetica-Italic p} [GPa]"
set ylabel "{/Helvetica-Italic V}/{/Helvetica-Italic V}_0 [-]"
set xrange [0.0 : 2.05]
set yrange [0.70 : 1.02]
unset grid
set xtics 0.5
set ytics 0.05

set key bottom left reverse Left font ",20"

DFT = '../data/figS5_g1_dft.dat'
SLC = '../data/figS5_g1_slc.dat'
NAS = '../data/figS5_g1_nasir.dat'
MAT = '../data/figS5_g1_matpes.dat'

plot \
    SLC u 1:2 w lp lw 2.0 pt 5  ps 1.4 lc rgb '#cb4335' title 'FF/SLC',\
    DFT u 1:2 w lp lw 2.0 pt 7  ps 1.6 lc rgb '#1b4f72' title 'DFT/r^2SCAN+rVV10',\
    NAS u 1:2 w lp lw 2.0 pt 9  ps 1.6 lc rgb '#e67e22' title 'MACE/Nasir',\
    MAT u 1:2 w lp lw 2.0 pt 13 ps 1.6 lc rgb '#117a65' title 'MACE-MATPES-0/r^2SCAN'
