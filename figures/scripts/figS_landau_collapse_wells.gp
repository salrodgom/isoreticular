#!/usr/bin/gnuplot
# Dimensionless Landau wells fbar(phi; pi, eta) = - pi phi^2 + (1/2) phi^4
# + (eta/12) phi^6 at fixed reduced pressure pi = p/p_c(k) - 1 = 0.5
# for the five members of the RHO isoreticular family. Mean-field reference
# (eta = 0) drawn dashed. The well minima are marked with filled circles.
#
# Input:
#   ../data/figS_landau_collapse_wells.dat
#     columns: phi  fbar(eta_1=5.32)  fbar(eta_2=0.12)  fbar(eta_3=0.05)
#              fbar(eta_4<0.02)  fbar(eta_5<0.02)  fbar(eta=0, mean-field)
#   ../data/figS_landau_collapse_wells_minima.dat
#     columns: eta_k  phi_eq  fbar_eq

set term postscript eps color enhanced blacktext 'Helvetica,26'
set output 'figS_landau_collapse_wells.eps'
set encoding iso_8859_1
set locale "en_GB.UTF-8"

set xlabel "{/Symbol f} = {/Symbol D}/{/Symbol d}({/Helvetica-Italic k}) [-]"
set ylabel "2 {/Helvetica-Italic f}_{eff} / [{/Helvetica-Italic A} {/Symbol d}^2({/Helvetica-Italic k})] [-]"
set xrange [0.0 : 1.4]
set yrange [-0.20 : 0.30]
set xtics 0.2
set ytics 0.05
unset grid

set key top left reverse Left font ",20" samplen 2.5

set arrow from 0,0 to 1.4,0 nohead lc rgb 'dark-grey' lw 1

set label "{/Symbol t} = 0.5" at 0.05,0.27 left font ",22" tc rgb '#333333'

DATA = '../data/figS_landau_collapse_wells.dat'
MIN  = '../data/figS_landau_collapse_wells_minima.dat'

# Same cool-to-warm palette as the trajectory panel.
plot \
    DATA u 1:7 w l lw 2.0 dt 2 lc rgb 'dark-grey' title 'mean-field ({/Symbol h}=0)', \
    DATA u 1:2 w l lw 3.0 lc rgb '#2166ac' title '{/Helvetica-Italic k}=1 ({/Symbol h}=5.3)', \
    DATA u 1:3 w l lw 3.0 lc rgb '#67a9cf' title '{/Helvetica-Italic k}=2 ({/Symbol h}=0.12)',   \
    DATA u 1:4 w l lw 3.0 lc rgb '#f4a582' title '{/Helvetica-Italic k}=3 ({/Symbol h}=0.05)',   \
    DATA u 1:5 w l lw 3.0 lc rgb '#d6604d' title '{/Helvetica-Italic k}=4 ({/Symbol h}<0.02)',   \
    DATA u 1:6 w l lw 3.0 dt 3 lc rgb '#7b241c' title '{/Helvetica-Italic k}=5 ({/Symbol h}<0.02)', \
    MIN  u 2:3 w p pt 7 ps 1.8 lc rgb '#333333' notitle
