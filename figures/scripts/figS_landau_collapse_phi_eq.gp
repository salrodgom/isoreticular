#!/usr/bin/gnuplot
# Dimensionless collapse of the broken-phase order parameter for the RHO
# isoreticular family. log-log of phi = Delta/delta(k) vs pi = p/p_c(k) - 1.
# Each curve corresponds to a single eta(k) value (Table tab:SI-landau):
#   eta_1 = 5.32 (G_1),  eta_2 = 0.12 (G_2),  eta_3 = 0.05 (G_3),
#   eta_4 < 0.02 (G_4), eta_5 < 0.02 (G_5).  Mean-field reference (eta = 0) shown dashed.
# Asymptotic slopes: 1/2 (mean-field) at small pi for all curves; 1/4 at large
# pi for G_1 (sextic-dominated). For G_2-G_5 the slope is essentially 1/2 over
# the full window.
#
# Input:
#   ../data/figS_landau_collapse_phi_eq.dat
#   columns: pi  phi(eta=0)  phi(eta_1=5.32)  phi(eta_2=0.12)  phi(eta_3=0.05)
#            phi(eta_4<0.02)  phi(eta_5<0.02)

set term postscript eps color enhanced blacktext 'Helvetica,26'
set output 'figS_landau_collapse_phi_eq.eps'
set encoding iso_8859_1
set locale "en_GB.UTF-8"

set xlabel "{/Symbol t} = {/Helvetica-Italic p}/{/Helvetica-Italic p}_c({/Helvetica-Italic k}) - 1 [-]"
set ylabel "{/Symbol f} = {/Symbol D}/{/Symbol d}({/Helvetica-Italic k}) [-]"
set logscale x
set logscale y
set xrange [0.01 : 10.0]
set yrange [0.08 : 4.0]
set xtics ("0.01" 0.01, "0.1" 0.1, "1" 1, "10" 10)
set ytics ("0.1" 0.1, "0.3" 0.3, "1" 1, "3" 3)
unset mytics
unset mxtics
unset grid

set key bottom right reverse Left font ",20" samplen 2.5

# Slope-reference lines (rendered to the left of each curve)
# slope 1/2 (mean-field): log(phi) = 0.5 log(pi) + const ; passes through (0.02, 0.14)
set arrow from 0.012,0.110 to 0.30,0.548 nohead lc rgb 'dark-grey' dt 2 lw 1.5
set label "slope 1/2 ({/Symbol b}=1/2)" at 0.30,0.560 left font ",18" tc rgb '#333333'
# slope 1/4 (tricritical) along the eta=160 high-pi asymptote
set arrow from 1.0,0.50 to 10.0,0.889 nohead lc rgb 'dark-grey' dt 2 lw 1.5
set label "slope 1/4 ({/Symbol b}=1/4)" at 1.05,0.46 left font ",18" tc rgb '#333333'

DATA = '../data/figS_landau_collapse_phi_eq.dat'

# Cool-to-warm gradient: G_1 deepest blue, G_5 deepest red.
plot \
    DATA u 1:2 w l lw 2.0  dt 2 lc rgb 'dark-grey' title 'mean-field ({/Symbol h}=0)', \
    DATA u 1:3 w l lw 3.0  lc rgb '#2166ac' title '{/Helvetica-Italic k}=1 ({/Symbol h}=5.3)', \
    DATA u 1:4 w l lw 3.0  lc rgb '#67a9cf' title '{/Helvetica-Italic k}=2 ({/Symbol h}=0.12)',   \
    DATA u 1:5 w l lw 3.0  lc rgb '#f4a582' title '{/Helvetica-Italic k}=3 ({/Symbol h}=0.05)',   \
    DATA u 1:6 w l lw 3.0  lc rgb '#d6604d' title '{/Helvetica-Italic k}=4 ({/Symbol h}<0.02)',   \
    DATA u 1:7 w l lw 3.0  dt 3 lc rgb '#7b241c' title '{/Helvetica-Italic k}=5 ({/Symbol h}<0.02)'
