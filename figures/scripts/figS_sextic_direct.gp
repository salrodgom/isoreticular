#!/usr/bin/gnuplot
# Direct sextic-Landau fit Delta^2(p) = C [sqrt(1+eta(p/p_c-1))-1] on the
# broken phase for G_1-G_5. Recovers eta(k) without going through the
# inversion of beta_obs of tab:SI-beta-crossover.
#
# Inputs (produced by prep_figS_sextic_direct.py):
#   ../data/figS_sextic_direct_G{1..5}_data.dat   (p, Delta) observed
#   ../data/figS_sextic_direct_G{1..5}_fit.dat    (p, Delta) sextic fit
#   ../data/figS_sextic_direct_pars.dat           fitted (C, eta, U_eff/A, W/A)

set term postscript eps color enhanced blacktext 'Helvetica,26'
set output 'figS_sextic_direct.eps'
set encoding iso_8859_1
set locale "en_GB.UTF-8"

set xlabel "Hydrostatic pressure, {/Helvetica-Italic p} [GPa]"
set ylabel "D8R distortion, {/Symbol D} [{\305}]"
set xrange [0 : 2]
set yrange [-0.05 : 2.5]
set xtics 0.5
set ytics 0.5
unset grid
set key bottom right reverse Left font ",20" samplen 1.8 spacing 1.1 box opaque

DATA = '../data'

plot \
    DATA.'/figS_sextic_direct_G1_fit.dat'  u 1:2 w l lw 2.0 lc rgb '#cb4335' notitle, \
    DATA.'/figS_sextic_direct_G2_fit.dat'  u 1:2 w l lw 2.0 lc rgb '#f1c40f' notitle, \
    DATA.'/figS_sextic_direct_G3_fit.dat'  u 1:2 w l lw 2.0 lc rgb '#27ae60' notitle, \
    DATA.'/figS_sextic_direct_G4_fit.dat'  u 1:2 w l lw 2.0 lc rgb '#2874a6' notitle, \
    DATA.'/figS_sextic_direct_G5_fit.dat'  u 1:2 w l lw 2.0 lc rgb '#7d3c98' notitle, \
    DATA.'/figS_sextic_direct_G1_data.dat' u 1:2 w p pt 7 ps 1.0 lc rgb '#cb4335' title '{/Helvetica-Italic G}_1 ({/Symbol h}=5.3)', \
    DATA.'/figS_sextic_direct_G2_data.dat' u 1:2 w p pt 7 ps 1.0 lc rgb '#f1c40f' title '{/Helvetica-Italic G}_2 ({/Symbol h}=0.12)', \
    DATA.'/figS_sextic_direct_G3_data.dat' u 1:2 w p pt 7 ps 1.0 lc rgb '#27ae60' title '{/Helvetica-Italic G}_3 ({/Symbol h}=0.05)', \
    DATA.'/figS_sextic_direct_G4_data.dat' u 1:2 w p pt 7 ps 1.0 lc rgb '#2874a6' title '{/Helvetica-Italic G}_4 ({/Symbol h}<0.01)', \
    DATA.'/figS_sextic_direct_G5_data.dat' u 1:2 w p pt 7 ps 1.0 lc rgb '#7d3c98' title '{/Helvetica-Italic G}_5 ({/Symbol h}<0.01)'
