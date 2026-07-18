#!/usr/bin/gnuplot
# Main-text Figure 5(a) - effective Landau potential f_eff/A for G_1 (RHO),
# quartic form with the coefficients of the manuscript:
#   U_eff/A = 0.205 A^-2  (delta = 2.210 A, p_c = 0.9418 GPa)
# Three pressures: 0.5 p_c, p_c, 1.5 p_c, drawn in shades of the G_1
# member colour of Figure 2 (#cb4335), light to dark with pressure.
# Axes common with fig5_landau_wells_G5.gp so the relative shape of the
# wells is directly comparable.
# Self-contained (analytic curves; no data files required).

set term postscript eps color enhanced blacktext 'Helvetica,26'
set output 'fig5_landau_wells_G1.eps'
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

# Quartic effective potential, U = U_eff/A for G_1
U  = 0.205
f(x,r) = 0.5*(1.0-r)*x**2 + 0.25*U*x**4
# Acentric minimum at r = 1.5: Delta_eq = sqrt(0.5/U), f_min = -0.25**2/U... (computed inline)
Deq = sqrt(0.5/U)
fmin = f(Deq,1.5)

set samples 1001

plot \
    f(x,0.5) w l lw 3.0 lc rgb '#f1948a' title '{/Helvetica-Italic p} = 0.5 {/Helvetica-Italic p}_c',\
    f(x,1.0) w l lw 3.0 lc rgb '#cb4335' title '{/Helvetica-Italic p} =      {/Helvetica-Italic p}_c',\
    f(x,1.5) w l lw 3.0 lc rgb '#7b241c' title '{/Helvetica-Italic p} = 1.5 {/Helvetica-Italic p}_c',\
    '+' u (Deq):(fmin) every ::0::0 w p pt 7 ps 2.0 lc rgb '#7b241c' notitle
