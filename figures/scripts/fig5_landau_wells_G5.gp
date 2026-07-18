#!/usr/bin/gnuplot
# Main-text Figure 5(b) - effective Landau potential f_eff/A for G_5 (PST-20),
# quartic form with the coefficients of the manuscript:
#   U_eff/A = 4.640 A^-2  (delta = 0.464 A, p_c = 0.1282 GPa)
# Three pressures: 0.5 p_c, p_c, 1.5 p_c, drawn in shades of the G_5
# member colour of Figure 2 (#7d3c98), light to dark with pressure.
# Axes common with fig5_landau_wells_G1.gp so the relative shape of the
# wells is directly comparable: on this shared scale the G_5 basin is
# shallow and close to the origin, which is the point of the figure.
# Self-contained (analytic curves; no data files required).

set term postscript eps color enhanced blacktext 'Helvetica,26'
set output 'fig5_landau_wells_G5.eps'
set encoding iso_8859_1
set locale "en_GB.UTF-8"

set xlabel "{/Symbol D} [{\305}]"
set ylabel "{/Helvetica-Italic f}_{eff}({/Symbol D}; {/Helvetica-Italic p}) / {/Helvetica-Italic A} [{\305}^2]"
set xrange [0.0 : 2.5]
set yrange [-0.5 : 1.0]
unset grid
set ytics 0.25
set xtics 0.5

unset key

set arrow from 0,0 to 2.5,0 nohead lc rgb 'dark-grey' lw 1

# Quartic effective potential, U = U_eff/A for G_5
U  = 4.640
f(x,r) = 0.5*(1.0-r)*x**2 + 0.25*U*x**4
Deq = sqrt(0.5/U)
fmin = f(Deq,1.5)

set samples 1001

plot \
    f(x,0.5) w l lw 3.0 lc rgb '#c39bd3' title '{/Helvetica-Italic p} = 0.5 {/Helvetica-Italic p}_c',\
    f(x,1.0) w l lw 3.0 lc rgb '#7d3c98' title '{/Helvetica-Italic p} =      {/Helvetica-Italic p}_c',\
    f(x,1.5) w l lw 3.0 lc rgb '#4a235a' title '{/Helvetica-Italic p} = 1.5 {/Helvetica-Italic p}_c',\
    '+' u (Deq):(fmin) every ::0::0 w p pt 7 ps 2.0 lc rgb '#4a235a' notitle
