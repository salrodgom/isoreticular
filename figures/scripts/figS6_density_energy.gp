#!/usr/bin/gnuplot
# Figure S6 - framework density rho vs cohesive enthalpy per T-atom relative
# to alpha-quartz, Delta h_quartz, for all SLC-optimised structures of
# G_1..G_5 across the full pressure scan 0 to 2 GPa. Each point is
# colour-coded by hydrostatic pressure.
#
# The Balestra et al. (CGD 2024) parabola E(rho) = E_quartz_balestra + k_f
# (rho - FD_0)^2 anchored at G_1 is overlaid (dashed grey).
# Black asterisks mark the equilibrium (p -> 0) point of each G_k.
#
# Inputs:
#   ../data/figS6_density_energy_G{1..5}.dat   (p, rho, E)
#   ../data/figS6_density_energy_equilibrium.dat

set term postscript eps color enhanced blacktext 'Helvetica,22'
set output 'figS6_density_energy.eps'
set encoding iso_8859_1
set locale "en_GB.UTF-8"

# Alpha-quartz reference (SLC, eV/Si). Same value used in gp_isoOvsPc and
# in panel (e) of Figure 5 for the Delta h_quartz bars.
E_quartz_alpha = -128.703350991

# Balestra et al. CGD 2024 parabola anchored at the G_1 equilibrium point
# (parabola vertex at FD_0). E_quartz_balestra is the absolute energy the
# parabola predicts at FD_0 given the G_1 anchor.
E_quartz_balestra = -128.739469
k_f      = 1.3652902936e-03
FD_0     = 27.735093
# Delta h with respect to alpha-quartz:
balestra_dh(x) = (E_quartz_balestra - E_quartz_alpha) + k_f * (x - FD_0)**2

set xlabel "Framework density {/Symbol r} [T-atoms / 10^3 {\305}^3]"
set ylabel "{/Symbol D}{/Helvetica-Italic h}_{quartz} [eV/Si]"
set xrange [14.5 : 19.6]
set yrange [0.0 : 0.90]
unset grid
set xtics 1.0
set ytics 0.1

# Legend outside the plot, below (right side is reserved for the colourbar).
set key below reverse Left font ",18" samplen 1.5 maxcols 4

# Pressure colourbar (viridis-like 5-stop palette).
set palette defined ( \
    0.0 '#440154', \
    0.25 '#3b528b', \
    0.50 '#21918c', \
    0.75 '#5ec962', \
    1.0 '#fde725')
set cbrange [0.0 : 2.0]
set cblabel "Hydrostatic pressure {/Helvetica-Italic p} [GPa]"

G1 = '../data/figS6_density_energy_G1.dat'
G2 = '../data/figS6_density_energy_G2.dat'
G3 = '../data/figS6_density_energy_G3.dat'
G4 = '../data/figS6_density_energy_G4.dat'
G5 = '../data/figS6_density_energy_G5.dat'
EQ = '../data/figS6_density_energy_equilibrium.dat'

plot \
    balestra_dh(x) w l lw 1.5 dt 2 lc rgb 'dark-grey' title 'Balestra et al. (CGD 2024) parabola',\
    G1 u 2:($3 - E_quartz_alpha):1 w p pt 7  ps 1.0 palette title 'RHO ({/Helvetica-Italic G}_1)',\
    G2 u 2:($3 - E_quartz_alpha):1 w p pt 5  ps 1.0 palette title 'PWN ({/Helvetica-Italic G}_2)',\
    G3 u 2:($3 - E_quartz_alpha):1 w p pt 9  ps 1.2 palette title 'PAU ({/Helvetica-Italic G}_3)',\
    G4 u 2:($3 - E_quartz_alpha):1 w p pt 13 ps 1.2 palette title 'MWF ({/Helvetica-Italic G}_4)',\
    G5 u 2:($3 - E_quartz_alpha):1 w p pt 11 ps 1.2 palette title '{/Helvetica-Italic G}_5',\
    EQ u 1:($2 - E_quartz_alpha) w p pt 3 ps 2.0 lc rgb 'black' lw 2 title 'Equilibrium ({/Helvetica-Italic p} {/Symbol \256} 0)'
