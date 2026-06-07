#!/usr/bin/gnuplot
# Figure S4 - G_1 (RHO) finite-size benchmark of the SLC pressure scan:
# D8R distortion Delta vs hydrostatic pressure p, computed in the cubic
# 1x1x1 unit cell (48 T-atoms, red) and in the 2x2x2 supercell (384 T-atoms,
# blue). After densifying the pressure grid in the critical region, both
# cell sizes converge to p_c = 0.94 GPa within 4 MPa.
#
# Fit form (manuscript Equation): Delta = delta_eq * (p/p_c - 1)^beta
#   1x1x1:   p_c = 0.944 GPa,  beta = 0.40,  delta_amp = 2.250 A
#   2x2x2:   p_c = 0.940 GPa,  beta = 0.358, delta_amp = 2.213 A
#
# Key (legend) shown without box border.
#
# Inputs:
#   ../data/figS4_g1_finite_size_1x1x1_stable.dat
#   ../data/figS4_g1_finite_size_1x1x1_unstable.dat
#   ../data/figS4_g1_finite_size_2x2x2.dat

set term postscript eps color enhanced blacktext 'Helvetica,24'
set output 'figS4_g1_finite_size.eps'
set encoding iso_8859_1
set locale "en_GB.UTF-8"

set xlabel "Hydrostatic pressure {/Helvetica-Italic p} [GPa]"
set ylabel "D8R distortion {/Symbol D} [{\305}]"
set xrange [0.0 : 2.05]
set yrange [-0.1 : 2.7]
unset grid
set xtics 0.5
set ytics 0.5

set key top left reverse Left font ",18"

# Fit curves; restricted to the broken branch p >= p_c.
f1(x) = (x >= 0.944) ? 2.250 * (x/0.944 - 1.0)**0.399 : 1/0
f2(x) = (x >= 0.940) ? 2.213 * (x/0.940 - 1.0)**0.358 : 1/0

# Light grey baseline at Delta = 0 to make the cubic branch points stand out.
set arrow from 0,0 to 2.05,0 nohead lc rgb 'dark-grey' lw 1

D11 = '../data/figS4_g1_finite_size_1x1x1_stable.dat'
D11U = '../data/figS4_g1_finite_size_1x1x1_unstable.dat'
D22 = '../data/figS4_g1_finite_size_2x2x2.dat'

plot \
    f1(x) w l lw 2.0 lc rgb '#c0392b' title '1{/Symbol \264}1{/Symbol \264}1 fit  ({/Helvetica-Italic p}_c = 0.944 GPa, {/Symbol b} = 0.40)',\
    f2(x) w l lw 2.0 dt 2 lc rgb '#1f4e79' title '2{/Symbol \264}2{/Symbol \264}2 fit  ({/Helvetica-Italic p}_c = 0.940 GPa, {/Symbol b} = 0.36)',\
    D11  u 1:2 w p pt 7 ps 1.2 lc rgb '#c0392b' title '1{/Symbol \264}1{/Symbol \264}1 data ({/Symbol w}_1 {/Symbol \263} 0)',\
    D11U u 1:2 w p pt 6 ps 1.2 lc rgb '#c0392b' title '1{/Symbol \264}1{/Symbol \264}1 data ({/Symbol w}_1 < 0)',\
    D22  u 1:2 w p pt 4 ps 1.4 lw 2 lc rgb '#1f4e79' title '2{/Symbol \264}2{/Symbol \264}2 data'
