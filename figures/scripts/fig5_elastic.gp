#!/usr/bin/gnuplot -persist
# Generates ElasticG1.eps, ElasticG2.eps, ElasticG3.eps for panels (a)-(c)
# of Figure 5 (fig:Mechanical) of the manuscript. Each EPS plots the cubic
# elastic constants C11, C12 and C44 as a function of hydrostatic pressure
# for the corresponding G_k zeolite, using the unit-cell calculations.
#
# Data columns in dir_RHO_isoreticular_G{k}_SG_P1/data_pressure_delta.txt:
#   $1  = k label
#   $2  = pressure [bar]                -> /1e4 -> GPa
#   $9  = first phonon frequency        -> filter $9>=0 keeps stable points
#   $18 = C11 [GPa]
#   $19 = C12 [GPa]
#   $20 = C44 [GPa]
#
# G_1 (RHO) is taken from the 2x2x2 SUPERCELL (dir_RHO_isoreticular_G1_222_SG_P1,
# 384 T-atoms) because the 48-atom unit cell exhibits enhanced numerical scatter
# in the broken-phase Δ near p_c (cf. Section S Finite-size in the manuscript).
# G_2 and G_3 use their unit cells, which already have hundreds of atoms.

set term postscript eps color blacktext enhanced 'Helvetica,26'
set encoding iso_8859_1
set locale "en_GB.UTF-8"
set xlabel 'Hydrostatic pressure, {/Helvetica-Italic p} / [GPa]'
set xrange [0:2]
set yrange [-10:100]
set ylabel 'Elastic constant, {/Helvetica-Italic C_{ij}} / [GPa]'
set key top right
unset grid

# Path to the raw pressure scans (relative to manuscript_figures/build/).
RAW = '../../initial_structures_RHO_isoreticular'

# Equilibrium row per pressure (single Δ per P): for each pressure keep the
# row with the smallest C_ii (the broken phase has lower C11 than the
# cubic-metastable cousin at the same P, so the line follows the physical
# branch instead of zig-zagging between the two branches).
# Pipeline (LC_ALL=C and sort -g enforce numeric sort across locales):
#   awk filters by ω₁≥0 (stable) and emits "P_GPa  Cij"
#   sort -g -k1,1 -k2,2g sorts by pressure ascending, then by Cij ascending
#   awk '!seen[$1]++' keeps only the FIRST occurrence per P (smallest Cij)
keep_eq = "| LC_ALL=C sort -g -k1,1 -k2,2g | awk '!seen[$1]++'"

do for [k=1:3] {
    set output sprintf('fig5_elastic_G%d.eps', k)
    # G_1 -> 2x2x2 supercell; G_2, G_3 -> unit cell.
    suffix = (k == 1 ? "_222_SG_P1" : "_SG_P1")
    file = sprintf("%s/dir_RHO_isoreticular_G%d%s/data_pressure_delta.txt", RAW, k, suffix)
    p \
      "< awk '{if($9>=0) print $2/10000,$18}' ".file.keep_eq w lp lw 1.5 pt 7  ps 1.2 lc rgb 'red'          t 'C_{11}',\
      "< awk '{if($9< 0) print $2/10000,$18}' ".file                  w p          pt 6  ps 1.2 lc rgb 'red'          notitle,\
      "< awk '{if($9>=0) print $2/10000,$19}' ".file.keep_eq w lp lw 1.5 pt 9  ps 1.2 lc rgb 'blue'         t 'C_{12}',\
      "< awk '{if($9< 0) print $2/10000,$19}' ".file                  w p          pt 8  ps 1.2 lc rgb 'blue'         notitle,\
      "< awk '{if($9>=0) print $2/10000,$20}' ".file.keep_eq w lp lw 1.5 pt 11 ps 1.2 lc rgb 'forest-green' t 'C_{44}',\
      "< awk '{if($9< 0) print $2/10000,$20}' ".file                  w p          pt 10 ps 1.2 lc rgb 'forest-green' notitle
}

