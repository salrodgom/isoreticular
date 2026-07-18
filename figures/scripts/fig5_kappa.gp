#!/usr/bin/gnuplot -persist
# Figure 5(d) of the manuscript: isothermal compressibility κ_T(p) for the
# unit cells of G_1 to G_5. Filename CompressibilityVsPressure.eps to match
# what \includegraphics{images/CompressibilityVsPressure} loads in the .tex.
# G_1 source switched from the 2x2x2 supercell to the 1x1x1 unit cell for
# consistency with k > 1.
#
# Multiplot: the main panel shows the full κ_T(p) trace for every member
# and the zoom inset in the top-right corner expands the pre-transition
# region (p < 0.6 GPa, κ_T ∈ [-0.08, 0.03] GPa^-1) where every G_k shares
# the same baseline κ_T ≈ 0.013 GPa^-1 and G_3 develops its narrow
# negative-compressibility excursion at p ≈ 0.35 GPa.
set term postscript eps color blacktext enhanced 'Helvetica,26'
set output 'CompressibilityVsPressure.eps'
set encoding iso_8859_1
set locale "en_GB.UTF-8"

# Data pre-processing pipeline shared by main panel and inset.
keep_eq = "| LC_ALL=C sort -g -k1,1 -k2,2gr | awk '!seen[$1]++'"
stable_pos  = keep_eq . " | awk '$2>0'"
mech_unstab = keep_eq . " | awk '$2<=0'"
# Column 8 of data_pressure_delta.txt is κ_T in 1/GPa; column 9 is the
# soft-mode flag (ω₁≥0 phonon-stable, ω₁<0 elastically unstable).
# Filled circles (pt 7) = phonon- and elastically-stable (ω₁≥0 AND κ_T>0).
# Open circles  (pt 6) = phonon OR elastically unstable (ω₁<0 or κ_T<0).

set multiplot

# -------------------------------------------------------------------- main
set origin 0,0
set size 1,1
set key bottom right font ",18" samplen 2 spacing 1.1 width -1
set xlabel 'Hydrostatic pressure, {/Helvetica-Italic p} [GPa]'
set xrange [0:2]
set ylabel '{/Symbol k}_T [GPa^{-1}]'
set yrange [-0.1:0.30]
set xtics 0,0.5,2
set ytics
set arrow 1 from 0,0 to 2,0 nohead lc rgb 'gray60' dt 2 lw 0.8
unset grid

p \
  "< awk '{if($9>=0) print $2/10000,$8}' ../../initial_structures_RHO_isoreticular/dir_RHO_isoreticular_G1_222_SG_P1/data_pressure_delta.txt".stable_pos  w lp lw 1.5 pt 7 lc rgb '#cb4335' ps 1.2 t 'G_1',\
  "< awk '{if($9>=0) print $2/10000,$8}' ../../initial_structures_RHO_isoreticular/dir_RHO_isoreticular_G1_222_SG_P1/data_pressure_delta.txt".mech_unstab w p           pt 6 lc rgb '#cb4335' ps 1.2 notitle,\
  "< awk '{if($9< 0) print $2/10000,$8}' ../../initial_structures_RHO_isoreticular/dir_RHO_isoreticular_G1_222_SG_P1/data_pressure_delta.txt"             w p           pt 6 lc rgb '#cb4335' ps 1.2 notitle,\
  "< awk '{if($9>=0) print $2/10000,$8}' ../../initial_structures_RHO_isoreticular/dir_RHO_isoreticular_G2_SG_P1/data_pressure_delta.txt".stable_pos      w lp lw 1.5 pt 7 lc rgb '#f1c40f' ps 1.2 t 'G_2',\
  "< awk '{if($9>=0) print $2/10000,$8}' ../../initial_structures_RHO_isoreticular/dir_RHO_isoreticular_G2_SG_P1/data_pressure_delta.txt".mech_unstab     w p           pt 6 lc rgb '#f1c40f' ps 1.2 notitle,\
  "< awk '{if($9< 0) print $2/10000,$8}' ../../initial_structures_RHO_isoreticular/dir_RHO_isoreticular_G2_SG_P1/data_pressure_delta.txt"                 w p           pt 6 lc rgb '#f1c40f' ps 1.2 notitle,\
  "< awk '{if($9>=0) print $2/10000,$8}' ../../initial_structures_RHO_isoreticular/dir_RHO_isoreticular_G3_SG_P1/data_pressure_delta.txt".stable_pos      w lp lw 1.5 pt 7 lc rgb '#27ae60' ps 1.2 t 'G_3',\
  "< awk '{if($9>=0) print $2/10000,$8}' ../../initial_structures_RHO_isoreticular/dir_RHO_isoreticular_G3_SG_P1/data_pressure_delta.txt".mech_unstab     w p           pt 6 lc rgb '#27ae60' ps 1.2 notitle,\
  "< awk '{if($9< 0) print $2/10000,$8}' ../../initial_structures_RHO_isoreticular/dir_RHO_isoreticular_G3_SG_P1/data_pressure_delta.txt"                 w p           pt 6 lc rgb '#27ae60' ps 1.2 notitle,\
  "< awk '{if($9>=0) print $2/10000,$8}' ../../initial_structures_RHO_isoreticular/dir_RHO_isoreticular_G4_SG_P1/data_pressure_delta.txt".stable_pos      w lp lw 1.5 pt 7 lc rgb '#2874a6' ps 1.2 t 'G_4',\
  "< awk '{if($9>=0) print $2/10000,$8}' ../../initial_structures_RHO_isoreticular/dir_RHO_isoreticular_G4_SG_P1/data_pressure_delta.txt".mech_unstab     w p           pt 6 lc rgb '#2874a6' ps 1.2 notitle,\
  "< awk '{if($9< 0) print $2/10000,$8}' ../../initial_structures_RHO_isoreticular/dir_RHO_isoreticular_G4_SG_P1/data_pressure_delta.txt"                 w p           pt 6 lc rgb '#2874a6' ps 1.2 notitle,\
  "< awk '{if($9>=0) print $2/10000,$8}' ../../initial_structures_RHO_isoreticular/dir_RHO_isoreticular_G5_SG_P1/data_pressure_delta.txt".stable_pos      w lp lw 1.5 pt 7 lc rgb '#7d3c98' ps 1.2 t 'G_5',\
  "< awk '{if($9>=0) print $2/10000,$8}' ../../initial_structures_RHO_isoreticular/dir_RHO_isoreticular_G5_SG_P1/data_pressure_delta.txt".mech_unstab     w p           pt 6 lc rgb '#7d3c98' ps 1.2 notitle,\
  "< awk '{if($9< 0) print $2/10000,$8}' ../../initial_structures_RHO_isoreticular/dir_RHO_isoreticular_G5_SG_P1/data_pressure_delta.txt"                 w p           pt 6 lc rgb '#7d3c98' ps 1.2 notitle

# ------------------------------------------------------------------- inset
# Zoom on the pre-transition baseline (p < 0.6 GPa). Same data, same colour
# code, no legend. Size ≈ 32% × 32% of the canvas (~80% of a "default"
# matplotlib-style inset), positioned in the top-right corner where the
# main κ_T trace is empty.
unset arrow 1
set origin 0.56, 0.5
set size 0.40, 0.45
unset key
set xlabel ""
set ylabel ""
set xrange [0.25:0.5]
set yrange [-0.1:0.25]
set xtics font ",16"
set ytics ("-0.1" -0.1, "0" 0, "0.25" 0.25) font ",16"
set arrow 2 from 0,0 to 0.6,0 nohead lc rgb 'gray60' dt 2 lw 0.5
p \
  "< awk '{if($9>=0) print $2/10000,$8}' ../../initial_structures_RHO_isoreticular/dir_RHO_isoreticular_G3_SG_P1/data_pressure_delta.txt".stable_pos      w lp lw 1.0 pt 7 lc rgb '#27ae60' ps 0.7 notitle,\
  "< awk '{if($9>=0) print $2/10000,$8}' ../../initial_structures_RHO_isoreticular/dir_RHO_isoreticular_G3_SG_P1/data_pressure_delta.txt".mech_unstab     w p           pt 6 lc rgb '#27ae60' ps 0.7 notitle,\
  "< awk '{if($9< 0) print $2/10000,$8}' ../../initial_structures_RHO_isoreticular/dir_RHO_isoreticular_G3_SG_P1/data_pressure_delta.txt"                 w p           pt 6 lc rgb '#27ae60' ps 0.7 notitle

unset multiplot
