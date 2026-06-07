#!/usr/bin/gnuplot
# Figure 3, panels (a)-(d) of the manuscript:
# Pressure evolution of the six lowest Gamma-point phonon branches
# for G_1 to G_4 in a 2x2 multipanel EPS (Nature Materials style).
#
# All branches are plotted as filled circles colour-coded by the D8R
# distortion Delta. The colour scale distinguishes centric branch
# (Delta = 0, dark purple) from acentric branch (Delta > 0,
# blue-green-yellow). Vertical dashed line at p_c; horizontal dashed
# line at omega = 0.
#
# Per-row data structure of data_pressure_delta.txt:
#   col  2 = pressure (bar)
#   col  4 = D8R distortion Delta (nm)
#   cols 9-11  = three lowest IMAGINARY frequencies (cm^-1, negative)
#   cols 12-17 = six lowest REAL frequencies (cm^-1, positive)
#
# Inputs (raw pressure scans):
#   ../../initial_structures_RHO_isoreticular/dir_RHO_isoreticular_G<k>_*/data_pressure_delta.txt

set term postscript eps color enhanced blacktext 'Helvetica,16' size 7.5,5.5
set output 'fig3_phonons.eps'
set encoding iso_8859_1
set locale "en_GB.UTF-8"

RAW = '../../initial_structures_RHO_isoreticular'

# --- p_c per member (soft-mode values, identical to Table 1) ---
pc_1 = 0.9418
pc_2 = 0.5338
pc_3 = 0.3627
pc_4 = 0.2241

# --- Awk filter: emit up to 6 (p, omega, Delta) points per data row,
#     pairing the three lowest "slots" with their imaginary counterparts
#     when the parent-symmetry calculation reports an imaginary mode.
F_ALL_MODES = "awk 'NF>=17 { p=$2/1e4; D=$4*10; \
    if ($12>0.5)      print p, $12, D; else if ($9<-0.5)  print p, $9,  D; \
    if ($13>0.5)      print p, $13, D; else if ($10<-0.5) print p, $10, D; \
    if ($14>0.5)      print p, $14, D; else if ($11<-0.5) print p, $11, D; \
    if ($15>0.5)      print p, $15, D; \
    if ($16>0.5)      print p, $16, D; \
    if ($17>0.5)      print p, $17, D }'"

# --- Shared aesthetics ---
# Viridis-style perceptually-uniform colormap, dark at low Delta for good
# contrast against the white background, bright at high Delta.
set palette defined (\
    0    '#440154',\
    0.25 '#3b528b',\
    0.50 '#21918c',\
    1.00 '#5ec962',\
    2.00 '#fde725')
set cbrange [0 : 2.5]
set cblabel "{/Symbol D} [{\305}]" offset 1.5,0
unset key

set style line 11 lt 1 lc rgb '#888888' dt 2 lw 1.2   # vertical p_c marker
set style line 12 lt 1 lc rgb '#bbbbbb' dt 3 lw 0.8   # horizontal omega=0
set tics nomirror
# Full rectangular box around each panel: bottom (1) + left (2) + top (4)
# + right (8) = 15. The Nature Materials style closes the plot frame.
set border 15 lw 1.0

# Opaque white textbox so the per-panel labels (a, b, c, d) stay legible
# even where data points cluster near the top of the panel.
set style textbox opaque fc rgb 'white' noborder margins 0.5,0.5

# Margins for tight Nature Mat. style 2x2 grid. The layout fills the EPS
# from LM=0.08 to the right side of the colorbar, and from BM=0.10 to
# TM=0.96, with no empty space at the bottom-right.
LM = 0.08
BM = 0.10
TM = 0.96
RM_PLOT = 0.87       # right edge of the plot panels
CB_LEFT = 0.90       # left edge of the colorbar
CB_WIDTH = 0.025
W  = (RM_PLOT - LM) * 0.485
H  = (TM - BM) * 0.43
GAP_X = 0.03
GAP_Y = 0.13

set multiplot

# ============================================================
# Panel (a) G_1 (RHO)
# ============================================================
set lmargin at screen LM
set rmargin at screen LM + W
set tmargin at screen TM
set bmargin at screen TM - H

set xrange [0 : 2.0]
set yrange [-15 : 55]
set xtics 0.5 format ""
set ytics 10 format "%g"
set ylabel "{/Symbol w} [cm^{-1}]" offset 1,0
unset xlabel
unset colorbox

set arrow 1 from pc_1, -15 to pc_1, 55 nohead ls 11
set arrow 2 from 0, 0 to 2.0, 0 nohead ls 12

set label 1 "{/Helvetica-Bold a}  RHO ({/Helvetica-Italic G}_1)" at graph 0.04, 0.92 font ",15" boxed front

plot \
    "< ".F_ALL_MODES." ".RAW."/dir_RHO_isoreticular_G1_222_SG_P1/data_pressure_delta.txt" u 1:2:3 w p pt 7 ps 0.85 lc palette

# ============================================================
# Panel (b) G_2 (PST-29)
# ============================================================
set lmargin at screen LM + W + GAP_X
set rmargin at screen LM + 2*W + GAP_X
set tmargin at screen TM
set bmargin at screen TM - H

set xrange [0 : 2.0]
set yrange [-15 : 55]
set xtics 0.5 format ""
set ytics 10 format ""
unset ylabel
unset xlabel

unset arrow 1
unset arrow 2
set arrow 1 from pc_2, -15 to pc_2, 55 nohead ls 11
set arrow 2 from 0, 0 to 2.0, 0 nohead ls 12

unset label 1
set label 1 "{/Helvetica-Bold b}  PST-29 ({/Helvetica-Italic G}_2)" at graph 0.04, 0.92 font ",15" boxed front

# Add the shared colorbar (full vertical span of both rows of panels)
set colorbox vertical user origin CB_LEFT, (BM + GAP_Y) size CB_WIDTH, (H + GAP_Y + H - GAP_Y)

plot \
    "< ".F_ALL_MODES." ".RAW."/dir_RHO_isoreticular_G2_SG_P1/data_pressure_delta.txt" u 1:2:3 w p pt 7 ps 0.85 lc palette

unset colorbox

# ============================================================
# Panel (c) G_3 (PAU)
# ============================================================
set lmargin at screen LM
set rmargin at screen LM + W
set tmargin at screen BM + H + GAP_Y
set bmargin at screen BM + GAP_Y

set xrange [0 : 2.0]
set yrange [-15 : 45]
set xtics 0.5 format "%.1f"
set ytics 10 format "%g"
set ylabel "{/Symbol w} [cm^{-1}]" offset 1,0
set xlabel "{/Helvetica-Italic p} [GPa]" offset 0,0.3

unset arrow 1
unset arrow 2
set arrow 1 from pc_3, -15 to pc_3, 45 nohead ls 11
set arrow 2 from 0, 0 to 2.0, 0 nohead ls 12

unset label 1
set label 1 "{/Helvetica-Bold c}  PAU ({/Helvetica-Italic G}_3)" at graph 0.04, 0.92 font ",15" boxed front

plot \
    "< ".F_ALL_MODES." ".RAW."/dir_RHO_isoreticular_G3_SG_P1/data_pressure_delta.txt" u 1:2:3 w p pt 7 ps 0.85 lc palette

# ============================================================
# Panel (d) G_4 (ZSM-25)
# ============================================================
set lmargin at screen LM + W + GAP_X
set rmargin at screen LM + 2*W + GAP_X
set tmargin at screen BM + H + GAP_Y
set bmargin at screen BM + GAP_Y

set xrange [0 : 2.0]
set yrange [-15 : 35]
set xtics 0.5 format "%.1f"
set ytics 10 format ""
unset ylabel
set xlabel "{/Helvetica-Italic p} [GPa]" offset 0,0.3

unset arrow 1
unset arrow 2
set arrow 1 from pc_4, -15 to pc_4, 35 nohead ls 11
set arrow 2 from 0, 0 to 2.0, 0 nohead ls 12

unset label 1
set label 1 "{/Helvetica-Bold d}  ZSM-25 ({/Helvetica-Italic G}_4)" at graph 0.04, 0.92 font ",15" boxed front

plot \
    "< ".F_ALL_MODES." ".RAW."/dir_RHO_isoreticular_G4_SG_P1/data_pressure_delta.txt" u 1:2:3 w p pt 7 ps 0.85 lc palette

unset multiplot
