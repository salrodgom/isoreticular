#!/usr/bin/gnuplot
# Figure S9(G_2) - histogram of per-isolated-8MR distortion Delta_8 (8MR
# rings that do NOT participate in a D8R pair) for PWN (G_2), at nine
# common hydrostatic pressures across the family (0, 0.1, 0.2, 0.4, 0.6,
# 0.8, 1.0, 1.5, 2.0 GPa).
#
# This panel is the only one of Figure S9 that carries the pressure
# legend; the other three panels (G_3-G_5) reuse the same colour code.

set term postscript eps color enhanced blacktext 'Helvetica,22'
set output 'figS9_isolated_hist_G2.eps'
set encoding iso_8859_1
set locale "en_GB.UTF-8"

set xlabel "Per-isolated-8MR distortion {/Symbol D}_8 [{\305}]"
set ylabel "Probability density [{\305}^{-1}]"
set xrange [0.0 : 2.5]
set yrange [0 : 6]
unset grid
set xtics 0.5
set ytics 1.0

set key top right reverse Left font ",14" spacing 1.0 samplen 1.5

binwidth = 0.04
bin(x) = binwidth * floor(x / binwidth + 0.5)

F0   = '../data/figS9d_isolated_G2_p0bar.dat'
F01  = '../data/figS9d_isolated_G2_p1000bar.dat'
F02  = '../data/figS9d_isolated_G2_p2000bar.dat'
F04  = '../data/figS9d_isolated_G2_p4000bar.dat'
F06  = '../data/figS9d_isolated_G2_p6000bar.dat'
F08  = '../data/figS9d_isolated_G2_p8000bar.dat'
F10  = '../data/figS9d_isolated_G2_p10000bar.dat'
F15  = '../data/figS9d_isolated_G2_p15000bar.dat'
F20  = '../data/figS9d_isolated_G2_p20000bar.dat'

stats F0  u 1 nooutput; N0  = STATS_records
stats F01 u 1 nooutput; N01 = STATS_records
stats F02 u 1 nooutput; N02 = STATS_records
stats F04 u 1 nooutput; N04 = STATS_records
stats F06 u 1 nooutput; N06 = STATS_records
stats F08 u 1 nooutput; N08 = STATS_records
stats F10 u 1 nooutput; N10 = STATS_records
stats F15 u 1 nooutput; N15 = STATS_records
stats F20 u 1 nooutput; N20 = STATS_records

plot \
    F0  u (bin($1)):(1.0/N0/binwidth)  smooth freq w steps lw 2.0 lc rgb '#2166ac' title '{/Helvetica-Italic p} = 0 GPa',\
    F01 u (bin($1)):(1.0/N01/binwidth) smooth freq w steps lw 2.0 lc rgb '#4393c3' title '0.1 GPa',\
    F02 u (bin($1)):(1.0/N02/binwidth) smooth freq w steps lw 2.0 lc rgb '#92c5de' title '0.2 GPa',\
    F04 u (bin($1)):(1.0/N04/binwidth) smooth freq w steps lw 2.0 lc rgb '#76d7c4' title '0.4 GPa',\
    F06 u (bin($1)):(1.0/N06/binwidth) smooth freq w steps lw 2.0 lc rgb '#f9e79f' title '0.6 GPa',\
    F08 u (bin($1)):(1.0/N08/binwidth) smooth freq w steps lw 2.0 lc rgb '#f4a582' title '0.8 GPa',\
    F10 u (bin($1)):(1.0/N10/binwidth) smooth freq w steps lw 2.0 lc rgb '#e67e22' title '1.0 GPa',\
    F15 u (bin($1)):(1.0/N15/binwidth) smooth freq w steps lw 2.0 lc rgb '#cb4335' title '1.5 GPa',\
    F20 u (bin($1)):(1.0/N20/binwidth) smooth freq w steps lw 2.0 lc rgb '#7b241c' title '2.0 GPa'
