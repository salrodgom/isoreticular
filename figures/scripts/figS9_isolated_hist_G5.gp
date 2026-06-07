#!/usr/bin/gnuplot
# Figure S9(G_5) - histogram of per-isolated-8MR distortion Delta_8 at
# nine common pressures (0, 0.1, 0.2, 0.4, 0.6, 0.8, 1.0, 1.5, 2.0 GPa).
# No legend: shared with panel G_2.

set term postscript eps color enhanced blacktext 'Helvetica,22'
set output 'figS9_isolated_hist_G5.eps'
set encoding iso_8859_1
set locale "en_GB.UTF-8"

set xlabel "Per-isolated-8MR distortion {/Symbol D}_8 [{\305}]"
set ylabel "Probability density [{\305}^{-1}]"
set xrange [0.0 : 2.5]
set yrange [0 : 6]
unset grid
set xtics 0.5
set ytics 1.0
unset key

binwidth = 0.04
bin(x) = binwidth * floor(x / binwidth + 0.5)

F0   = '../data/figS9d_isolated_G5_p0bar.dat'
F01  = '../data/figS9d_isolated_G5_p1000bar.dat'
F02  = '../data/figS9d_isolated_G5_p2000bar.dat'
F04  = '../data/figS9d_isolated_G5_p4000bar.dat'
F06  = '../data/figS9d_isolated_G5_p6000bar.dat'
F08  = '../data/figS9d_isolated_G5_p8000bar.dat'
F10  = '../data/figS9d_isolated_G5_p10000bar.dat'
F15  = '../data/figS9d_isolated_G5_p15000bar.dat'
F20  = '../data/figS9d_isolated_G5_p20000bar.dat'

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
    F0  u (bin($1)):(1.0/N0/binwidth)  smooth freq w steps lw 2.0 lc rgb '#2166ac',\
    F01 u (bin($1)):(1.0/N01/binwidth) smooth freq w steps lw 2.0 lc rgb '#4393c3',\
    F02 u (bin($1)):(1.0/N02/binwidth) smooth freq w steps lw 2.0 lc rgb '#92c5de',\
    F04 u (bin($1)):(1.0/N04/binwidth) smooth freq w steps lw 2.0 lc rgb '#76d7c4',\
    F06 u (bin($1)):(1.0/N06/binwidth) smooth freq w steps lw 2.0 lc rgb '#f9e79f',\
    F08 u (bin($1)):(1.0/N08/binwidth) smooth freq w steps lw 2.0 lc rgb '#f4a582',\
    F10 u (bin($1)):(1.0/N10/binwidth) smooth freq w steps lw 2.0 lc rgb '#e67e22',\
    F15 u (bin($1)):(1.0/N15/binwidth) smooth freq w steps lw 2.0 lc rgb '#cb4335',\
    F20 u (bin($1)):(1.0/N20/binwidth) smooth freq w steps lw 2.0 lc rgb '#7b241c'
