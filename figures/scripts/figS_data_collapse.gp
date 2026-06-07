#!/usr/bin/gnuplot
# Empirical data collapse for the RHO isoreticular family,
# Figure fig:isocriticalpressure2 of the SI:
#     x = |p/p_c(k) - 1|^beta(k)
#     y = Delta(p) / delta(k)
# Each member is rescaled by its own (p_c, delta, beta) taken from
# Table tab:SI-landau / tab:SI-beta-crossover of the manuscript:
#
#    G_k    p_c [GPa]    delta [A]    beta
#     1     0.9418        2.210       0.355
#     2     0.5338        1.236       0.455
#     3     0.3627        1.014       0.429
#     4     0.2241        0.744       0.480
#     5     0.1282        0.464       0.554
#
# If the Heaviside-power phenomenology of Equation \ref{eq:fit} is exact for
# every member, the broken-phase points fall on the diagonal y = x and the
# different frameworks collapse onto a single straight line. The diagonal
# is drawn in grey for reference.
#
# Data: raw GULP/PLUMED pressure scans under
#   ../../initial_structures_RHO_isoreticular/dir_RHO_isoreticular_G<k>_*/data_pressure_delta.txt
# Columns: $1=k label, $2=p [bar], $4=Delta [nm] (->*10 [A]), $9=omega_1.
# Stable rows (omega_1>=0 or NF<9 short relaxations) are kept; rows with
# p == p_c are dropped (singular log point).

set term postscript eps color enhanced blacktext 'Helvetica,26'
set output 'figS_data_collapse.eps'
set encoding iso_8859_1
set locale "en_GB.UTF-8"

# Member-specific scales (must match tab:SI-landau / tab:SI-beta-crossover)
x01 = 0.9418 ; d01 = 2.210 ; b01 = 0.355
x02 = 0.5338 ; d02 = 1.236 ; b02 = 0.455
x03 = 0.3627 ; d03 = 1.014 ; b03 = 0.429
x04 = 0.2241 ; d04 = 0.744 ; b04 = 0.480
x05 = 0.1282 ; d05 = 0.464 ; b05 = 0.554

RAW = '../../initial_structures_RHO_isoreticular'
STABLE_FILT = "(NF<9 || $9>=0)"

set xlabel "{/Symbol \174}{/Helvetica-Italic p}/{/Helvetica-Italic p}_c({/Helvetica-Italic k}) - 1{/Symbol \174}^{ {/Symbol b}({/Helvetica-Italic k})} [-]"
set ylabel "{/Symbol D}/{/Symbol d}({/Helvetica-Italic k}) [-]"
set xrange [0.05 : 3.0]
set yrange [0.05 : 3.0]
set logscale x
set logscale y
set xtics ("0.1" 0.1, "0.3" 0.3, "1" 1, "3" 3)
set ytics ("0.1" 0.1, "0.3" 0.3, "1" 1, "3" 3)
unset mxtics
unset mytics
unset grid

set key bottom right reverse Left font ",20" samplen 1.8 spacing 1.1 box opaque

# Diagonal y = x reference: Delta/delta = |p/p_c-1|^beta is the
# Heaviside-power prediction. Members that obey it land on the diagonal.
set samples 1001
diag(x) = x

# Helper functions for the awk filters (returned as strings).
xfun(pc, b) = sprintf("( ((($2/10000)/%.5f)-1)>0 ? ((($2/10000)/%.5f)-1)**%.4f : (1-(($2/10000)/%.5f))**%.4f )", pc, pc, b, pc, b)
yfun(d) = sprintf("$4*10/%.4f", d)

plot \
    diag(x) w l lw 1.5 dt 2 lc rgb 'dark-grey' title 'diagonal {/Symbol D}/{/Symbol d} = {/Symbol \174}{/Helvetica-Italic p/p}_c-1{/Symbol \174}^{/Symbol b}', \
    "< awk '{if(".STABLE_FILT." && ($2/10000)!=".sprintf("%.5f", x01).") print ".xfun(x01,b01)." , ".yfun(d01)." }' ".RAW."/dir_RHO_isoreticular_G1_222_SG_P1/data_pressure_delta.txt" \
        w p pt 7  ps 1.2 lc rgb '#cb4335' title '{/Helvetica-Italic G}_1', \
    "< awk '{if(".STABLE_FILT." && ($2/10000)!=".sprintf("%.5f", x02).") print ".xfun(x02,b02)." , ".yfun(d02)." }' ".RAW."/dir_RHO_isoreticular_G2_SG_P1/data_pressure_delta.txt" \
        w p pt 9  ps 1.4 lc rgb '#f1c40f' title '{/Helvetica-Italic G}_2', \
    "< awk '{if(".STABLE_FILT." && ($2/10000)!=".sprintf("%.5f", x03).") print ".xfun(x03,b03)." , ".yfun(d03)." }' ".RAW."/dir_RHO_isoreticular_G3_SG_P1/data_pressure_delta.txt" \
        w p pt 11 ps 1.4 lc rgb '#27ae60' title '{/Helvetica-Italic G}_3', \
    "< awk '{if(".STABLE_FILT." && ($2/10000)!=".sprintf("%.5f", x04).") print ".xfun(x04,b04)." , ".yfun(d04)." }' ".RAW."/dir_RHO_isoreticular_G4_SG_P1/data_pressure_delta.txt" \
        w p pt 13 ps 1.4 lc rgb '#2874a6' title '{/Helvetica-Italic G}_4', \
    "< awk '{if(".STABLE_FILT." && ($2/10000)!=".sprintf("%.5f", x05).") print ".xfun(x05,b05)." , ".yfun(d05)." }' ".RAW."/dir_RHO_isoreticular_G5_SG_P1/data_pressure_delta.txt" \
        w p pt 5  ps 1.2 lc rgb '#7d3c98' title '{/Helvetica-Italic G}_5'
