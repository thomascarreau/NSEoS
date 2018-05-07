set terminal postscript eps size 16.0,13.0 enhanced color \
font 'Helvetica,40' linewidth 3

set output 'matrix.eps'

set title "Correlation matrix for the 325 MM" font 'Helvetica, 60'
set xtics ("n_t" 0, "P_t" 1, "n_{sat}" 2, "{/Symbol l}_{sat}" 3, "K_{sat}" 4, "Q_{sat}" 5, "Z_{sat}" 6,"J_{sym}" 7, "L_{sym}" 8, "K_{sym}" 9, "Q_{sym}" 10, "Z_{sym}" 11, "m*/m" 12, "{/Symbol D}m*/m" 13, "b" 14)
set ytics ("n_t" 0, "P_t" 1, "n_{sat}" 2, "{/Symbol l}_{sat}" 3, "K_{sat}" 4, "Q_{sat}" 5, "Z_{sat}" 6,"J_{sym}" 7, "L_{sym}" 8, "K_{sym}" 9, "Q_{sym}" 10, "Z_{sym}" 11, "m*/m" 12, "{/Symbol D}m*/m" 13, "b" 14)

p 'matrix.out' matrix w image

reset

set terminal postscript eps size 9.0,6.0 enhanced color \
font 'Helvetica,35' linewidth 3

set output 'nt.eps'

set xlabel "n_t [fm^{-3}]"
set ylabel "Counts"

p [0.079:0.091][0:25] 'nt.histo' u ($1+$2)/2.:3 w histeps lc rgb "red" lw 3 notitle

set output 'pt.eps'

set xlabel "P_t [MeV/fm^3]"
set ylabel "Counts"

p [0.395:0.530][0:25] 'pt.histo' u ($1+$2)/2.:3 w histeps lc rgb "red" lw 3 notitle
