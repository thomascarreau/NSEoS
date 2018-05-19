set terminal postscript eps size 16.0,13.0 enhanced color \
font 'Helvetica,40' linewidth 3

set output 'matrix.eps'

set title "Correlation matrix for the 325 MM" font 'Helvetica, 60'

set xtics ("n_t" 0, "P_t" 1, "n_{sat}" 2, "{/Symbol l}_{sat}" 3, \
"K_{sat}" 4, "Q_{sat}" 5, "Z_{sat}" 6,"J_{sym}" 7, "L_{sym}" 8, \
"K_{sym}" 9, "Q_{sym}" 10, "Z_{sym}" 11, "m*/m" 12, "{/Symbol D}m*/m" 13, \
"b" 14, "{/Symbol s}_0" 15, "b_s" 16) font 'Helvetica,35'

set ytics ("n_t" 0, "P_t" 1, "n_{sat}" 2, "{/Symbol l}_{sat}" 3, \
"K_{sat}" 4, "Q_{sat}" 5, "Z_{sat}" 6,"J_{sym}" 7, "L_{sym}" 8, \
"K_{sym}" 9, "Q_{sym}" 10, "Z_{sym}" 11, "m*/m" 12, "{/Symbol D}m*/m" 13, \
"b" 14, "{/Symbol s}_0" 15, "b_s" 16) font 'Helvetica,35'

p 'matrix.out' matrix w image
