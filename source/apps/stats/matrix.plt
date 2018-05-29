set terminal postscript eps size 16.0,13.0 enhanced color \
font 'Helvetica,40' linewidth 3

set xtics ("n_{sat}" 0, "{/Symbol l}_{sat}" 1, \
"K_{sat}" 2, "Q_{sat}" 3, "Z_{sat}" 4,"J_{sym}" 5, "L_{sym}" 6, \
"K_{sym}" 7, "Q_{sym}" 8, "Z_{sym}" 9, "m*/m" 10, "{/Symbol D}m*/m" 11, \
"b" 12, "{/Symbol s}_0" 13, "b_s" 14) font 'Helvetica,35'

set ytics ("n_{sat}" 0, "{/Symbol l}_{sat}" 1, \
"K_{sat}" 2, "Q_{sat}" 3, "Z_{sat}" 4,"J_{sym}" 5, "L_{sym}" 6, \
"K_{sym}" 7, "Q_{sym}" 8, "Z_{sym}" 9, "m*/m" 10, "{/Symbol D}m*/m" 11, \
"b" 12, "{/Symbol s}_0" 13, "b_s" 14) font 'Helvetica,35'

# set xtics ("n_t" 0, "P_t" 1, "n_{sat}" 2, "{/Symbol l}_{sat}" 3, \
# "K_{sat}" 4, "Q_{sat}" 5, "Z_{sat}" 6,"J_{sym}" 7, "L_{sym}" 8, \
# "K_{sym}" 9, "Q_{sym}" 10, "Z_{sym}" 11, "m*/m" 12, "{/Symbol D}m*/m" 13, \
# "b" 14, "{/Symbol s}_0" 15, "b_s" 16) font 'Helvetica,35'

# set ytics ("n_t" 0, "P_t" 1, "n_{sat}" 2, "{/Symbol l}_{sat}" 3, \
# "K_{sat}" 4, "Q_{sat}" 5, "Z_{sat}" 6,"J_{sym}" 7, "L_{sym}" 8, \
# "K_{sym}" 9, "Q_{sym}" 10, "Z_{sym}" 11, "m*/m" 12, "{/Symbol D}m*/m" 13, \
# "b" 14, "{/Symbol s}_0" 15, "b_s" 16) font 'Helvetica,35'

set palette defined (0 'white', 1 'green')

set output 'matrix.eps'
p 'matrix.out' matrix w image notitle, \
'' matrix using 1:2:(sprintf('%.2f', $3)) with labels font ',30' notitle
