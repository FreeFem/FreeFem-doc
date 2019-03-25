set terminal svg
set output "arccos.svg"
set grid

set xrange [-1:1]
set yrange [0:pi]

plot acos(x) lw 2 lt rgb "red" title "arccos(x)
