set terminal svg
set output "cosh.svg"
set grid

set xrange [-pi:pi]
set yrange [0:10]

plot cosh(x) lw 2 lt rgb "red" title "cosh(x)
