set terminal svg
set output "sinh.svg"
set grid

set xrange [-pi:pi]
set yrange [-2.*pi:2.*pi]

plot sinh(x) lw 2 lt rgb "red" title "sinh(x)
