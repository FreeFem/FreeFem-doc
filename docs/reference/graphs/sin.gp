set terminal svg
set output "sin.svg"
set grid

set xrange [0:2.*pi]
set yrange [-1.:1.]

plot sin(x) lw 2 lt rgb "red" title "sin(x)
