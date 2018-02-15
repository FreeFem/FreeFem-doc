set terminal svg
set output "arcsin.svg"
set grid

set xrange [-1:1]
set yrange [-pi/2.:pi/2.]

plot asin(x) lw 2 lt rgb "red" title "arcsin(x)
