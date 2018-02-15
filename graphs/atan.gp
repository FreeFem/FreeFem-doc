set terminal svg
set output "arctan.svg"
set grid

set xrange [-10:10]
set yrange [-pi/2.:pi/2.]

plot atan(x) lw 2 lt rgb "red" title "arctan(x)
