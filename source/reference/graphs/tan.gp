set terminal svg
set output "tan.svg"
set grid

set xrange [-pi/2.:pi/2.]
set yrange [-pi:pi]

plot tan(x) lw 2 lt rgb "red" title "tan(x)
