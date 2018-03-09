set terminal svg
set output "cos.svg"
set grid

set xrange [0.:2*pi]
set yrange [-1:1]

plot cos(x) lw 2 lt rgb "red" title "cos(x)
