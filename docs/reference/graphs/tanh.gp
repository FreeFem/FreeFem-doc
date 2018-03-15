set terminal svg
set output "tanh.svg"
set grid

set xrange [-pi:pi]
set yrange [-1:1]

plot tanh(x) lw 2 lt rgb "red" title "tanh(x)
