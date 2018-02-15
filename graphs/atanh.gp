set terminal svg
set output "arctanh.svg"
set grid

set xrange [-1:1]

plot atanh(x) lw 2 lt rgb "red" title "arctanh(x)
