set terminal svg
set output "arcsinh.svg"
set grid

set xrange [-100:100]

plot asinh(x) lw 2 lt rgb "red" title "arcsinh(x)
