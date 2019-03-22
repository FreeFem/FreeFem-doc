set terminal svg
set output "arccosh.svg"
set grid

set xrange [1:100]

plot acosh(x) lw 2 lt rgb "red" title "arccosh(x)
