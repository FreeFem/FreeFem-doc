set terminal png
set output "OptimalControl_J.png"

set grid
set yrange [0: 45]

set style line 2 lc rgb 'red' pt 7
plot "J.txt" title "J" ls 2

