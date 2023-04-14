set terminal png
set title "Temperature distribution over space" 
set cblabel "Temperature"
unset key
set tic scale 0
set palette rgb 33, 13, 10
set xrange [0.0:1.0]
set yrange [0.0:1.0]
set size square

set output 'plot/result.png'
plot 'data/solution.dat' with image
