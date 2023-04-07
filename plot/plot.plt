set title "Final temperature distribution" 
set cblabel "Temperature"
unset key

set tic scale 0
set palette rgb 33, 13, 10

set xrange [0.0:1.0]
set yrange [0.0:1.0]
set size square

set terminal png
set output 'plot/result.png'
plot 'plot/solution.dat' with image
