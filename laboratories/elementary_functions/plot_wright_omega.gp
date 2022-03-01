
gnuplot

set termoption enhanced
set xzeroaxis
set yzeroaxis
set grid

set palette model RGB defined (0.0 "blue", 0.3 "green", 0.7 "yellow", 1 "red")
set colorbox
set pm3d at sb corners2color geomean

mag(x,y) = sqrt(x*x+y*y)

set title "Wright omega function W_0(z)"
set xlabel "p"
splot [-10.0:10.0][-10.0:10.0][-14.0:10.0] "test_wright_omega.txt" index 0 using 1:2:3 with pm3d title "W_0(z)"

set title "Wright omega function W_0(z)"
set xlabel "p"
splot [-10.0:10.0][-10.0:10.0][-14.0:10.0] "test_wright_omega.txt" index 0 using 1:2:4 with pm3d title "W_0(z)"

set title "Wright omega function W_0(z)"
set xlabel "p"
splot [-10.0:10.0][-10.0:10.0][0.0:15.0] "test_wright_omega.txt" index 0 using 1:2:(mag($3,$4)) with pm3d title "W_0(z)"
