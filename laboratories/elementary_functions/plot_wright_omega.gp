
gnuplot

set termoption enhanced
set xzeroaxis
set yzeroaxis
set grid

set palette model RGB defined (0.0 "blue", 0.3 "green", 0.7 "yellow", 1 "red")
set colorbox
set pm3d at sb corners2color geomean

mag(x,y) = sqrt(x*x+y*y)

set title "Wright omega function {/Symbol w}(z)"
set xlabel "x"
set xlabel "y"
splot [-10.0:10.0][-10.0:10.0][-14.0:10.0] "test_wright_omega.txt" index 0 using 1:2:3 with pm3d title "{/Symbol w}(z)"

set title "Wright omega function {/Symbol w}(z)"
set xlabel "x"
set xlabel "y"
splot [-10.0:10.0][-10.0:10.0][-14.0:10.0] "test_wright_omega.txt" index 0 using 1:2:4 with pm3d title "{/Symbol w}(z)"

set title "Wright omega function {/Symbol w}(z)"
set xlabel "x"
set xlabel "y"
splot [-10.0:10.0][-10.0:10.0][0.0:15.0] "test_wright_omega.txt" index 0 using 1:2:(mag($3,$4)) with pm3d title "{/Symbol w}(z)"


set title "Lambert Wm(x)"
set xlabel "x"
plot [-0.4:0.0][-7.0:-0.8] "test_lambert_wm.txt" using 1:2 with lines title "Wm(x)"

set title "Lambert Wp(x)"
set xlabel "x"
plot [-0.4:4.0][-1.2:1.2] "test_lambert_wp.txt" using 1:2 with lines title "Wp(x)"
