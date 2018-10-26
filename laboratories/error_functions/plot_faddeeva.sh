gnuplot
set termoption enhanced
set xzeroaxis
set yzeroaxis
set grid

h = 10

set cbrange [-h:+h]
set palette model RGB defined (-1 "black", -0.20 "red", -0.10 "yellow", -0.05 "green", 0.0 "blue", 0.05 "green", 0.10 "yellow", 0.20 "red", 1 "black")
set colorbox
set pm3d clip4in at bs corners2color geomean

clamp(z, lo, hi) = (z < lo) ? lo : ((z > hi) ? hi : z)

set xlabel "Re[z]"
set ylabel "Im[z]"

splot [-5.0:5.0][-5.0:5.0][-h:+h] \
  "test_faddeeva.txt" using 1:2:(clamp($3,-h,+h)) with pm3d title "Re[w(z)]"

splot [-5.0:5.0][-5.0:5.0][-h:+h] \
  "test_faddeeva.txt" using 1:2:(clamp($4,-h,+h)) with pm3d title "Im[w(z)]"
