
/usr/local/bin/gnuplot
set hidden3d

set palette model RGB defined (0.0 "blue", 0.3 "green", 0.7 "yellow", 1 "red")
set colorbox
set pm3d at sb corners2color geomean

#splot [-20:5][-5:5][-20:20] "plot/riemann_zeta_*.txt" index * with lines title "*"

splot [-20:5][-5:5][-20:20] "plot/riemann_zeta_float.txt" index 0 with pm3d title "Re(zeta(s))"
splot [-20:5][-5:5][-20:20] "plot/riemann_zeta_float.txt" index 1 with pm3d title "Im(zeta(s))"
splot [-20:5][-5:5][0:80] "plot/riemann_zeta_float.txt" index 2 with pm3d title "|zeta(s)|"
splot [-20:5][-5:5][-180:180] "plot/riemann_zeta_float.txt" index 3 with pm3d title "Arg(zeta(s))"

splot [-20:5][-5:5][-20:20] "plot/riemann_zeta_double.txt" index 0 with pm3d title "Re(zeta(s))"
splot [-20:5][-5:5][-20:20] "plot/riemann_zeta_double.txt" index 1 with pm3d title "Im(zeta(s))"
splot [-20:5][-5:5][0:80] "plot/riemann_zeta_double.txt" index 2 with pm3d title "|zeta(s)|"
splot [-20:5][-5:5][-180:180] "plot/riemann_zeta_double.txt" index 3 with pm3d title "Arg(zeta(s))"

splot [-20:5][-5:5][-20:20] "plot/riemann_zeta_long_double.txt" index 0 with pm3d title "Re(zeta(s))"
splot [-20:5][-5:5][-20:20] "plot/riemann_zeta_long_double.txt" index 1 with pm3d title "Im(zeta(s))"
splot [-20:5][-5:5][0:80] "plot/riemann_zeta_long_double.txt" index 2 with pm3d title "|zeta(s)|"
splot [-20:5][-5:5][-180:180] "plot/riemann_zeta_long_double.txt" index 3 with pm3d title "Arg(zeta(s))"
