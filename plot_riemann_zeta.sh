
/usr/local/bin/gnuplot
set hidden3d

splot [-20:5][-5:5][-20:20] "plot/riemann_zeta_float.txt" index 0 with lines title "Re(zeta(s))"
splot [-20:5][-5:5][-20:20] "plot/riemann_zeta_float.txt" index 1 with lines title "Im(zeta(s))"
splot [-20:5][-5:5][0:80] "plot/riemann_zeta_float.txt" index 2 with lines title "|zeta(s)|"
splot [-20:5][-5:5][-180:180] "plot/riemann_zeta_float.txt" index 3 with lines title "Arg(zeta(s))"
