
/usr/local/bin/gnuplot
set hidden3d

splot [-20:5][-5:5][-20:20] "plot/gamma_lanczos_float.txt" index 0 with lines title "Re(logGamma(z))"
splot [-20:5][-5:5][-20:20] "plot/gamma_lanczos_float.txt" index 1 with lines title "Im(logGamma(z))"
splot [-20:5][-5:5][0:40] "plot/gamma_lanczos_float.txt" index 2 with lines title "|logGamma(z)|"
splot [-20:5][-5:5][-180:180] "plot/gamma_lanczos_float.txt" index 3 with lines title "Arg(logGamma(z))"

splot [-20:5][-5:5][-20:20] "plot/gamma_lanczos_double.txt" index 0 with lines title "Re(logGamma(z))"
splot [-20:5][-5:5][-20:20] "plot/gamma_lanczos_double.txt" index 1 with lines title "Im(logGamma(z))"
splot [-20:5][-5:5][0:40] "plot/gamma_lanczos_double.txt" index 2 with lines title "|logGamma(z)|"
splot [-20:5][-5:5][-180:180] "plot/gamma_lanczos_double.txt" index 3 with lines title "Arg(logGamma(z))"

splot [-20:5][-5:5][-8:20] "plot/gamma_lanczos_long_double.txt" index 0 with lines title "Re(logGamma(z))"
splot [-20:5][-5:5][-8:20] "plot/gamma_lanczos_long_double.txt" index 1 with lines title "Im(logGamma(z))"
splot [-20:5][-5:5][0:40] "plot/gamma_lanczos_long_double.txt" index 2 with lines title "|logGamma(z)|"
splot [-20:5][-5:5][-180:180] "plot/gamma_lanczos_long_double.txt" index 3 with lines title "Arg(logGamma(z))"

#splot [-20:5][-5:5][-8:20] "plot/gamma_lanczos__float128.txt" index 0 with lines title "Re(logGamma(z))"
#splot [-20:5][-5:5][-8:20] "plot/gamma_lanczos__float128.txt" index 1 with lines title "Im(logGamma(z))"
#splot [-20:5][-5:5][0:40] "plot/gamma_lanczos__float128.txt" index 2 with lines title "|logGamma(z)|"
#splot [-20:5][-5:5][-180:180] "plot/gamma_lanczos__float128.txt" index 3 with lines title "Arg(logGamma(z))"



splot [-20:5][-5:5][-20:20] "plot/gamma_spouge_float.txt" index 0 with lines title "Re(logGamma(z))"
splot [-20:5][-5:5][-20:20] "plot/gamma_spouge_float.txt" index 1 with lines title "Im(logGamma(z))"
splot [-20:5][-5:5][0:40] "plot/gamma_spouge_float.txt" index 2 with lines title "|logGamma(z)|"
splot [-20:5][-5:5][-180:180] "plot/gamma_spouge_float.txt" index 3 with lines title "Arg(logGamma(z))"

splot [-20:5][-5:5][-20:20] "plot/gamma_spouge_double.txt" index 0 with lines title "Re(logGamma(z))"
splot [-20:5][-5:5][-20:20] "plot/gamma_spouge_double.txt" index 1 with lines title "Im(logGamma(z))"
splot [-20:5][-5:5][0:40] "plot/gamma_spouge_double.txt" index 2 with lines title "|logGamma(z)|"
splot [-20:5][-5:5][-180:180] "plot/gamma_spouge_double.txt" index 3 with lines title "Arg(logGamma(z))"

splot [-20:5][-5:5][-20:20] "plot/gamma_spouge_long_double.txt" index 0 with lines title "Re(logGamma(z))"
splot [-20:5][-5:5][-20:20] "plot/gamma_spouge_long_double.txt" index 1 with lines title "Im(logGamma(z))"
splot [-20:5][-5:5][0:40] "plot/gamma_spouge_long_double.txt" index 2 with lines title "|logGamma(z)|"
splot [-20:5][-5:5][-180:180] "plot/gamma_spouge_long_double.txt" index 3 with lines title "Arg(logGamma(z))"

#splot [-20:5][-5:5][-20:20] "plot/gamma_spouge__float128.txt" index 0 with lines title "Re(logGamma(z))"
#splot [-20:5][-5:5][-20:20] "plot/gamma_spouge__float128.txt" index 1 with lines title "Im(logGamma(z))"
#splot [-20:5][-5:5][0:40] "plot/gamma_spouge__float128.txt" index 2 with lines title "|logGamma(z)|"
#splot [-20:5][-5:5][-180:180] "plot/gamma_spouge__float128.txt" index 3 with lines title "Arg(logGamma(z))"
