
gnuplot

set hidden3d

set palette model RGB defined (0.0 "blue", 0.3 "green", 0.7 "yellow", 1 "red")
set colorbox
set pm3d at sb corners2color geomean

splot [-50:50][-50:50][-40:120] "plot/beta_double.txt" index 0 with pm3d title "Beta(a,b)"

splot [-20:5][-5:5][-40:20] "plot/log_gamma_lanczos_float.txt" index 0 with pm3d title "Re(logGamma(z))"
splot [-20:5][-5:5][-20:20] "plot/log_gamma_lanczos_float.txt" index 1 with pm3d title "Im(logGamma(z))"
splot [-20:5][-5:5][0:40] "plot/log_gamma_lanczos_float.txt" index 2 with pm3d title "|logGamma(z)|"
splot [-20:5][-5:5][-180:180] "plot/log_gamma_lanczos_float.txt" index 3 with pm3d title "Arg(logGamma(z))"

splot [-20:5][-5:5][-40:20] "plot/log_gamma_lanczos_double.txt" index 0 with pm3d title "Re(logGamma(z))"
splot [-20:5][-5:5][-20:20] "plot/log_gamma_lanczos_double.txt" index 1 with pm3d title "Im(logGamma(z))"
splot [-20:5][-5:5][0:40] "plot/log_gamma_lanczos_double.txt" index 2 with pm3d title "|logGamma(z)|"
splot [-20:5][-5:5][-180:180] "plot/log_gamma_lanczos_double.txt" index 3 with pm3d title "Arg(logGamma(z))"

splot [-20:5][-5:5][-40:20] "plot/log_gamma_lanczos_long_double.txt" index 0 with pm3d title "Re(logGamma(z))"
splot [-20:5][-5:5][-20:20] "plot/log_gamma_lanczos_long_double.txt" index 1 with pm3d title "Im(logGamma(z))"
splot [-20:5][-5:5][0:40] "plot/log_gamma_lanczos_long_double.txt" index 2 with pm3d title "|logGamma(z)|"
splot [-20:5][-5:5][-180:180] "plot/log_gamma_lanczos_long_double.txt" index 3 with pm3d title "Arg(logGamma(z))"

#splot [-20:5][-5:5][-8:20] "plot/log_gamma_lanczos__float128.txt" index 0 with pm3d title "Re(logGamma(z))"
#splot [-20:5][-5:5][-8:20] "plot/log_gamma_lanczos__float128.txt" index 1 with pm3d title "Im(logGamma(z))"
#splot [-20:5][-5:5][0:40] "plot/log_gamma_lanczos__float128.txt" index 2 with pm3d title "|logGamma(z)|"
#splot [-20:5][-5:5][-180:180] "plot/log_gamma_lanczos__float128.txt" index 3 with pm3d title "Arg(logGamma(z))"



splot [-20:5][-5:5][-40:20] "plot/log_gamma_spouge_float.txt" index 0 with pm3d title "Re(logGamma(z))"
splot [-20:5][-5:5][-20:20] "plot/log_gamma_spouge_float.txt" index 1 with pm3d title "Im(logGamma(z))"
splot [-20:5][-5:5][0:40] "plot/log_gamma_spouge_float.txt" index 2 with pm3d title "|logGamma(z)|"
splot [-20:5][-5:5][-180:180] "plot/log_gamma_spouge_float.txt" index 3 with pm3d title "Arg(logGamma(z))"

splot [-20:5][-5:5][-40:20] "plot/log_gamma_spouge_double.txt" index 0 with pm3d title "Re(logGamma(z))"
splot [-20:5][-5:5][-20:20] "plot/log_gamma_spouge_double.txt" index 1 with pm3d title "Im(logGamma(z))"
splot [-20:5][-5:5][0:40] "plot/log_gamma_spouge_double.txt" index 2 with pm3d title "|logGamma(z)|"
splot [-20:5][-5:5][-180:180] "plot/log_gamma_spouge_double.txt" index 3 with pm3d title "Arg(logGamma(z))"

splot [-20:5][-5:5][-40:20] "plot/log_gamma_spouge_long_double.txt" index 0 with pm3d title "Re(logGamma(z))"
splot [-20:5][-5:5][-20:20] "plot/log_gamma_spouge_long_double.txt" index 1 with pm3d title "Im(logGamma(z))"
splot [-20:5][-5:5][0:40] "plot/log_gamma_spouge_long_double.txt" index 2 with pm3d title "|logGamma(z)|"
splot [-20:5][-5:5][-180:180] "plot/log_gamma_spouge_long_double.txt" index 3 with pm3d title "Arg(logGamma(z))"

#splot [-20:5][-5:5][-20:20] "plot/log_gamma_spouge__float128.txt" index 0 with pm3d title "Re(logGamma(z))"
#splot [-20:5][-5:5][-20:20] "plot/log_gamma_spouge__float128.txt" index 1 with pm3d title "Im(logGamma(z))"
#splot [-20:5][-5:5][0:40] "plot/log_gamma_spouge__float128.txt" index 2 with pm3d title "|logGamma(z)|"
#splot [-20:5][-5:5][-180:180] "plot/log_gamma_spouge__float128.txt" index 3 with pm3d title "Arg(logGamma(z))"
