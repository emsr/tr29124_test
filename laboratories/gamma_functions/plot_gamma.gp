
gnuplot

set termoption enhanced

set hidden3d

set palette model RGB defined (0.0 "blue", 0.3 "green", 0.7 "yellow", 1 "red")
set colorbox
set pm3d at bs corners2color geomean

set title "Beta function B(a, b)"
splot [-3:3][-3:3][-40:120] "plot_data/beta_double.txt" index 0 with pm3d title "B(a,b)"


set title "log Gamma function Re{ln[{/Symbol G}(a)} - Lanczos method"
splot [-20:5][-5:5][-40:20] "plot_data/log_gamma_lanczos_float.txt" index 0 with pm3d title "Re(log{/Symbol G}(z))"

set title "log Gamma function Im{ln[{/Symbol G}(a)} - Lanczos method"
splot [-20:5][-5:5][-20:20] "plot_data/log_gamma_lanczos_float.txt" index 1 with pm3d title "Im(log{/Symbol G}(z))"

set title "log Gamma function |ln[{/Symbol G}(a)| - Lanczos method"
splot [-20:5][-5:5][0:40] "plot_data/log_gamma_lanczos_float.txt" index 2 with pm3d title "|log{/Symbol G}(z)|"

set title "log Gamma function Arg(ln[{/Symbol G}(a)) - Lanczos method"
splot [-20:5][-5:5][-180:180] "plot_data/log_gamma_lanczos_float.txt" index 3 with pm3d title "Arg(log{/Symbol G}(z))"


set title "log Gamma function Re{ln[{/Symbol G}(a)} - Lanczos method"
splot [-20:5][-5:5][-40:20] "plot_data/log_gamma_lanczos_double.txt" index 0 with pm3d title "Re(log{/Symbol G}(z))"

set title "log Gamma function Im{ln[{/Symbol G}(a)} - Lanczos method"
splot [-20:5][-5:5][-20:20] "plot_data/log_gamma_lanczos_double.txt" index 1 with pm3d title "Im(log{/Symbol G}(z))"

set title "log Gamma function |ln[{/Symbol G}(a)| - Lanczos method"
splot [-20:5][-5:5][0:40] "plot_data/log_gamma_lanczos_double.txt" index 2 with pm3d title "|log{/Symbol G}(z)|"

set title "log Gamma function Arg(ln[{/Symbol G}(a)) - Lanczos method"
splot [-20:5][-5:5][-180:180] "plot_data/log_gamma_lanczos_double.txt" index 3 with pm3d title "Arg(log{/Symbol G}(z))"


set title "log Gamma function Re{ln[{/Symbol G}(a)} - Lanczos method"
splot [-20:5][-5:5][-40:20] "plot_data/log_gamma_lanczos_long_double.txt" index 0 with pm3d title "Re(log{/Symbol G}(z))"

set title "log Gamma function Im{ln[{/Symbol G}(a)} - Lanczos method"
splot [-20:5][-5:5][-20:20] "plot_data/log_gamma_lanczos_long_double.txt" index 1 with pm3d title "Im(log{/Symbol G}(z))"

set title "log Gamma function |ln[{/Symbol G}(a)| - Lanczos method"
splot [-20:5][-5:5][0:40] "plot_data/log_gamma_lanczos_long_double.txt" index 2 with pm3d title "|log{/Symbol G}(z)|"

set title "log Gamma function Arg(ln[{/Symbol G}(a)) - Lanczos method"
splot [-20:5][-5:5][-180:180] "plot_data/log_gamma_lanczos_long_double.txt" index 3 with pm3d title "Arg(log{/Symbol G}(z))"


set title "log Gamma function Re{ln[{/Symbol G}(a)} - Lanczos method"
splot [-20:5][-5:5][-8:20] "plot_data/log_gamma_lanczos__float128.txt" index 0 with pm3d title "Re(log{/Symbol G}(z))"

set title "log Gamma function Im{ln[{/Symbol G}(a)} - Lanczos method"
splot [-20:5][-5:5][-8:20] "plot_data/log_gamma_lanczos__float128.txt" index 1 with pm3d title "Im(log{/Symbol G}(z))"

set title "log Gamma function |ln[{/Symbol G}(a)| - Lanczos method"
splot [-20:5][-5:5][0:40] "plot_data/log_gamma_lanczos__float128.txt" index 2 with pm3d title "|log{/Symbol G}(z)|"

set title "log Gamma function Arg(ln[{/Symbol G}(a)) - Lanczos method"
splot [-20:5][-5:5][-180:180] "plot_data/log_gamma_lanczos__float128.txt" index 3 with pm3d title "Arg(log{/Symbol G}(z))"



set title "log Gamma function Re{ln[{/Symbol G}(a)} - Spouge method"
splot [-20:5][-5:5][-40:20] "plot_data/log_gamma_spouge_float.txt" index 0 with pm3d title "Re(log{/Symbol G}(z))"

set title "log Gamma function Im{ln[{/Symbol G}(a)} - Lanczos method"
splot [-20:5][-5:5][-20:20] "plot_data/log_gamma_spouge_float.txt" index 1 with pm3d title "Im(log{/Symbol G}(z))"

set title "log Gamma function |ln[{/Symbol G}(a)| - Lanczos method"
splot [-20:5][-5:5][0:40] "plot_data/log_gamma_spouge_float.txt" index 2 with pm3d title "|log{/Symbol G}(z)|"

splot [-20:5][-5:5][-180:180] "plot_data/log_gamma_spouge_float.txt" index 3 with pm3d title "Arg(log{/Symbol G}(z))"


set title "log Gamma function Re{ln[{/Symbol G}(a)} - Spouge method"
splot [-20:5][-5:5][-40:20] "plot_data/log_gamma_spouge_double.txt" index 0 with pm3d title "Re(log{/Symbol G}(z))"

set title "log Gamma function Im{ln[{/Symbol G}(a)} - Lanczos method"
splot [-20:5][-5:5][-20:20] "plot_data/log_gamma_spouge_double.txt" index 1 with pm3d title "Im(log{/Symbol G}(z))"

set title "log Gamma function |ln[{/Symbol G}(a)| - Lanczos method"
splot [-20:5][-5:5][0:40] "plot_data/log_gamma_spouge_double.txt" index 2 with pm3d title "|log{/Symbol G}(z)|"

set title "log Gamma function Arg(ln[{/Symbol G}(a)) - Lanczos method"
splot [-20:5][-5:5][-180:180] "plot_data/log_gamma_spouge_double.txt" index 3 with pm3d title "Arg(log{/Symbol G}(z))"


set title "log Gamma function Re{ln[{/Symbol G}(a)} - Spouge method"
splot [-20:5][-5:5][-40:20] "plot_data/log_gamma_spouge_long_double.txt" index 0 with pm3d title "Re(log{/Symbol G}(z))"

set title "log Gamma function Im{ln[{/Symbol G}(a)} - Lanczos method"
splot [-20:5][-5:5][-20:20] "plot_data/log_gamma_spouge_long_double.txt" index 1 with pm3d title "Im(log{/Symbol G}(z))"

set title "log Gamma function |ln[{/Symbol G}(a)| - Lanczos method"
splot [-20:5][-5:5][0:40] "plot_data/log_gamma_spouge_long_double.txt" index 2 with pm3d title "|log{/Symbol G}(z)|"

set title "log Gamma function Arg(ln[{/Symbol G}(a)) - Lanczos method"
splot [-20:5][-5:5][-180:180] "plot_data/log_gamma_spouge_long_double.txt" index 3 with pm3d title "Arg(log{/Symbol G}(z))"


set title "log Gamma function Re{ln[{/Symbol G}(a)} - Spouge method"
splot [-20:5][-5:5][-20:20] "plot_data/log_gamma_spouge__float128.txt" index 0 with pm3d title "Re(log{/Symbol G}(z))"

set title "log Gamma function Im{ln[{/Symbol G}(a)} - Lanczos method"
splot [-20:5][-5:5][-20:20] "plot_data/log_gamma_spouge__float128.txt" index 1 with pm3d title "Im(log{/Symbol G}(z))"

set title "log Gamma function |ln[{/Symbol G}(a)| - Lanczos method"
splot [-20:5][-5:5][0:40] "plot_data/log_gamma_spouge__float128.txt" index 2 with pm3d title "|log{/Symbol G}(z)|"

set title "log Gamma function Arg(ln[{/Symbol G}(a)) - Lanczos method"
splot [-20:5][-5:5][-180:180] "plot_data/log_gamma_spouge__float128.txt" index 3 with pm3d title "Arg(log{/Symbol G}(z))"
