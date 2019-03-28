
gnuplot

set termoption enhanced

set xzeroaxis
set yzeroaxis
set grid

set title "{/Symbol z}(s)"
set xlabel "s"
plot [-25:25][-20:20] \
    "test_riemann_zeta.txt" index 0 using 1:2 with lines title "{/Symbol z}(s)", \
                         "" index 0 using 1:3 with lines title "{/Symbol z}(s) GSL"

set title "{/Symbol z}(1/2 +i{/Symbol t})"
set xlabel "{/Symbol t}"
plot [-25:25][-2.5:2.5] \
    "test_riemann_zeta.txt" index 2 using 1:2 with lines title "Re[{/Symbol z}(1/2 +i{/Symbol t})]", \
                         "" index 2 using 1:3 with lines title "Im[{/Symbol z}(1/2 +i{/Symbol t})]"

gnuplot

set termoption enhanced

set hidden3d

set palette model RGB defined (0.0 "blue", 0.3 "green", 0.7 "yellow", 1 "red")
set colorbox
set pm3d at bs corners2color geomean

splot [-20:5][-5:5][-20:20] "../plot_data/riemann_zeta_float.txt" index 0 with pm3d title "Re({/Symbol z}(s))"
splot [-20:5][-5:5][-20:20] "../plot_data/riemann_zeta_float.txt" index 1 with pm3d title "Im({/Symbol z}(s))"
splot [-20:5][-5:5][0:80] "../plot_data/riemann_zeta_float.txt" index 2 with pm3d title "|{/Symbol z}(s)|"
splot [-20:5][-5:5][-180:180] "../plot_data/riemann_zeta_float.txt" index 3 with pm3d title "Arg({/Symbol z}(s))"


splot [-20:5][-5:5][-20:20] "../plot_data/riemann_zeta_double.txt" index 0 with pm3d title "Re({/Symbol z}(s))"
splot [-20:5][-5:5][-20:20] "../plot_data/riemann_zeta_double.txt" index 1 with pm3d title "Im({/Symbol z}(s))"
splot [-20:5][-5:5][0:80] "../plot_data/riemann_zeta_double.txt" index 2 with pm3d title "|{/Symbol z}(s)|"
splot [-20:5][-5:5][-180:180] "../plot_data/riemann_zeta_double.txt" index 3 with pm3d title "Arg({/Symbol z}(s))"


splot [-20:5][-5:5][-20:20] "../plot_data/riemann_zeta_long_double.txt" index 0 with pm3d title "Re({/Symbol z}(s))"
splot [-20:5][-5:5][-20:20] "../plot_data/riemann_zeta_long_double.txt" index 1 with pm3d title "Im({/Symbol z}(s))"
splot [-20:5][-5:5][0:80] "../plot_data/riemann_zeta_long_double.txt" index 2 with pm3d title "|{/Symbol z}(s)|"
splot [-20:5][-5:5][-180:180] "../plot_data/riemann_zeta_long_double.txt" index 3 with pm3d title "Arg({/Symbol z}(s))"
