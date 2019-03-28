
gnuplot

set termoption enhanced

set hidden3d

set colorbox
set pm3d at bs corners2color geomean


set palette model RGB defined (-1.0 "black", -0.2 "violet", -0.1 "blue", 0.0 "green", 0.1 "yellow", 0.2 "orange", 1.0 "red")

set title "Cylindrical Bessel function J_{/Symbol n}(x) - float"
set xlabel "{/Symbol n}"
set ylabel "x"
splot [0:100][0:100][-1:1] "../plot_data/cyl_bessel_j_float.txt" index 0 with pm3d title "J_{/Symbol n}(x)"

set title "Cylindrical Bessel function J_{/Symbol n}(x) - double"
set xlabel "{/Symbol n}"
set ylabel "x"
splot [0:100][0:100][-1:1] "../plot_data/cyl_bessel_j_double.txt" index 0 with pm3d title "J_{/Symbol n}(x)"


set palette model RGB defined (-1.0 "black", -0.2 "violet", -0.1 "blue", 0.0 "green", 0.1 "yellow", 0.2 "orange", 1.0 "red")


set title "Cylindrical Neumann function N_{/Symbol n}(x) - float"
set xlabel "{/Symbol n}"
set ylabel "x"
splot [0:100][0:100][-10:1] "../plot_data/cyl_neumann_float.txt" index 0 with pm3d title "N_{/Symbol n}(x)"

set title "Cylindrical Neumann function N_{/Symbol n}(x) - double"
set xlabel "{/Symbol n}"
set ylabel "x"
splot [0:100][0:100][-10:1] "../plot_data/cyl_neumann_double.txt" index 0 with pm3d title "N_{/Symbol n}(x)"

