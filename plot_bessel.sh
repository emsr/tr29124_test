
/usr/local/bin/gnuplot

set hidden3d

set palette model RGB defined (0.0 "blue", 0.3 "green", 0.7 "yellow", 1 "red")
set palette model RGB defined (-1.0 "black", -0.2 "blue", 0.0 "green", 0.2 "yellow", 1.0 "red")
set colorbox
set pm3d at bs corners2color geomean

splot [0:100][0:100][-1:1] "plot/cyl_bessel_j_float.txt" index 0 with pm3d title "J_\nu(x)"
splot [0:100][0:100][-1:1] "plot/cyl_bessel_j_double.txt" index 0 with pm3d title "J_\nu(x)"

splot [0:100][0:100][-10:1] "plot/cyl_neumann_float.txt" index 0 with pm3d title "N_\nu(x)"
splot [0:100][0:100][-10:1] "plot/cyl_neumann_double.txt" index 0 with pm3d title "N_\nu(x)"
