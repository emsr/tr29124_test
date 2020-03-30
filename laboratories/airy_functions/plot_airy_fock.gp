gnuplot
set termoption enhanced
set xzeroaxis
set yzeroaxis
set grid

set palette model RGB defined (-1 "black", -0.50 "red", -0.20 "yellow", -0.1 "green", 0.0 "blue", 0.1 "green", 0.20 "yellow", 0.50 "red", 1 "black")
set colorbox
set pm3d clip4in at bs corners2color geomean

set xlabel "Re[z]"
set ylabel "Im[z]"

splot [-20:20][-20:20] "test_airy_fock.txt" using 1:2:(log(sqrt($3*$3+$4*$4))) with pm3d title "log[|w_1(z)|]"

splot [-20:20][-20:20] "test_airy_fock.txt" using 1:2:(log(sqrt($5*$5+$6*$6))) with pm3d title "log[|w_2(z)|]"

splot [-20:20][-20:20] "test_airy_fock.txt" using 1:2:(log(sqrt($7*$7+$8*$8))) with pm3d title "log[|w'_1(z)|]"

splot [-20:20][-20:20] "test_airy_fock.txt" using 1:2:(log(sqrt($9*$9+$10*$10))) with pm3d title "log[|w'_2(z)|]"
