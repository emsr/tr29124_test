gnuplot


set title "Appell function F_1(x,y)"
set xlabel "x"
set ylabel "y"
splot [-1:1][-1:1][-1000:1000] "test_appell.txt" index 0 with pm3d title "F_1(x,y)"
