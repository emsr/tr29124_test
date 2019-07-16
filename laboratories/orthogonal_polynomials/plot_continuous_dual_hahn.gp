# http://pool-serv1.mathematik.uni-kassel.de/caop/CAOP/plot_contdualhahn.html

gnuplot

load '../plotting_tools/settings.gp'

set title "Continuous Dual Hahn Polynomial S_n(x^2; a, b, c)"
set xlabel "x"
plot [0.0:11.0][-5.0:10.0] \
    "test_continuous_dual_hahn.txt" index 5 using 1:2 with lines ls 1 title "S_5(x^2; a, b, c)"

set title "Continuous Dual Hahn Polynomial S_n(x^2; a, b, c)"
set xlabel "x"
plot [0.0:11.0][-10.0:50.0] \
    "test_continuous_dual_hahn.txt" index 0 using 1:2 with lines ls 1 title "S_0(x^2; a, b, c)", \
                   "" index 1 using 1:2 with lines ls 2 title "S_1(x^2; a, b, c)", \
                   "" index 2 using 1:2 with lines ls 3 title "S_2(x^2; a, b, c)", \
                   "" index 3 using 1:2 with lines ls 4 title "S_3(x^2; a, b, c)", \
                   "" index 4 using 1:2 with lines ls 5 title "S_4(x^2; a, b, c)", \
                   "" index 5 using 1:2 with lines ls 6 title "S_5(x^2; a, b, c)", \
                   "" index 6 using 1:2 with lines ls 7 title "S_6(x^2; a, b, c)", \
                   "" index 7 using 1:2 with lines ls 8 title "S_7(x^2; a, b, c)", \
                   "" index 8 using 1:2 with lines ls 9 title "S_8(x^2; a, b, c)", \
                   "" index 9 using 1:2 with lines ls 10 title "S_9(x^2; a, b, c)", \
                   "" index 10 using 1:2 with lines ls 11 title "S_{10}(x^2; a, b, c)"

set term push
set term png
set output "continuous_dual_hahn.png"
replot
set term pop
