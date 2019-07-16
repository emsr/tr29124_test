# http://pool-serv1.mathematik.uni-kassel.de/caop/CAOP/plot_wilson.html

gnuplot

load '../plotting_tools/settings.gp'

set title "Wilson polynomial W_n(x^2; a, b, c, d)"
set xlabel "x"
plot [0.0:20.0][-10.0:20.0] \
    "test_wilson.txt" index 5 using 1:2 with lines ls 1 title "W_5(x^2; a, b, c, d)"

set title "Wilson polynomial W_n(x^2; a, b, c, d)"
set xlabel "x"
plot [0.0:20.0][-125.0:125.0] \
    "test_wilson.txt" index 0 using 1:2 with lines ls 1 title "W_0(x^2; a, b, c, d)", \
                   "" index 1 using 1:2 with lines ls 2 title "W_1(x^2; a, b, c, d)", \
                   "" index 2 using 1:2 with lines ls 3 title "W_2(x^2; a, b, c, d)", \
                   "" index 3 using 1:2 with lines ls 4 title "W_3(x^2; a, b, c, d)", \
                   "" index 4 using 1:2 with lines ls 5 title "W_4(x^2; a, b, c, d)", \
                   "" index 5 using 1:2 with lines ls 6 title "W_5(x^2; a, b, c, d)", \
                   "" index 6 using 1:2 with lines ls 7 title "W_6(x^2; a, b, c, d)", \
                   "" index 7 using 1:2 with lines ls 8 title "W_7(x^2; a, b, c, d)", \
                   "" index 8 using 1:2 with lines ls 9 title "W_8(x^2; a, b, c, d)", \
                   "" index 9 using 1:2 with lines ls 10 title "W_9(x^2; a, b, c, d)", \
                   "" index 10 using 1:2 with lines ls 11 title "W_{10}(x^2; a, b, c, d)"

set term push
set term png
set output "wilson.png"
replot
set term pop
