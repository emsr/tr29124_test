# http://pool-serv1.mathematik.uni-kassel.de/caop/CAOP/plot_conthahn.html

gnuplot

load 'settings.gp'

set xlabel "x"
plot [-2.0:2.0][-6500.0:6500.0] \
    "test_continuous_hahn.txt" index 5 using 1:2 with lines ls 1 title "p_5(x; a, b, c, d)"

set xlabel "x"
plot [-2.0:2.0][-75000.0:75000.0] \
    "test_continuous_hahn.txt" index 0 using 1:2 with lines ls 1 title "p_0(x; a, b, c, d)", \
                   "" index 1 using 1:2 with lines ls 2 title "p_1(x; a, b, c, d)", \
                   "" index 2 using 1:2 with lines ls 3 title "p_2(x; a, b, c, d)", \
                   "" index 3 using 1:2 with lines ls 4 title "p_3(x; a, b, c, d)", \
                   "" index 4 using 1:2 with lines ls 5 title "p_4(x; a, b, c, d)", \
                   "" index 5 using 1:2 with lines ls 6 title "p_5(x; a, b, c, d)", \
                   "" index 6 using 1:2 with lines ls 7 title "p_6(x; a, b, c, d)", \
                   "" index 7 using 1:2 with lines ls 8 title "p_7(x; a, b, c, d)", \
                   "" index 8 using 1:2 with lines ls 9 title "p_8(x; a, b, c, d)", \
                   "" index 9 using 1:2 with lines ls 10 title "p_9(x; a, b, c, d)", \
                   "" index 10 using 1:2 with lines ls 11 title "p_{10}(x; a, b, c, d)"
