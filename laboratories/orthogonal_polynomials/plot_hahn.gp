# http://pool-serv1.mathematik.uni-kassel.de/caop/CAOP/plot_hahn.html

gnuplot

load '../plotting_tools/settings.gp'

set title "Hahn Polynomial Q_n(x; {/Symbol a}, {/Symbol b}, N)"
set xlabel "x"
plot [0.0:10.0][-1.4:1.4] \
    "test_hahn.txt" index 5 using 1:2 with lines ls 1 title "Q_5(x; {/Symbol a}, {/Symbol b}, N)"

set title "Hahn Polynomials Q_n(x; {/Symbol a}, {/Symbol b}, N)"
set xlabel "x"
plot [0.0:10.0][-65.0:65.0] \
    "test_hahn.txt" index 0 using 1:2 with lines ls 1 title "Q_0(x; {/Symbol a}, {/Symbol b}, N)", \
                 "" index 1 using 1:2 with lines ls 2 title "Q_1(x; {/Symbol a}, {/Symbol b}, N)", \
                 "" index 2 using 1:2 with lines ls 3 title "Q_2(x; {/Symbol a}, {/Symbol b}, N)", \
                 "" index 3 using 1:2 with lines ls 4 title "Q_3(x; {/Symbol a}, {/Symbol b}, N)", \
                 "" index 4 using 1:2 with lines ls 5 title "Q_4(x; {/Symbol a}, {/Symbol b}, N)", \
                 "" index 5 using 1:2 with lines ls 6 title "Q_5(x; {/Symbol a}, {/Symbol b}, N)", \
                 "" index 6 using 1:2 with lines ls 7 title "Q_6(x; {/Symbol a}, {/Symbol b}, N)", \
                 "" index 7 using 1:2 with lines ls 8 title "Q_7(x; {/Symbol a}, {/Symbol b}, N)", \
                 "" index 8 using 1:2 with lines ls 9 title "Q_8(x; {/Symbol a}, {/Symbol b}, N)", \
                 "" index 9 using 1:2 with lines ls 10 title "Q_9(x; {/Symbol a}, {/Symbol b}, N)", \
                 "" index 10 using 1:2 with lines ls 11 title "Q_{10}(x; {/Symbol a}, {/Symbol b}, N)"
