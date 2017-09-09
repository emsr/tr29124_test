# http://pool-serv1.mathematik.uni-kassel.de/caop/CAOP/plot_hahn.html

gnuplot
set termoption enhanced
set xzeroaxis
set yzeroaxis
set grid

set xlabel "x"
plot [0.0:15.0][-20.0:20.0] \
    "test_hahn.txt" index 0 using 1:2 with lines title "Q_0(x; {/Symbol a}, {/Symbol b}, N)", \
                 "" index 1 using 1:2 with lines title "Q_1(x; {/Symbol a}, {/Symbol b}, N)", \
                 "" index 2 using 1:2 with lines title "Q_2(x; {/Symbol a}, {/Symbol b}, N)", \
                 "" index 3 using 1:2 with lines title "Q_3(x; {/Symbol a}, {/Symbol b}, N)", \
                 "" index 4 using 1:2 with lines title "Q_4(x; {/Symbol a}, {/Symbol b}, N)", \
                 "" index 5 using 1:2 with lines title "Q_5(x; {/Symbol a}, {/Symbol b}, N)", \
                 "" index 6 using 1:2 with lines title "Q_6(x; {/Symbol a}, {/Symbol b}, N)", \
                 "" index 7 using 1:2 with lines title "Q_7(x; {/Symbol a}, {/Symbol b}, N)", \
                 "" index 8 using 1:2 with lines title "Q_8(x; {/Symbol a}, {/Symbol b}, N)", \
                 "" index 9 using 1:2 with lines title "Q_9(x; {/Symbol a}, {/Symbol b}, N)", \
                 "" index 10 using 1:2 with lines title "Q_10(x; {/Symbol a}, {/Symbol b}, N)"
