# http://pool-serv1.mathematik.uni-kassel.de/caop/CAOP/plot_hahn.html

gnuplot
set termoption enhanced
set xzeroaxis
set yzeroaxis
set grid

set style line  1 lw 1.2 lc  1
set style line  2 lw 1.2 lc  2
set style line  3 lw 1.2 lc  3
set style line  4 lw 1.2 lc  4
set style line  5 lw 1.2 lc  5
set style line  6 lw 1.2 lc  6
set style line  7 lw 1.2 lc  7
set style line  8 lw 1.2 lc  8
set style line  9 lw 1.2 lc  9 dt 2
set style line 10 lw 1.2 lc 10 dt 2
set style line 11 lw 1.2 lc 11 dt 2
set style line 12 lw 1.2 lc 12 dt 2

set xlabel "x"
plot [0.0:10.0][-1.4:1.4] \
    "test_hahn.txt" index 5 using 1:2 with lines ls 1 title "Q_5(x; {/Symbol a}, {/Symbol b}, N)"

set xlabel "x"
plot [0.0:10.0][-20.0:20.0] \
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
