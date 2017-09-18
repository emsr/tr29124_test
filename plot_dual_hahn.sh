# http://pool-serv1.mathematik.uni-kassel.de/caop/CAOP/plot_dualhahn.html

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

set xlabel "{/Symbol l}(x)"
plot [0.0:100.0][-2.0:12.0] \
    "test_dual_hahn.txt" index 5 using 1:2 with lines ls 1 title "R_5({/Symbol l}(x); {/Symbol g}, {/Symbol d}, N)"

set xlabel "{/Symbol l}(x)"
plot [0.0:100.0][-20.0:200.0] \
    "test_dual_hahn.txt" index 0 using 1:2 with lines ls 1 title "R_0({/Symbol l}(x); {/Symbol g}, {/Symbol d}, N)", \
                 "" index 1 using 1:2 with lines ls 2 title "R_1({/Symbol l}(x); {/Symbol g}, {/Symbol d}, N)", \
                 "" index 2 using 1:2 with lines ls 3 title "R_2({/Symbol l}(x); {/Symbol g}, {/Symbol d}, N)", \
                 "" index 3 using 1:2 with lines ls 4 title "R_3({/Symbol l}(x); {/Symbol g}, {/Symbol d}, N)", \
                 "" index 4 using 1:2 with lines ls 5 title "R_4({/Symbol l}(x); {/Symbol g}, {/Symbol d}, N)", \
                 "" index 5 using 1:2 with lines ls 6 title "R_5({/Symbol l}(x); {/Symbol g}, {/Symbol d}, N)", \
                 "" index 6 using 1:2 with lines ls 7 title "R_6({/Symbol l}(x); {/Symbol g}, {/Symbol d}, N)", \
                 "" index 7 using 1:2 with lines ls 8 title "R_7({/Symbol l}(x); {/Symbol g}, {/Symbol d}, N)", \
                 "" index 8 using 1:2 with lines ls 9 title "R_8({/Symbol l}(x); {/Symbol g}, {/Symbol d}, N)", \
                 "" index 9 using 1:2 with lines ls 10 title "R_9({/Symbol l}(x); {/Symbol g}, {/Symbol d}, N)", \
                 "" index 10 using 1:2 with lines ls 11 title "R_{10}({/Symbol l}(x); {/Symbol g}, {/Symbol d}, N)"
