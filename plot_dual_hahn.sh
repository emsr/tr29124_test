# http://pool-serv1.mathematik.uni-kassel.de/caop/CAOP/plot_dualhahn.html

gnuplot

load 'settings.gp'

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
