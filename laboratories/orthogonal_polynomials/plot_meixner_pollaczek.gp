# http://pool-serv1.mathematik.uni-kassel.de/caop/CAOP/plot_meixpol.html

gnuplot

load '../plotting_tools/settings.gp'

set title "Meixner-Pollaczek Polynomial P^{({/Symbol l})}_0(x; {/Symbol f})"
set xlabel "x"
plot [-3.2:3.2][-12.0:12.0] \
    "test_meixner_pollaczek.txt" index 5 using 1:2 with lines ls 1 title "P^{({/Symbol l})}_5(x; {/Symbol f})"

set title "Meixner-Pollaczek Polynomial P^{({/Symbol l})}_0(x; {/Symbol f})"
set xlabel "x"
plot [-3.2:3.2][-12.0:12.0] \
    "test_meixner_pollaczek.txt" index 0 using 1:2 with lines ls 1 title "P^{({/Symbol l})}_0(x; {/Symbol f})", \
                   "" index 1 using 1:2 with lines ls 2 title "P^{({/Symbol l})}_1(x; {/Symbol f})", \
                   "" index 2 using 1:2 with lines ls 3 title "P^{({/Symbol l})}_2(x; {/Symbol f})", \
                   "" index 3 using 1:2 with lines ls 4 title "P^{({/Symbol l})}_3(x; {/Symbol f})", \
                   "" index 4 using 1:2 with lines ls 5 title "P^{({/Symbol l})}_4(x; {/Symbol f})", \
                   "" index 5 using 1:2 with lines ls 6 title "P^{({/Symbol l})}_5(x; {/Symbol f})", \
                   "" index 6 using 1:2 with lines ls 7 title "P^{({/Symbol l})}_6(x; {/Symbol f})", \
                   "" index 7 using 1:2 with lines ls 8 title "P^{({/Symbol l})}_7(x; {/Symbol f})", \
                   "" index 8 using 1:2 with lines ls 9 title "P^{({/Symbol l})}_8(x; {/Symbol f})", \
                   "" index 9 using 1:2 with lines ls 10 title "P^{({/Symbol l})}_9(x; {/Symbol f})", \
                   "" index 10 using 1:2 with lines ls 11 title "P^{({/Symbol l})}_{10}(x; {/Symbol f})"
