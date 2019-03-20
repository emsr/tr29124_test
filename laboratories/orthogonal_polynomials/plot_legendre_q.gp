gnuplot

load 'setings.gp'

set title "Legendre function Q_l(x)"
set xlabel "x"
plot [-1.2:1.2][-3.0:3.0] \
    "test_legendre_q.txt" index 0 using 1:2 with lines title "Q_0(x)", \
                       "" index 1 using 1:2 with lines title "Q_1(x)", \
                       "" index 2 using 1:2 with lines title "Q_2(x)", \
                       "" index 3 using 1:2 with lines title "Q_3(x)", \
                       "" index 4 using 1:2 with lines title "Q_4(x)", \
                       "" index 5 using 1:2 with lines title "Q_5(x)"
