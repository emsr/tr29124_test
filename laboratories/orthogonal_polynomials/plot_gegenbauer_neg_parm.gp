
gnuplot

load 'settings.gp'

set title "gegenbauer polynomial P_n^{({/Symbol l})}(x)"
set xlabel "x"
plot [-2.0:2.0][-100.0:100.0] \
    "test_gegenbauer_neg_parm.txt" index 0 using 1:2 with lines ls 1 title "C_n^{({/Symbol l}-)}(x)", \
                                "" index 0 using 1:3 with lines ls 2 title "C_n^{({/Symbol l})}(x)", \
                                "" index 0 using 1:4 with lines ls 3 title "CW_n^{({/Symbol l}+)}(x)"

set title "gegenbauer polynomial P_n^{({/Symbol l})}(x)"
set xlabel "x"
plot [-2.0:2.0][-100.0:100.0] \
    "test_gegenbauer_neg_parm.txt" index 1 using 1:2 with lines ls 1 title "C_n^{({/Symbol l}-)}(x)", \
                                "" index 1 using 1:3 with lines ls 2 title "C_n^{({/Symbol l})}(x)", \
                                "" index 1 using 1:4 with lines ls 3 title "CW_n^{({/Symbol l}+)}(x)"

set title "gegenbauer polynomial P_n^{({/Symbol l})}(x)"
set xlabel "x"
plot [-2.0:2.0][-100.0:100.0] \
    "test_gegenbauer_neg_parm.txt" index 2 using 1:2 with lines ls 1 title "C_n^{({/Symbol l}-)}(x)", \
                                "" index 2 using 1:3 with lines ls 2 title "C_n^{({/Symbol l})}(x)", \
                                "" index 2 using 1:4 with lines ls 3 title "CW_n^{({/Symbol l}+)}(x)"
