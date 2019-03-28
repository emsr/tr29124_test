gnuplot

load '../plotting_tools/settings.gp'

set title "Legendre function P_l^{(m)}(x)"
set xlabel "x"
plot [-1.0:1.0][-3.0:3.0] \
    "test_assoc_legendre.txt" index 0 using 1:2 with lines title "P_0(x)", \
                           "" index 1 using 1:2 with lines title "P_1(x)", \
                           "" index 2 using 1:2 with lines title "P_2(x)", \
                           "" index 3 using 1:2 with lines title "P_3(x)", \
                           "" index 4 using 1:2 with lines title "P_4(x)", \
                           "" index 5 using 1:2 with lines title "P_5(x)"
