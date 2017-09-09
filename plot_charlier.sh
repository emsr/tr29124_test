# http://pool-serv1.mathematik.uni-kassel.de/caop/CAOP/plot_charlier.html

gnuplot
set termoption enhanced
set xzeroaxis
set yzeroaxis
set grid

set xlabel "x"
plot [0.0:15.0][-20.0:20.0] \
    "test_charlier.txt" index 0 using 1:2 with lines title "C_0(x;a)", \
                     "" index 1 using 1:2 with lines title "C_1(x;a)", \
                     "" index 2 using 1:2 with lines title "C_2(x;a)", \
                     "" index 3 using 1:2 with lines title "C_3(x;a)", \
                     "" index 4 using 1:2 with lines title "C_4(x;a)", \
                     "" index 5 using 1:2 with lines title "C_5(x;a)", \
                     "" index 6 using 1:2 with lines title "C_6(x;a)", \
                     "" index 7 using 1:2 with lines title "C_7(x;a)", \
                     "" index 8 using 1:2 with lines title "C_8(x;a)", \
                     "" index 9 using 1:2 with lines title "C_9(x;a)", \
                     "" index 10 using 1:2 with lines title "C_10(x;a)"
