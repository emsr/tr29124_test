# http://pool-serv1.mathematik.uni-kassel.de/caop/CAOP/plot_krawtchouk.html

gnuplot
set termoption enhanced
set xzeroaxis
set yzeroaxis
set grid

set xlabel "x"
plot [0.0:15.0][-20.0:20.0] \
    "test_krawtchouk.txt" index 0 using 1:2 with lines title "K_0(x; p, N)", \
                     "" index 1 using 1:2 with lines title "K_1(x; p, N)", \
                     "" index 2 using 1:2 with lines title "K_2(x; p, N)", \
                     "" index 3 using 1:2 with lines title "K_3(x; p, N)", \
                     "" index 4 using 1:2 with lines title "K_4(x; p, N)", \
                     "" index 5 using 1:2 with lines title "K_5(x; p, N)", \
                     "" index 6 using 1:2 with lines title "K_6(x; p, N)", \
                     "" index 7 using 1:2 with lines title "K_7(x; p, N)", \
                     "" index 8 using 1:2 with lines title "K_8(x; p, N)", \
                     "" index 9 using 1:2 with lines title "K_9(x; p, N)", \
                     "" index 10 using 1:2 with lines title "K_{10}(x; p, N)"