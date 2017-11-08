# http://pool-serv1.mathematik.uni-kassel.de/caop/CAOP/plot_charlier.html

gnuplot

set termoption enhanced
set xzeroaxis
set yzeroaxis
set grid

set title "tgamma(x) - {/Symbol G}(x)"
set xlabel "x"
plot [-5.0:+5.0][-100.0:+100.0] \
    "tgammaq-patched.dat" index 0 using 1:2 with lines title "{/Symbol G}(x) - Patched", \
    "tgammaq-buggy.dat"   index 0 using 1:2 with lines title "{/Symbol G}(x) - Buggy"
