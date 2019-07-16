# http://pool-serv1.mathematik.uni-kassel.de/caop/CAOP/plot_charlier.html

gnuplot

load 'setings.gp'

set title "Bernoulli Polynomial B_n(x)"
set xlabel "x"
plot [0.0:9.5][-25.0:16.0] \
    "test_bernoulli.txt" index 6 using 1:2 with lines title "B_0(x)"

set term push
set term png
set output "bernoulli.png"
replot
set term pop
