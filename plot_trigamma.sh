
/usr/local/bin/gnuplot

set xzeroaxis
set yzeroaxis
set grid

plot [0.0:10.0][0.0:1.0] "test_trigamma.txt" using 1:2 with lines title "phi_1(x)"

plot [0.0:10.0][0.0:1.0] "test_trigamma.txt" using 1:2 with lines title "phi_2(x)"
