
gnuplot

set xzeroaxis
set yzeroaxis
set grid

plot [0.0:20.0][0.0:1.0] "test_trigamma.txt" index 0 using 1:2 with lines title "phi_1(x)"

plot [0.0:20.0][-0.5:0.0] "test_trigamma.txt" index 1 using 1:2 with lines title "phi_2(x)"
