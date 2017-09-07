
gnuplot
set termoption enhanced
set xzeroaxis
set yzeroaxis
set grid

set xlabel "p"
plot [-1.0:1.0][-2.0:2.0] "test_inv_erf.txt" index 0 using 1:2 with lines title "erf^{-1}(p)"

set xlabel "p"
plot [2.0:0.0][-2.0:2.0] "test_inv_erf.txt" index 1 using 1:2 with lines title "erfc^{-1}(p)"
