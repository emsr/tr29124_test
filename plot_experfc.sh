
gnuplot
set termoption enhanced
set xzeroaxis
set yzeroaxis
set grid

set xlabel "x"
plot [-2.0:15.0][0.0:2.0] "test_experfc.txt" index 0 using 1:2 with lines title "experfc(x)"
