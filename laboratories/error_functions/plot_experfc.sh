
gnuplot
set termoption enhanced
set xzeroaxis
set yzeroaxis
set grid

set title "Scaled complementary error function experfc(x) = e^{x^2}erfc(x)"
set xlabel "x"
plot [-2.0:15.0][0.0:2.0] "test_experfc.txt" index 0 using 1:2 with lines title "experfc(x)"
