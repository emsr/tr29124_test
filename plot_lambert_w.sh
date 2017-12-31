
gnuplot
set termoption enhanced
set xzeroaxis
set yzeroaxis
set grid

set title "Lambert W function W(z)"
set xlabel "p"
plot [-0.34:2.8][-1.1:1.1] "test_lambert_w.txt" index 0 using 1:2 with lines title "W(z)"

set title "Lambert W function W(z)"
set xlabel "p"
plot [2.0:0.0][-2.0:2.0] "test_lambert_w.txt" index 1 using 1:2 with lines title "W(z)"
