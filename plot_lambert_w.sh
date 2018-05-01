
gnuplot

set termoption enhanced
set xzeroaxis
set yzeroaxis
set grid

set title "Lambert W function W_0(z)"
set xlabel "p"
plot [-0.5:2.8][-1.1:1.1] "test_lambert_w.txt" index 0 using 1:2 with lines title "W_0(z)"

set title "Lambert W function W_{-1}(z)"
set xlabel "p"
plot [-0.5:2.8][-5.0:-0.9] "test_lambert_w.txt" index 1 using 1:2 with lines title "W_{-1}(z)"

set title "Lambert W function W(z)"
set xlabel "p"
plot [-0.5:2.8][-5.0:1.1] \
    "test_lambert_w.txt" index 0 using 1:2 with lines title "W_0(z)", \
                      "" index 1 using 1:2 with lines title "W_{-1}(z)"
