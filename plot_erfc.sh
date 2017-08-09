
gnuplot
set termoption enhanced
set xzeroaxis
set yzeroaxis
set grid

set xlabel "x"
plot [-2.0:2.0][-1.0:2.0] \
    "test_erfc.txt" index 0 using 1:2 with lines title "erfc(x)", \
                 "" index 1 using 1:2 with lines title "erf(x)"
