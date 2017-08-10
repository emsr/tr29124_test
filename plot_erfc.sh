
gnuplot
set termoption enhanced
set xzeroaxis
set yzeroaxis
set grid

set xlabel "x"
plot [-2.0:2.0][-1.0:2.0] \
    "test_erfc.txt" index 0 using 1:2 with lines title "erfc(x)", \
                 "" index 1 using 1:2 with lines title "erf(x)"

set xlabel "x"
plot [-2.0:5.0][0.0:2.0] \
  "test_erfc.txt" index 2 using 1:2 with lines title "i^{-1}erfc(x)", \
               "" index 2 using 1:3 with lines title "i^0erfc(x)", \
               "" index 2 using 1:4 with lines title "i^1erfc(x)", \
               "" index 2 using 1:5 with lines title "i^2erfc(x)", \
               "" index 2 using 1:6 with lines title "i^3erfc(x)", \
               "" index 2 using 1:7 with lines title "i^4erfc(x)", \
               "" index 2 using 1:8 with lines title "i^5erfc(x)"
