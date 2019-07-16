
gnuplot
set termoption enhanced
set xzeroaxis
set yzeroaxis
set grid

set title "Chebyshev polynomials and their derivatives, n = 1"
set xlabel "x"
plot [-1.0:1.0][-3.0:3.0] \
    "test_chebyshev.txt" index 1 using 1:2 with lines title "T(x)", \
                      "" index 1 using 1:3 with lines title "T'(x)", \
                      "" index 1 using 1:4 with lines title "U(x)", \
                      "" index 1 using 1:5 with lines title "U'(x)", \
                      "" index 1 using 1:6 with lines title "V(x)", \
                      "" index 1 using 1:7 with lines title "V'(x)", \
                      "" index 1 using 1:8 with lines title "W(x)", \
                      "" index 1 using 1:9 with lines title "W'(x)"

set term push
set term png
set output "chebyshev.png"
replot
set term pop
