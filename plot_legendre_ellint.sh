
gnuplot

set termoption enhanced

set xzeroaxis
set yzeroaxis
set grid

set title "Complete elliptic function of the first kind K(k)"
set xlabel "x"
plot [0.0:1.0][0.0:4.0] \
                    "test_legendre_ellint.txt" index 0 using 1:2 with lines title "K(k)", \
                    "" index 0 using 1:3 with lines title "K(k) GSL"

set title "Complete elliptic function of the first kind K(k)"
set xlabel "x"
plot [0.0:1.0][0.0:4.0] \
                    "test_legendre_ellint.txt" index 1 using 1:2 with lines title "K(k)", \
                    "" index 1 using 1:3 with lines title "K(k) GSL"
