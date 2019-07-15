
gnuplot

set termoption enhanced

set xzeroaxis
set yzeroaxis
set grid

set title "Complete elliptic function of the first kind K(k)"
set xlabel "x"
plot [0.0:1.0][1.5:3.5] \
                    "output/test_legendre_ellint.txt" index 0 using 1:2 with lines title "K(k)", \
                    "" index 0 using 1:3 with lines title "K(k) GSL"

set title "Incomplete elliptic function of the first kind F(k,{/Symbol f})"
set xlabel "{/Symbol f}"
plot [0.0:1.6][0.0:3.5] \
                    "output/test_legendre_ellint.txt" index 1 using 1:2 with lines title "F(0.0,{/Symbol f})", \
                    "" index 1 using 1:5 with lines title "F(0.5,{/Symbol f})", \
                    "" index 1 using 1:8 with lines title "F(0.75,{/Symbol f})", \
                    "" index 1 using 1:11 with lines title "F(0.9,{/Symbol f})", \
                    "" index 1 using 1:14 with lines title "F(0.95,{/Symbol f})", \
                    "" index 1 using 1:17 with lines title "F(0.99,{/Symbol f})"

set title "Complete elliptic function of the second kind E(k)"
set xlabel "x"
plot [0.0:1.0][1.0:1.6] \
                    "output/test_legendre_ellint.txt" index 2 using 1:2 with lines title "E(k)", \
                    "" index 2 using 2:3 with lines title "E(k) GSL"

set title "Incomplete elliptic function of the second kind E(k,{/Symbol f})"
set xlabel "{/Symbol f}"
plot [0.0:1.6][0.0:1.6] \
                    "output/test_legendre_ellint.txt" index 3 using 1:2 with lines title "E(0.0,{/Symbol f})", \
                    "" index 3 using 1:5 with lines title "E(0.5,{/Symbol f})", \
                    "" index 3 using 1:8 with lines title "E(0.75,{/Symbol f})", \
                    "" index 3 using 1:11 with lines title "E(0.9,{/Symbol f})", \
                    "" index 3 using 1:14 with lines title "E(0.95,{/Symbol f})", \
                    "" index 3 using 1:17 with lines title "E(0.99,{/Symbol f})"
