
gnuplot

set termoption enhanced

set xzeroaxis
set yzeroaxis
set grid

# We're right on with Boost.
#                             "" index 0 using 2:4 with lines title "Z(k, {/Symbol f}) Boost"

# k = 0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 0.99, 1.0

set title "Jazobi zeta function Z(k, {/Symbol f})"
set xlabel "{/Symbol f}"
plot [-1.57079:+1.57079][-.5:.5] \
  "output/test_jacobi_zeta.txt" index  0 using 2:3 with lines title "Z(k=0, {/Symbol f})", \
                             "" index  1 using 2:3 with lines title "Z(k=0.1, {/Symbol f})", \
                             "" index  2 using 2:3 with lines title "Z(k=0.2, {/Symbol f})", \
                             "" index  3 using 2:3 with lines title "Z(k=0.3, {/Symbol f})", \
                             "" index  4 using 2:3 with lines title "Z(k=0.4, {/Symbol f})", \
                             "" index  5 using 2:3 with lines title "Z(k=0.5, {/Symbol f})", \
                             "" index  6 using 2:3 with lines title "Z(k=0.6, {/Symbol f})", \
                             "" index  7 using 2:3 with lines title "Z(k=0.7, {/Symbol f})", \
                             "" index  8 using 2:3 with lines title "Z(k=0.8, {/Symbol f})", \
                             "" index  9 using 2:3 with lines title "Z(k=0.9, {/Symbol f})", \
                             "" index 10 using 2:3 with lines title "Z(k=0.99, {/Symbol f})", \
                             "" index 11 using 2:3 with lines title "Z(k=1, {/Symbol f})"

set term push
set term png
set output "jacobi_zeta.png"
replot
set term pop
