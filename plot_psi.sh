
gnuplot

set termoption enhanced

set xzeroaxis
set yzeroaxis
set grid

set title "Digamma {/Symbol y}(x)"
set xlabel "x"
plot [-10.0:10.0][-10.0:10.0] \
                    "test_psi.txt" index 0 using 1:2 with lines title "{/Symbol y}(x)", \
                    "" index 0 using 1:3 with lines title "{/Symbol y}(x) GSL"

# The harmonic numbers don't plot :-(
set title "Digamma {/Symbol y}(x)"
set xlabel "x"
plot [0.0:50.0][0.0:4.0] \
                    "test_psi.txt" index 1 using 1:2 with lines title "{/Symbol y}(x)", \
                    "" index 1 using 1:3 with lines title "{/Symbol y}(x) GSL", \
                    "" index 1 using 1:4 with lines title "H_n"
