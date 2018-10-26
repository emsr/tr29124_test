
gnuplot

set termoption enhanced

set xzeroaxis
set yzeroaxis
set grid

set title "Digamma  {/Symbol y}(x)"
set xlabel "x"
plot [-10.0:10.0][-10.0:10.0] \
                    "test_digamma.txt" index 0 using 1:2 with lines title "{/Symbol y}(x)", \
                    "" index 0 using 1:3 with lines title "{/Symbol y}(x) GSL"

# The harmonic numbers don't plot :-(
set title "Digamma  {/Symbol y}(x)"
set xlabel "x"
plot [0.0:50.0][0.0:4.0] \
                    "test_digamma.txt" index 1 using 1:2 with lines title "{/Symbol y}(x)", \
                    "" index 1 using 1:3 with lines title "{/Symbol y}(x) GSL", \
                    "" index 1 using 1:4 with lines title "H_n"

# 
set title "Polygamma  {/Symbol y}^{(m)}(x)"
set xlabel "x"
plot [-4.0:4.0][-40.0:40.0] \
    "test_polygamma.txt" index 3 using 1:2 with lines title "{/Symbol y}^{(0)}(x)", \
                "" index 3 using 1:3 with lines title "{/Symbol y}^{(1)}(x)", \
                "" index 3 using 1:4 with lines title "{/Symbol y}^{(2)}(x)", \
                "" index 3 using 1:5 with lines title "{/Symbol y}^{(3)}(x)", \
                "" index 3 using 1:6 with lines title "{/Symbol y}^{(4)}(x)", \
                "" index 3 using 1:7 with lines title "{/Symbol y}^{(5)}(x)"
