
gnuplot

load 'setings.gp'

set title "Associated laguerre polynomial L_1^{({/Symbol a})}(x)"
set xlabel "x"
plot [-10.0:20.0][-60.0:60.0] \
    "test_assoc_laguerre.txt" index 5 using 1:2 with lines title "L_1^{(0)}(x)", \
                           "" index 6 using 1:2 with lines title "L_1^{(1/3)}(x)", \
                           "" index 7 using 1:2 with lines title "L_1^{(1/2)}(x)", \
                           "" index 8 using 1:2 with lines title "L_1^{(2/3)}(x)", \
                           "" index 9 using 1:2 with lines title "L_1^{(1)}(x)"

set title "Associated laguerre polynomial L_2^{({/Symbol a})}(x)"
set xlabel "x"
plot [-10.0:20.0][-60.0:60.0] \
    "test_assoc_laguerre.txt" index 10 using 1:2 with lines title "L_2^{(0)}(x)", \
                           "" index 11 using 1:2 with lines title "L_2^{(1/3)}(x)", \
                           "" index 12 using 1:2 with lines title "L_2^{(1/2)}(x)", \
                           "" index 13 using 1:2 with lines title "L_2^{(2/3)}(x)", \
                           "" index 14 using 1:2 with lines title "L_2^{(1)}(x)"

set title "Associated laguerre polynomial L_5^{({/Symbol a})}(x)"
set xlabel "x"
plot [-10.0:20.0][-60.0:60.0] \
    "test_assoc_laguerre.txt" index 15 using 1:2 with lines title "L_5^{(0)}(x)", \
                           "" index 16 using 1:2 with lines title "L_5^{(1/3)}(x)", \
                           "" index 17 using 1:2 with lines title "L_5^{(1/2)}(x)", \
                           "" index 18 using 1:2 with lines title "L_5^{(2/3)}(x)", \
                           "" index 19 using 1:2 with lines title "L_5^{(1)}(x)"

set title "Associated laguerre polynomial L_1^{({/Symbol a})}(x)"
set xlabel "x"
plot [-10.0:20.0][-60.0:60.0] \
    "test_assoc_laguerre.txt" index 23 using 1:2 with lines title "L_1^{(-1)}(x)", \
                           "" index 24 using 1:2 with lines title "L_1^{(-2)}(x)", \
                           "" index 25 using 1:2 with lines title "L_1^{(-3)}(x)"

set title "Associated laguerre polynomial L_2^{({/Symbol a})}(x)"
set xlabel "x"
plot [-10.0:20.0][-60.0:60.0] \
    "test_assoc_laguerre.txt" index 26 using 1:2 with lines title "L_2^{(-1)}(x)", \
                           "" index 27 using 1:2 with lines title "L_2^{(-2)}(x)", \
                           "" index 28 using 1:2 with lines title "L_2^{(-3)}(x)"

set title "Associated laguerre polynomial L_5^{({/Symbol a})}(x)"
set xlabel "x"
plot [-10.0:20.0][-60.0:60.0] \
    "test_assoc_laguerre.txt" index 29 using 1:2 with lines title "L_5^{(-1)}(x)", \
                           "" index 30 using 1:2 with lines title "L_5^{(-2)}(x)", \
                           "" index 31 using 1:2 with lines title "L_5^{(-3)}(x)"
