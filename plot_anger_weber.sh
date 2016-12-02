
gnuplot
set xzeroaxis
set yzeroaxis
set grid

plot [-8.0:8.0][-1.0:1.0] "test_anger_weber.txt" index 0 using 1:2 with lines title "J_0(x)"
plot [-8.0:8.0][-1.0:1.0] "test_anger_weber.txt" index 0 using 1:3 with lines title "E_0(x)"


plot [-8.0:8.0][-1.0:1.0] \
"test_anger_weber.txt" index 0 using 1:2 with lines title "J_0(x)", \
                    "" index 0 using 1:3 with lines title "H_0(x)", \
                    "" index 1 using 1:2 with lines title "J_0.5(x)", \
                    "" index 1 using 1:3 with lines title "H_0.5(x)", \
                    "" index 2 using 1:2 with lines title "J_1(x)", \
                    "" index 2 using 1:3 with lines title "H_1(x)", \
                    "" index 3 using 1:2 with lines title "J_1.5(x)", \
                    "" index 3 using 1:3 with lines title "H_1.5(x)"

plot [-8.0:8.0][-1.0:1.0] \
"test_anger_weber.txt" index 0 using 1:2 with lines title "J_0(x)", \
                    "" index 1 using 1:2 with lines title "J_0.5(x)", \
                    "" index 2 using 1:2 with lines title "J_1(x)", \
                    "" index 3 using 1:2 with lines title "J_1.5(x)"

plot [-8.0:8.0][-1.0:1.0] \
"test_anger_weber.txt" index 0 using 1:3 with lines title "H_0(x)", \
                    "" index 1 using 1:3 with lines title "H_0.5(x)", \
                    "" index 2 using 1:3 with lines title "H_1(x)", \
                    "" index 3 using 1:3 with lines title "H_1.5(x)"
