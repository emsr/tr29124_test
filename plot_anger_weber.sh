
gnuplot
set termoption enhanced
set xzeroaxis
set yzeroaxis
set grid

plot [-8.0:8.0][-1.0:1.0] "test_anger_weber.txt" index 0 using 1:2 with lines title "J_0(x)"
plot [-8.0:8.0][-1.0:1.0] "test_anger_weber.txt" index 0 using 1:3 with lines title "E_0(x)"


plot [-8.0:8.0][-1.0:1.0] \
"test_anger_weber.txt" index 0 using 1:2 with lines title "{/:Bold J}_0(x)", \
                    "" index 0 using 1:3 with lines title "{/:Bold H}_0(x)", \
                    "" index 1 using 1:2 with lines title "{/:Bold J}_{1/2}(x)", \
                    "" index 1 using 1:3 with lines title "{/:Bold H}_{1/2}(x)", \
                    "" index 2 using 1:2 with lines title "{/:Bold J}_1(x)", \
                    "" index 2 using 1:3 with lines title "{/:Bold H}_1(x)", \
                    "" index 3 using 1:2 with lines title "{/:Bold J}_{3/2}(x)", \
                    "" index 3 using 1:3 with lines title "{/:Bold H}_{3/2}(x)"

plot [-8.0:8.0][-1.0:1.0] \
"test_anger_weber.txt" index 0 using 1:2 with lines title "{/:Bold J}_0(x)", \
                    "" index 1 using 1:2 with lines title "{/:Bold J}_{1/2}(x)", \
                    "" index 2 using 1:2 with lines title "{/:Bold J}_1(x)", \
                    "" index 3 using 1:2 with lines title "{/:Bold J}_{3/2}(x)"

plot [-8.0:8.0][-1.0:1.0] \
"test_anger_weber.txt" index 0 using 1:2 with lines title "{/:Bold J}_0(x)", \
                    "" index 1 using 1:2 with lines title "{/:Bold J}_{1/2}(x)", \
                    "" index 2 using 1:2 with lines title "{/:Bold J}_1(x)", \
                    "" index 3 using 1:2 with lines title "{/:Bold J}_{3/2}(x)", \
                    "" index 4 using 1:2 with lines title "{/:Bold J}_{1.999}(x)", \
                    "" index 5 using 1:2 with lines title "{/:Bold J}_{2}(x)", \
                    "" index 6 using 1:2 with lines title "{/:Bold J}_{2.999}(x)", \
                    "" index 7 using 1:2 with lines title "{/:Bold J}_{3}(x)"

plot [-8.0:8.0][-1.0:1.0] \
"test_anger_weber.txt" index 0 using 1:3 with lines title "{/:Bold H}_0(x)", \
                    "" index 1 using 1:3 with lines title "{/:Bold H}_{1/2}(x)", \
                    "" index 2 using 1:3 with lines title "{/:Bold H}_1(x)", \
                    "" index 3 using 1:3 with lines title "{/:Bold H}_{3/2}(x)", \
                    "" index 4 using 1:3 with lines title "{/:Bold H}_{1.999}(x)", \
                    "" index 5 using 1:3 with lines title "{/:Bold H}_{2}(x)", \
                    "" index 6 using 1:3 with lines title "{/:Bold H}_{2.999}(x)", \
                    "" index 7 using 1:3 with lines title "{/:Bold H}_{3}(x)"
