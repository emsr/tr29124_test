
gnuplot

set termoption enhanced

set xzeroaxis
set yzeroaxis
set grid

asymp(x, c) = c

# Figure 1
set title "Clausen function Cl_2(x)"
set xlabel "x"
plot [-10.0:10.0][-1.2:1.2] \
                    "test_clausen.txt" index 0 using 1:2 with lines title "Re[Cl_1(x)]", \
                    "" index 0 using 1:3 with lines title "Im[Cl_1(x)]", \
                    "" index 0 using 1:4 with lines title "Re[Cl_2(x)]", \
                    "" index 0 using 1:5 with lines title "Im[Cl_2(x)]", \
                    "" index 0 using 1:6 with lines title "Im[C_2(x) GSL]"
