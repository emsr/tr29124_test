
gnuplot

set termoption enhanced

set xzeroaxis
set yzeroaxis
set grid

set title "Bose-Einstein integral G_s(x)"
set xlabel "-x"
plot [-10.0:10.0][-5.0:10.0] \
                    "test_bose_einstein.txt" index 0 using 1:2 with lines title "G_{0}(x)", \
                    "" index 1 using 1:2 with lines title "G_{1/2}(x)", \
                    "" index 2 using 1:2 with lines title "G_{1}(x)", \
                    "" index 3 using 1:2 with lines title "G_{3/2}(x)", \
                    "" index 4 using 1:2 with lines title "G_{2}(x)", \
                    "" index 5 using 1:2 with lines title "G_{3}(x)", \
                    "" index 6 using 1:2 with lines title "G_{4}(x)", \
                    "" index 7 using 1:2 with lines title "G_{5}(x)"

set title "Fermi-Dirac integral F_s(x)"
set xlabel "-x"
plot [-10.0:10.0][0.0:10.0] \
                    "test_fermi_dirac.txt" index 0 using 1:2 with lines title "F_{0}(x)", \
                    "" index 1 using 1:2 with lines title "F_{1/2}(x)", \
                    "" index 2 using 1:2 with lines title "F_{1}(x)", \
                    "" index 3 using 1:2 with lines title "F_{3/2}(x)", \
                    "" index 4 using 1:2 with lines title "F_{2}(x)", \
                    "" index 5 using 1:2 with lines title "F_{3}(x)", \
                    "" index 6 using 1:2 with lines title "F_{4}(x)", \
                    "" index 7 using 1:2 with lines title "F_{5}(x)"
