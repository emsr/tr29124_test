
gnuplot

set termoption enhanced

set xzeroaxis
set yzeroaxis
set grid

plot [-5.0:20.0][0.0:125.0] "test_debye.txt" index 0 using 1:2 with lines title "Debye_1(x)", \
                                         "" index 0 using 1:3 with lines title "Debye_2(x)", \
                                         "" index 0 using 1:4 with lines title "Debye_3(x)", \
                                         "" index 0 using 1:5 with lines title "Debye_4(x)", \
                                         "" index 0 using 1:6 with lines title "Debye_5(x)"

