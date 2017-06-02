
gnuplot

set termoption enhanced

set xzeroaxis
set yzeroaxis
set grid

set title "Coulomb function"
set xlabel "x"
plot [0.0:20.0][-5.0:10.0] \
    "test_coulomb.txt" index 0 using 1:2 with lines title "F_{L}({/Symbol h}, {/Symbol r})", \
                    "" index 0 using 1:3 with lines title "G_{L}({/Symbol h}, {/Symbol r})", \
                    "" index 0 using 1:4 with lines title "F'_{L}({/Symbol h}, {/Symbol r})", \
                    "" index 0 using 1:5 with lines title "G'_{L}({/Symbol h}, {/Symbol r})"

set title "Coulomb function"
set xlabel "x"
plot [0.0:20.0][-5.0:10.0] \
    "test_coulomb.txt" index 1 using 1:2 with lines title "F_{L}({/Symbol h}, {/Symbol r})", \
                    "" index 1 using 1:3 with lines title "G_{L}({/Symbol h}, {/Symbol r})", \
                    "" index 1 using 1:4 with lines title "F'_{L}({/Symbol h}, {/Symbol r})", \
                    "" index 1 using 1:5 with lines title "G'_{L}({/Symbol h}, {/Symbol r})"

set title "Coulomb function"
set xlabel "x"
plot [0.0:20.0][-5.0:10.0] \
    "test_coulomb.txt" index 2 using 1:2 with lines title "F_{L}({/Symbol h}, {/Symbol r})", \
                    "" index 2 using 1:3 with lines title "G_{L}({/Symbol h}, {/Symbol r})", \
                    "" index 2 using 1:4 with lines title "F'_{L}({/Symbol h}, {/Symbol r})", \
                    "" index 2 using 1:5 with lines title "G'_{L}({/Symbol h}, {/Symbol r})"

set title "Coulomb function"
set xlabel "x"
plot [0.0:20.0][-5.0:10.0] \
    "test_coulomb.txt" index 3 using 1:2 with lines title "F_{L}({/Symbol h}, {/Symbol r})", \
                    "" index 3 using 1:3 with lines title "G_{L}({/Symbol h}, {/Symbol r})", \
                    "" index 3 using 1:4 with lines title "F'_{L}({/Symbol h}, {/Symbol r})", \
                    "" index 3 using 1:5 with lines title "G'_{L}({/Symbol h}, {/Symbol r})"
