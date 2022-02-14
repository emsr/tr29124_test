
gnuplot

set termoption enhanced

set xzeroaxis
set yzeroaxis
set grid

# lambda = 0

#filename = "test_coulomb.txt"
filename = "test_coulomb_new.txt"
#filename = "test_coulomb_gsl.txt"

set title "Coulomb functions and derivatives, {/Symbol l} = 0, {/Symbol h} = -2"
set xlabel "{/Symbol r}"
plot [0.0:20.0][-2.0:3.0] \
    filename index 0 using 1:2 with lines title "F_{L}({/Symbol h}, {/Symbol r})", \
                    "" index 0 using 1:3 with lines title "G_{L}({/Symbol h}, {/Symbol r})", \
                    "" index 0 using 1:4 with lines title "F'_{L}({/Symbol h}, {/Symbol r})", \
                    "" index 0 using 1:5 with lines title "G'_{L}({/Symbol h}, {/Symbol r})"

set title "Coulomb functions and derivatives, {/Symbol l} = 0, {/Symbol h} = 0"
set xlabel "{/Symbol r}"
plot [0.0:20.0][-1.0:1.5] \
    filename index 1 using 1:2 with lines title "F_{L}({/Symbol h}, {/Symbol r})", \
                    "" index 1 using 1:3 with lines title "G_{L}({/Symbol h}, {/Symbol r})", \
                    "" index 1 using 1:4 with lines title "F'_{L}({/Symbol h}, {/Symbol r})", \
                    "" index 1 using 1:5 with lines title "G'_{L}({/Symbol h}, {/Symbol r})"

set title "Coulomb functions and derivatives, {/Symbol l} = 0, {/Symbol h} = 2"
set xlabel "{/Symbol r}"
plot [0.0:20.0][-1.0:3.0] \
    filename index 2 using 1:2 with lines title "F_{L}({/Symbol h}, {/Symbol r})", \
                    "" index 2 using 1:3 with lines title "G_{L}({/Symbol h}, {/Symbol r})", \
                    "" index 2 using 1:4 with lines title "F'_{L}({/Symbol h}, {/Symbol r})", \
                    "" index 2 using 1:5 with lines title "G'_{L}({/Symbol h}, {/Symbol r})"

set title "Coulomb functions and derivatives, {/Symbol l} = 0, {/Symbol h} = 10"
set xlabel "{/Symbol r}"
plot [0.0:20.0][-1.0:3.0] \
    filename index 3 using 1:2 with lines title "F_{L}({/Symbol h}, {/Symbol r})", \
                    "" index 3 using 1:3 with lines title "G_{L}({/Symbol h}, {/Symbol r})", \
                    "" index 3 using 1:4 with lines title "F'_{L}({/Symbol h}, {/Symbol r})", \
                    "" index 3 using 1:5 with lines title "G'_{L}({/Symbol h}, {/Symbol r})"

# lambda = 1/2

set title "Coulomb functions and derivatives, {/Symbol l} = 1/2, {/Symbol h} = -2"
set xlabel "{/Symbol r}"
plot [0.0:20.0][-2.0:3.0] \
    filename index 4 using 1:2 with lines title "F_{L}({/Symbol h}, {/Symbol r})", \
                    "" index 4 using 1:3 with lines title "G_{L}({/Symbol h}, {/Symbol r})", \
                    "" index 4 using 1:4 with lines title "F'_{L}({/Symbol h}, {/Symbol r})", \
                    "" index 4 using 1:5 with lines title "G'_{L}({/Symbol h}, {/Symbol r})"

set title "Coulomb functions and derivatives, {/Symbol l} = 1/2, {/Symbol h} = 0"
set xlabel "{/Symbol r}"
plot [0.0:20.0][-1.0:1.5] \
    filename index 5 using 1:2 with lines title "F_{L}({/Symbol h}, {/Symbol r})", \
                    "" index 5 using 1:3 with lines title "G_{L}({/Symbol h}, {/Symbol r})", \
                    "" index 5 using 1:4 with lines title "F'_{L}({/Symbol h}, {/Symbol r})", \
                    "" index 5 using 1:5 with lines title "G'_{L}({/Symbol h}, {/Symbol r})"

set title "Coulomb functions and derivatives, {/Symbol l} = 1/2, {/Symbol h} = 2"
set xlabel "{/Symbol r}"
plot [0.0:20.0][-1.0:3.0] \
    filename index 6 using 1:2 with lines title "F_{L}({/Symbol h}, {/Symbol r})", \
                    "" index 6 using 1:3 with lines title "G_{L}({/Symbol h}, {/Symbol r})", \
                    "" index 6 using 1:4 with lines title "F'_{L}({/Symbol h}, {/Symbol r})", \
                    "" index 6 using 1:5 with lines title "G'_{L}({/Symbol h}, {/Symbol r})"

set title "Coulomb functions and derivatives, {/Symbol l} = 1/2, {/Symbol h} = 10"
set xlabel "{/Symbol r}"
plot [0.0:20.0][-1.0:3.0] \
    filename index 7 using 1:2 with lines title "F_{L}({/Symbol h}, {/Symbol r})", \
                    "" index 7 using 1:3 with lines title "G_{L}({/Symbol h}, {/Symbol r})", \
                    "" index 7 using 1:4 with lines title "F'_{L}({/Symbol h}, {/Symbol r})", \
                    "" index 7 using 1:5 with lines title "G'_{L}({/Symbol h}, {/Symbol r})"

# lambda = 1

set title "Coulomb functions and derivatives, {/Symbol l} = 1, {/Symbol h} = -2"
set xlabel "{/Symbol r}"
plot [0.0:20.0][-2.0:3.0] \
    filename index 8 using 1:2 with lines title "F_{L}({/Symbol h}, {/Symbol r})", \
                    "" index 8 using 1:3 with lines title "G_{L}({/Symbol h}, {/Symbol r})", \
                    "" index 8 using 1:4 with lines title "F'_{L}({/Symbol h}, {/Symbol r})", \
                    "" index 8 using 1:5 with lines title "G'_{L}({/Symbol h}, {/Symbol r})"

set title "Coulomb functions and derivatives, {/Symbol l} = 1, {/Symbol h} = 0"
set xlabel "{/Symbol r}"
plot [0.0:20.0][-1.0:1.5] \
    filename index 9 using 1:2 with lines title "F_{L}({/Symbol h}, {/Symbol r})", \
                    "" index 9 using 1:3 with lines title "G_{L}({/Symbol h}, {/Symbol r})", \
                    "" index 9 using 1:4 with lines title "F'_{L}({/Symbol h}, {/Symbol r})", \
                    "" index 9 using 1:5 with lines title "G'_{L}({/Symbol h}, {/Symbol r})"

set title "Coulomb functions and derivatives, {/Symbol l} = 1, {/Symbol h} = 2"
set xlabel "{/Symbol r}"
plot [0.0:20.0][-1.0:3.0] \
    filename index 10 using 1:2 with lines title "F_{L}({/Symbol h}, {/Symbol r})", \
                    "" index 10 using 1:3 with lines title "G_{L}({/Symbol h}, {/Symbol r})", \
                    "" index 10 using 1:4 with lines title "F'_{L}({/Symbol h}, {/Symbol r})", \
                    "" index 10 using 1:5 with lines title "G'_{L}({/Symbol h}, {/Symbol r})"

set title "Coulomb functions and derivatives, {/Symbol l} = 1, {/Symbol h} = 10"
set xlabel "{/Symbol r}"
plot [0.0:20.0][-1.0:3.0] \
    filename index 11 using 1:2 with lines title "F_{L}({/Symbol h}, {/Symbol r})", \
                    "" index 11 using 1:3 with lines title "G_{L}({/Symbol h}, {/Symbol r})", \
                    "" index 11 using 1:4 with lines title "F'_{L}({/Symbol h}, {/Symbol r})", \
                    "" index 11 using 1:5 with lines title "G'_{L}({/Symbol h}, {/Symbol r})"
