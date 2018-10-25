gnuplot

load 'setings.gp'

set title "Coulomb Functions and Derivatives {/Symbol l} = 1, {/Symbol h} = -2"
set xlabel "{/Symbol r}"
plot [0.0:20.0][-1.6:1.6] \
    "run_coulfg.txt" index 0 using 1:2 with lines title "F_{1}(-2, {/Symbol r})", \
                  "" index 0 using 1:3 with lines title "G_{1}(-2, {/Symbol r})", \
                  "" index 0 using 1:4 with lines title "F'_{1}(-2, {/Symbol r})", \
                  "" index 0 using 1:5 with lines title "G'_{1}(-2, {/Symbol r})"

set title "Coulomb Functions and Derivatives {/Symbol l} = 1, {/Symbol h} = 0"
set xlabel "{/Symbol r}"
plot [0.0:20.0][-1.6:1.6] \
    "run_coulfg.txt" index 1 using 1:2 with lines title "F_{1}(0, {/Symbol r})", \
                  "" index 1 using 1:3 with lines title "G_{1}(0, {/Symbol r})", \
                  "" index 1 using 1:4 with lines title "F'_{1}(0, {/Symbol r})", \
                  "" index 1 using 1:5 with lines title "G'_{1}(0, {/Symbol r})"

set title "Coulomb Functions and Derivatives {/Symbol l} = 1, {/Symbol h} = 2"
set xlabel "{/Symbol r}"
plot [0.0:20.0][-1.6:1.6] \
    "run_coulfg.txt" index 2 using 1:2 with lines title "F_{1}(2, {/Symbol r})", \
                  "" index 2 using 1:3 with lines title "G_{1}(2, {/Symbol r})", \
                  "" index 2 using 1:4 with lines title "F'_{1}(2, {/Symbol r})", \
                  "" index 2 using 1:5 with lines title "G'_{1}(2, {/Symbol r})"

set title "Coulomb Functions and Derivatives {/Symbol l} = 1, {/Symbol h} = 4"
set xlabel "{/Symbol r}"
plot [0.0:20.0][-1.6:1.6] \
    "run_coulfg.txt" index 3 using 1:2 with lines title "F_{1}(4, {/Symbol r})", \
                  "" index 3 using 1:3 with lines title "G_{1}(4, {/Symbol r})", \
                  "" index 3 using 1:4 with lines title "F'_{1}(4, {/Symbol r})", \
                  "" index 3 using 1:5 with lines title "G'_{1}(4, {/Symbol r})"

set title "Coulomb Functions and Derivatives {/Symbol l} = 1, {/Symbol h} = 6"
set xlabel "{/Symbol r}"
plot [0.0:20.0][-1.6:1.6] \
    "run_coulfg.txt" index 4 using 1:2 with lines title "F_{1}(6, {/Symbol r})", \
                  "" index 4 using 1:3 with lines title "G_{1}(6, {/Symbol r})", \
                  "" index 4 using 1:4 with lines title "F'_{1}(6, {/Symbol r})", \
                  "" index 4 using 1:5 with lines title "G'_{1}(6, {/Symbol r})"

set title "Coulomb Functions and Derivatives {/Symbol l} = 1, {/Symbol h} = 8"
set xlabel "x"
plot [0.0:20.0][-1.6:1.6] \
    "run_coulfg.txt" index 5 using 1:2 with lines title "F_{1}({/Symbol h}, {/Symbol r})", \
                  "" index 5 using 1:3 with lines title "G_{1}({/Symbol h}, {/Symbol r})", \
                  "" index 5 using 1:4 with lines title "F'_{1}({/Symbol h}, {/Symbol r})", \
                  "" index 5 using 1:5 with lines title "G'_{1}({/Symbol h}, {/Symbol r})"

set title "Coulomb Functions and Derivatives {/Symbol l} = 1, {/Symbol h} = 10"
set xlabel "{/Symbol r}"
plot [0.0:20.0][-1.6:1.6] \
    "run_coulfg.txt" index 6 using 1:2 with lines title "F_{1}(10, {/Symbol r})", \
                  "" index 6 using 1:3 with lines title "G_{1}(10, {/Symbol r})", \
                  "" index 6 using 1:4 with lines title "F'_{1}(10, {/Symbol r})", \
                  "" index 6 using 1:5 with lines title "G'_{1}(10, {/Symbol r})"



set title "Coulomb Function F_{{/Symbol l}=1}({/Symbol h}, {/Symbol r})"
set xlabel "{/Symbol r}"
plot [0.0:20.0][-1.6:1.6] \
    "run_coulfg.txt" index 0 using 1:2 with lines title "F_{1}(-2, {/Symbol r})", \
                  "" index 1 using 1:2 with lines title "F_{1}(0, {/Symbol r})", \
                  "" index 2 using 1:2 with lines title "F_{1}(2, {/Symbol r})", \
                  "" index 3 using 1:2 with lines title "F_{1}(4, {/Symbol r})", \
                  "" index 4 using 1:2 with lines title "F_{1}(6, {/Symbol r})", \
                  "" index 5 using 1:2 with lines title "F_{1}(8, {/Symbol r})", \
                  "" index 6 using 1:2 with lines title "F_{1}(10, {/Symbol r})"

set title "Coulomb Function G_{{/Symbol l}=1}({/Symbol h}, {/Symbol r})"
set xlabel "{/Symbol r}"
plot [0.0:20.0][-1.6:1.6] \
    "run_coulfg.txt" index 0 using 1:3 with lines title "G_{1}(-2, {/Symbol r})", \
                  "" index 1 using 1:3 with lines title "G_{1}(0, {/Symbol r})", \
                  "" index 2 using 1:3 with lines title "G_{1}(2, {/Symbol r})", \
                  "" index 3 using 1:3 with lines title "G_{1}(4, {/Symbol r})", \
                  "" index 4 using 1:3 with lines title "G_{1}(6, {/Symbol r})", \
                  "" index 5 using 1:3 with lines title "G_{1}(8, {/Symbol r})", \
                  "" index 6 using 1:3 with lines title "G_{1}(10, {/Symbol r})"

set title "Coulomb Function F'_{{/Symbol l}=1}({/Symbol h}, {/Symbol r})"
set xlabel "{/Symbol r}"
plot [0.0:20.0][-1.6:1.6] \
    "run_coulfg.txt" index 0 using 1:2 with lines title "F'_{1}(-2, {/Symbol r})", \
                  "" index 1 using 1:2 with lines title "F'_{1}(0, {/Symbol r})", \
                  "" index 2 using 1:2 with lines title "F'_{1}(2, {/Symbol r})", \
                  "" index 3 using 1:2 with lines title "F'_{1}(4, {/Symbol r})", \
                  "" index 4 using 1:2 with lines title "F'_{1}(6, {/Symbol r})", \
                  "" index 5 using 1:2 with lines title "F'_{1}(8, {/Symbol r})", \
                  "" index 6 using 1:2 with lines title "F'_{1}(10, {/Symbol r})"

set title "Coulomb Function G'_{{/Symbol l}=1}({/Symbol h}, {/Symbol r})"
set xlabel "{/Symbol r}"
plot [0.0:20.0][-1.6:1.6] \
    "run_coulfg.txt" index 0 using 1:3 with lines title "G'_{1}(-2, {/Symbol r})", \
                  "" index 1 using 1:3 with lines title "G'_{1}(0, {/Symbol r})", \
                  "" index 2 using 1:3 with lines title "G'_{1}(2, {/Symbol r})", \
                  "" index 3 using 1:3 with lines title "G'_{1}(4, {/Symbol r})", \
                  "" index 4 using 1:3 with lines title "G'_{1}(6, {/Symbol r})", \
                  "" index 5 using 1:3 with lines title "G'_{1}(8, {/Symbol r})", \
                  "" index 6 using 1:3 with lines title "G'_{1}(10, {/Symbol r})"
