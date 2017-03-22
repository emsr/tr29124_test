
gnuplot

set termoption enhanced

set xzeroaxis
set yzeroaxis
set grid

set title "The Riemann Zeta Function and Friends"
set xlabel "s"
plot [-5.0:5.0][-5.0:5.0] "test_dirichlet_eta.txt" index 0 using 1:2 with lines title "{/Symbol z}(s)", \
                                         "" index 0 using 1:3 with lines title "{/Symbol h}(s)", \
                                         "" index 0 using 1:4 with lines title "{/Symbol b}(s)", \
                                         "" index 0 using 1:5 with lines title "{/Symbol l}(s)"
