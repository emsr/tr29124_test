# http://pool-serv1.mathematik.uni-kassel.de/caop/CAOP/plot_racah.html

gnuplot

load 'settings.gp'

set title "Racah polynomial R_n({/Symbol l}(x); a, b, c, d)"
set xlabel "x"
plot [0.0:180.0][-0.1:1.1] \
    "test_racah.txt" index 5 using 1:2 with lines ls 1 title "R_5({/Symbol l}(x); a, b, c, d)"

set title "Racah polynomial R_n({/Symbol l}(x); a, b, c, d)"
set xlabel "x"
plot [0.0:200.0][-20.0:20.0] \
    "test_racah.txt" index 0 using 1:2 with lines ls 1 title "R_0({/Symbol l}(x); a, b, c, d)", \
                  "" index 1 using 1:2 with lines ls 2 title "R_1({/Symbol l}(x); a, b, c, d)", \
                  "" index 2 using 1:2 with lines ls 3 title "R_2({/Symbol l}(x); a, b, c, d)", \
                  "" index 3 using 1:2 with lines ls 4 title "R_3({/Symbol l}(x); a, b, c, d)", \
                  "" index 4 using 1:2 with lines ls 5 title "R_4({/Symbol l}(x); a, b, c, d)", \
                  "" index 5 using 1:2 with lines ls 6 title "R_5({/Symbol l}(x); a, b, c, d)", \
                  "" index 6 using 1:2 with lines ls 7 title "R_6({/Symbol l}(x); a, b, c, d)", \
                  "" index 7 using 1:2 with lines ls 8 title "R_7({/Symbol l}(x); a, b, c, d)", \
                  "" index 8 using 1:2 with lines ls 9 title "R_8({/Symbol l}(x); a, b, c, d)", \
                  "" index 9 using 1:2 with lines ls 10 title "R_9({/Symbol l}(x); a, b, c, d)", \
                  "" index 10 using 1:2 with lines ls 11 title "R_{10}({/Symbol l}(x); a, b, c, d)"
