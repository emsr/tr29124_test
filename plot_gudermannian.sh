
gnuplot

load 'settings.gp'

set title "Gudermannian gd(x)"
set xlabel "x"
set ylabel "{/Symbol f}"
plot [-3.1:3.1][-1.5:1.5] \
    "test_gudermannian.txt" index 0 using 1:2 with lines ls 1 title "gd(x) - series", \
                         "" index 0 using 1:3 with lines ls 2 title "gd(x) - trig", \
                         "" index 0 using 1:5 with lines ls 3 title "invgd(x) - series", \
                         "" index 0 using 1:6 with lines ls 4 title "invgd(x) - trig"
