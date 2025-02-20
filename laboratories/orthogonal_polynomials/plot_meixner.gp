# http://pool-serv1.mathematik.uni-kassel.de/caop/CAOP/plot_meixner.html

gnuplot

load 'setings.gp'

set title "Meixner Polynomial M_n(x; {/Symbol b}, c)"
set xlabel "x"
plot [0.0:16.0][-13.0:7.5] \
    "test_meixner.txt" index 5 using 1:2 with lines title "M_5(x; {/Symbol b}, c)"

set title "Meixner Polynomial M_n(x; {/Symbol b}, c)"
set xlabel "x"
plot [0.0:16.0][-20.0:45.0] \
    "test_meixner.txt" index 0 using 1:2 with lines title "M_0(x; {/Symbol b}, c)", \
                    "" index 1 using 1:2 with lines title "M_1(x; {/Symbol b}, c)", \
                    "" index 2 using 1:2 with lines title "M_2(x; {/Symbol b}, c)", \
                    "" index 3 using 1:2 with lines title "M_3(x; {/Symbol b}, c)", \
                    "" index 4 using 1:2 with lines title "M_4(x; {/Symbol b}, c)", \
                    "" index 5 using 1:2 with lines title "M_5(x; {/Symbol b}, c)", \
                    "" index 6 using 1:2 with lines title "M_6(x; {/Symbol b}, c)", \
                    "" index 7 using 1:2 with lines title "M_7(x; {/Symbol b}, c)", \
                    "" index 8 using 1:2 with lines title "M_8(x; {/Symbol b}, c)", \
                    "" index 9 using 1:2 with lines title "M_9(x; {/Symbol b}, c)", \
                    "" index 10 using 1:2 with lines title "M_{10}(x; {/Symbol b}, c)"

set term push
set term png
set output "meixner.png"
replot
set term pop
