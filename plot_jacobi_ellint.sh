
gnuplot

set termoption enhanced

set xzeroaxis
set yzeroaxis
set grid

# Figure 1
set title "Jacobi elliptic sine amplitude function sn(k,u)"
set xlabel "u"
plot [-10.0:10.0][-1.1:1.1] \
                    "test_jacobi_ellint.txt" index 0 using 1:2 with lines title "k=0.0", \
                    "" index 0 using 1:3 with lines title "k=0.5", \
                    "" index 0 using 1:4 with lines title "k=0.75", \
                    "" index 0 using 1:5 with lines title "k=0.95", \
                    "" index 0 using 1:6 with lines title "k=1.0"

# Figure 2
set title "Jacobi elliptic cosine amplitude function cn(k,u)"
set xlabel "u"
plot [-10.0:10.0][-1.1:1.1] \
                    "test_jacobi_ellint.txt" index 1 using 1:2 with lines title "k=0.0", \
                    "" index 1 using 1:3 with lines title "k=0.5", \
                    "" index 1 using 1:4 with lines title "k=0.75", \
                    "" index 1 using 1:5 with lines title "k=0.95", \
                    "" index 1 using 1:6 with lines title "k=1.0"

# Figure 2
set title "Jacobi elliptic delta amplitude function dn(k,u)"
set xlabel "u"
plot [-10.0:10.0][-0.1:1.1] \
                    "test_jacobi_ellint.txt" index 2 using 1:2 with lines title "k=0.0", \
                    "" index 2 using 1:3 with lines title "k=0.5", \
                    "" index 2 using 1:4 with lines title "k=0.75", \
                    "" index 2 using 1:5 with lines title "k=0.95", \
                    "" index 2 using 1:6 with lines title "k=1.0"
