
gnuplot

set termoption enhanced
set xzeroaxis
set yzeroaxis
set grid

plot [0.0:pi][-0.0001:+0.0001] \
  "rational_fit.txt" index 0 using 1:4 with lines title "iter 1", \
                  "" index 1 using 1:4 with lines title "iter 2", \
                  "" index 2 using 1:4 with lines title "iter 3", \
                  "" index 3 using 1:4 with lines title "iter 4", \
                  "" index 4 using 1:4 with lines title "iter 5"

plot [0.0:pi][-0.01:+0.01] \
  "rational_fit_new.txt" index 0 using 1:4 with lines title "iter 1", \
                  "" index 1 using 1:4 with lines title "iter 2", \
                  "" index 2 using 1:4 with lines title "iter 3", \
                  "" index 3 using 1:4 with lines title "iter 4", \
                  "" index 4 using 1:4 with lines title "iter 5"

plot [0.0:pi][-0.0001:+0.0001] \
  "rational_fit_old.txt" index 0 using 1:4 with lines title "iter 1", \
                  "" index 1 using 1:4 with lines title "iter 2", \
                  "" index 2 using 1:4 with lines title "iter 3", \
                  "" index 3 using 1:4 with lines title "iter 4", \
                  "" index 4 using 1:4 with lines title "iter 5"
