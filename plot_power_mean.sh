
gnuplot

set termoption enhanced

set xzeroaxis
set yzeroaxis
set grid

plot [0.0:10.0][1.5:6.7] \
                    "test_power_mean.txt" index 0 using 1:2 with lines title "powmean"
