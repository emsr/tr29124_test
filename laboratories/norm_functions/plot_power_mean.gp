
gnuplot

set termoption enhanced

set xzeroaxis
set yzeroaxis
set grid

set title "Power mean ..."
plot [0.0:10.0][1.5:6.7] \
                    "test_power_mean.txt" index 0 using 1:2 with lines title "powmean"

set term push
set term png
set output "power_mean.png"
replot
set term pop
