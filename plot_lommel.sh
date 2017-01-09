
gnuplot
set termoption enhanced
set xzeroaxis
set yzeroaxis
set grid

plot [0.0:20.0][-60.0:80.0] "test_lommel.txt" index 0 using 1:2 with lines title "s_{6/5,1/5}(x)"

plot [0.0:20.0][-60.0:80.0] "test_lommel.txt" index 1 using 1:2 with lines title "S_{6/5,1/5}(x)"


