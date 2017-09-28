
gnuplot

set termoption enhanced

set xzeroaxis
set yzeroaxis
set grid

set title "Trigamma function {/Symbol y}^{(1)}(x)"
plot [0.0:20.0][0.0:1.0] "test_polygamma.txt" index 0 using 1:2 with lines title "{/Symbol y}^{(1)}(x)"

set title "Tetragamma function {/Symbol y}^{(2)}(x)"
plot [0.0:20.0][-0.5:0.0] "test_polygamma.txt" index 1 using 1:2 with lines title "{/Symbol y}^{(2)}(x)"
