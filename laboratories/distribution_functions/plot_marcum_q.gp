
gnuplot

set termoption enhanced

set xzeroaxis
set yzeroaxis
set grid

set title "Marcum Q function Q_M(a,b)"
set xlabel "b"
plot [0.0:10.0][-0.1:1.1] \
   "test_marcum_q.txt" index 0 using 1:2 with lines title "Q_M(1,b)"
