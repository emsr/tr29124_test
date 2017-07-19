gnuplot
set termoption enhanced

# 
set title "Fibonacci function"
set xlabel "x"
plot [-10.0:10.0][-60:60] \
  "test_fibonacci.txt" index 1 using 1:2 with lines title "F_{/Symbol n}"

set title "Lucas function"
set xlabel "x"
plot [-10.0:10.0][-80:120] \
  "test_fibonacci.txt" index 55 using 1:2 with lines title "L_{/Symbol n}"
