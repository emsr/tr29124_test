gnuplot
set termoption enhanced

# 
set title "Fibonacci function"
set xlabel "x"
plot [-10.0:10.0][-60:60] \
  "test_fibonacci.txt" index 0 using 1:2 with points title "F_n", \
                    "" index 1 using 1:2 with lines title "F_{/Symbol n}"

set term push
set term png
set output "fibonacci_1.png"
replot
set term pop

# 
set title "Fibonacci polynomials"
set xlabel "x"
plot [-5.0:5.0][-60:60] \
  "test_fibonacci.txt" index 2 using 1:2 with lines title "F_{0}(x)", \
                    "" index 3 using 1:2 with lines title "F_{1}(x)", \
                    "" index 4 using 1:2 with lines title "F_{2}(x)", \
                    "" index 5 using 1:2 with lines title "F_{3}(x)", \
                    "" index 6 using 1:2 with lines title "F_{4}(x)", \
                    "" index 7 using 1:2 with lines title "F_{5}(x)"

set term push
set term png
set output "fibonacci_2.png"
replot
set term pop

# 
set title "Lucas function"
set xlabel "x"
plot [-10.0:10.0][-80:120] \
  "test_fibonacci.txt" index 53 using 1:2 with points title "L_n", \
                    "" index 54 using 1:2 with lines title "L_{/Symbol n}"

set term push
set term png
set output "lucas_1.png"
replot
set term pop

# 
set title "Lucas polynomials"
set xlabel "x"
plot [-5.0:5.0][-60:60] \
  "test_fibonacci.txt" index 56 using 1:2 with lines title "L_{0}(x)", \
                    "" index 57 using 1:2 with lines title "L_{1}(x)", \
                    "" index 58 using 1:2 with lines title "L_{2}(x)", \
                    "" index 59 using 1:2 with lines title "L_{3}(x)", \
                    "" index 60 using 1:2 with lines title "L_{4}(x)", \
                    "" index 61 using 1:2 with lines title "L_{5}(x)"

set term push
set term png
set output "lucas_2.png"
replot
set term pop
