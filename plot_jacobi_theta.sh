
gnuplot

set termoption enhanced

set hidden3d

set palette model RGB defined (-1 "black", -0.1 "blue", 0 "white", +0.1 "red", 1 "white")
set colorbox
set pm3d at bs corners2color geomean

# 
set title "Jacobi theta functions "
set xlabel "x"
plot [0.0:6.3][-1.1:1.1] \
                    "test_jacobi_theta.txt" index 0 using 1:2 with lines title "{/Symbol t}_1", \
                    "" index 0 using 1:3 with lines title "{/Symbol t}_2", \
                    "" index 0 using 1:4 with lines title "{/Symbol t}_3", \
                    "" index 0 using 1:5 with lines title "{/Symbol t}_4"

# 
set title "Jacobi theta functions "
set xlabel "x"
plot [0.0:6.3][-1.1:1.1] \
                    "test_jacobi_theta.txt" index 1 using 1:2 with lines title "{/Symbol t}_1", \
                    "" index 0 using 1:3 with lines title "{/Symbol t}_2", \
                    "" index 0 using 1:4 with lines title "{/Symbol t}_3", \
                    "" index 0 using 1:5 with lines title "{/Symbol t}_4"

# 
set title "Jacobi theta functions "
set xlabel "x"
plot [0.0:6.3][-1.1:1.1] \
                    "test_jacobi_theta.txt" index 2 using 1:2 with lines title "{/Symbol t}_1", \
                    "" index 0 using 1:3 with lines title "{/Symbol t}_2", \
                    "" index 0 using 1:4 with lines title "{/Symbol t}_3", \
                    "" index 0 using 1:5 with lines title "{/Symbol t}_4"

# 
set title "Jacobi theta functions "
set xlabel "x"
plot [0.0:6.3][-1.1:1.1] \
                    "test_jacobi_theta.txt" index 3 using 1:2 with lines title "", \
                    "" index 0 using 1:3 with lines title "", \
                    "" index 0 using 1:4 with lines title "", \
                    "" index 0 using 1:5 with lines title ""

# 
set title "Jacobi theta functions "
set xlabel "x"
plot [0.0:5.0][-1.1:1.1] \
                    "test_jacobi_theta.txt" index 4 using 1:2 with lines title "", \
                    "" index 0 using 1:3 with lines title "", \
                    "" index 0 using 1:4 with lines title "", \
                    "" index 0 using 1:5 with lines title ""

# 
set title "Jacobi theta functions "
set xlabel "x"
plot [0.0:1.0][-1.1:1.1] \
                    "test_jacobi_theta.txt" index 5 using 1:2 with lines title "", \
                    "" index 0 using 1:3 with lines title "", \
                    "" index 0 using 1:4 with lines title "", \
                    "" index 0 using 1:5 with lines title ""
