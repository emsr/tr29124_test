
gnuplot

load '../plotting_tools/settings.gp'

set title "Jacobi polynomial P_n^{({/Symbol a},{/Symbol b})}(x)"
set xlabel "x"
plot [-2.0:2.0][-100.0:100.0] \
    "test_jacobi_neg_parm.txt" index 0 using 1:2 with lines ls 1 title "P_n^{({/Symbol a},{/Symbol b-})}(x)", \
                            "" index 0 using 1:3 with lines ls 2 title "P_n^{({/Symbol a},{/Symbol b})}(x)", \
                            "" index 0 using 1:4 with lines ls 3 title "P_n^{({/Symbol a},{/Symbol b+})}(x)"

set title "Jacobi polynomial P_n^{({/Symbol a},{/Symbol b})}(x)"
set xlabel "x"
plot [-2.0:2.0][-100.0:100.0] \
    "test_jacobi_neg_parm.txt" index 1 using 1:2 with lines ls 1 title "P_n^{({/Symbol a},{/Symbol b-})}(x)", \
                            "" index 1 using 1:3 with lines ls 2 title "P_n^{({/Symbol a},{/Symbol b})}(x)", \
                            "" index 1 using 1:4 with lines ls 3 title "P_n^{({/Symbol a},{/Symbol b+})}(x)"

set title "Jacobi polynomial P_n^{({/Symbol a},{/Symbol b})}(x)"
set xlabel "x"
plot [-2.0:2.0][-100.0:100.0] \
    "test_jacobi_neg_parm.txt" index 2 using 1:2 with lines ls 1 title "P_n^{({/Symbol a},{/Symbol b-})}(x)", \
                            "" index 2 using 1:3 with lines ls 2 title "P_n^{({/Symbol a},{/Symbol b})}(x)", \
                            "" index 2 using 1:4 with lines ls 3 title "P_n^{({/Symbol a},{/Symbol b+})}(x)"


plot "jacobi_roots.gp" index 0  using 1:2 with points
