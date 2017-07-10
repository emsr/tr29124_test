
gnuplot

set termoption enhanced

#set hidden3d

#set palette model RGB defined (-1 "black", -0.1 "blue", 0 "white", +0.1 "red", 1 "white")
#set colorbox
#set pm3d at bs corners2color geomean

# These plots should follow DLMF.

# 
set title "Jacobi theta functions q = 0.15"
set xlabel "x"
plot [0.0:2.0][-2:2] \
  "test_jacobi_theta.txt" index 0 using 1:2 with lines title "{/Symbol q}_1(q,{/Symbol p}x)", \
                       "" index 0 using 1:3 with lines title "{/Symbol q}_2(q,{/Symbol p}x)", \
                       "" index 0 using 1:4 with lines title "{/Symbol q}_3(q,{/Symbol p}x)", \
                       "" index 0 using 1:5 with lines title "{/Symbol q}_4(q,{/Symbol p}x)"

# 
set title "Jacobi {/Symbol q}_1 function; q = 0.05, 0.5, 0.7, 0.9"
set xlabel "x"
plot [0.0:2.0][-6:6] \
  "test_jacobi_theta.txt" index 1 using 1:2 with lines title "{/Symbol q}_1(0.05,{/Symbol p}x)", \
                       "" index 1 using 1:3 with lines title "{/Symbol q}_1( 0.5,{/Symbol p}x)", \
                       "" index 1 using 1:4 with lines title "{/Symbol q}_1( 0.7,{/Symbol p}x)", \
                       "" index 1 using 1:5 with lines title "{/Symbol q}_1( 0.9,{/Symbol p}x)"

# 
set title "Jacobi {/Symbol q}_2 function; q = 0.05, 0.5, 0.7, 0.9"
set xlabel "x"
plot [0.0:2.0][-6:6] \
  "test_jacobi_theta.txt" index 2 using 1:2 with lines title "{/Symbol q}_2(0.05,{/Symbol p}x)", \
                       "" index 2 using 1:3 with lines title "{/Symbol q}_2( 0.5,{/Symbol p}x)", \
                       "" index 2 using 1:4 with lines title "{/Symbol q}_2( 0.7,{/Symbol p}x)", \
                       "" index 2 using 1:5 with lines title "{/Symbol q}_2( 0.9,{/Symbol p}x)"

# 
set title "Jacobi {/Symbol q}_3 function; q = 0.05, 0.5, 0.7, 0.9"
set xlabel "x"
plot [0.0:2.0][0:6] \
  "test_jacobi_theta.txt" index 3 using 1:2 with lines title "{/Symbol q}_3(0.05,{/Symbol p}x)", \
                       "" index 3 using 1:3 with lines title "{/Symbol q}_3( 0.5,{/Symbol p}x)", \
                       "" index 3 using 1:4 with lines title "{/Symbol q}_3( 0.7,{/Symbol p}x)", \
                       "" index 3 using 1:5 with lines title "{/Symbol q}_3( 0.9,{/Symbol p}x)"

# 
set title "Jacobi {/Symbol q}_4 function; q = 0.05, 0.5, 0.7, 0.9"
set xlabel "x"
plot [0.0:2.0][0:6] \
  "test_jacobi_theta.txt" index 4 using 1:2 with lines title "{/Symbol q}_4(0.05,{/Symbol p}x)", \
                       "" index 4 using 1:3 with lines title "{/Symbol q}_4( 0.5,{/Symbol p}x)", \
                       "" index 4 using 1:4 with lines title "{/Symbol q}_4( 0.7,{/Symbol p}x)", \
                       "" index 4 using 1:5 with lines title "{/Symbol q}_4( 0.9,{/Symbol p}x)"

# 
set title "Jacobi {/Symbol q}_1 function; x = 0, 0.4, 5, 10, 40"
set xlabel "x"
plot [0.0:1.0][-3:1] \
  "test_jacobi_theta.txt" index 5 using 1:2 with lines title "{/Symbol q}_1(q,  0)", \
                       "" index 5 using 1:3 with lines title "{/Symbol q}_1(q,0.4)", \
                       "" index 5 using 1:4 with lines title "{/Symbol q}_1(q,  5)", \
                       "" index 5 using 1:5 with lines title "{/Symbol q}_1(q, 10)", \
                       "" index 5 using 1:6 with lines title "{/Symbol q}_1(q, 40)"

# 
set title "Jacobi {/Symbol q}_2 function; x = 0, 0.4, 5, 10, 40"
set xlabel "x"
plot [0.0:1.0][-2:3] \
  "test_jacobi_theta.txt" index 6 using 1:2 with lines title "{/Symbol q}_2(q,  0)", \
                       "" index 6 using 1:3 with lines title "{/Symbol q}_2(q,0.4)", \
                       "" index 6 using 1:4 with lines title "{/Symbol q}_2(q,  5)", \
                       "" index 6 using 1:5 with lines title "{/Symbol q}_2(q, 10)", \
                       "" index 6 using 1:6 with lines title "{/Symbol q}_2(q, 40)"

# 
set title "Jacobi {/Symbol q}_3 function; x = 0, 0.4, 5, 10, 40"
set xlabel "x"
plot [0.0:1.0][0:3] \
  "test_jacobi_theta.txt" index 7 using 1:2 with lines title "{/Symbol q}_3(q,  0)", \
                       "" index 7 using 1:3 with lines title "{/Symbol q}_3(q,0.4)", \
                       "" index 7 using 1:4 with lines title "{/Symbol q}_3(q,  5)", \
                       "" index 7 using 1:5 with lines title "{/Symbol q}_3(q, 10)", \
                       "" index 7 using 1:6 with lines title "{/Symbol q}_3(q, 40)"

# 
set title "Jacobi {/Symbol q}_4 function; x = 0, 0.4, 5, 10, 40"
set xlabel "x"
plot [0.0:1.0][0:3] \
  "test_jacobi_theta.txt" index 8 using 1:2 with lines title "{/Symbol q}_4(q,  0)", \
                       "" index 8 using 1:3 with lines title "{/Symbol q}_4(q,0.4)", \
                       "" index 8 using 1:4 with lines title "{/Symbol q}_4(q,  5)", \
                       "" index 8 using 1:5 with lines title "{/Symbol q}_4(q, 10)", \
                       "" index 8 using 1:6 with lines title "{/Symbol q}_4(q, 40)"
