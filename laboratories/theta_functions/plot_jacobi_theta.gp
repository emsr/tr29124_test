
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

set term push
set term png
set output "jacobi_theta_a.png"
replot
set term pop

# 
set title "Jacobi {/Symbol q}_1 function; q = 0.05, 0.5, 0.7, 0.9"
set xlabel "x"
plot [0.0:2.0][-6:6] \
  "test_jacobi_theta.txt" index 1 using 1:2 with lines title "{/Symbol q}_1(0.05,{/Symbol p}x)", \
                       "" index 1 using 1:3 with lines title "{/Symbol q}_1( 0.5,{/Symbol p}x)", \
                       "" index 1 using 1:4 with lines title "{/Symbol q}_1( 0.7,{/Symbol p}x)", \
                       "" index 1 using 1:5 with lines title "{/Symbol q}_1( 0.9,{/Symbol p}x)"

set term push
set term png
set output "jacobi_theta_b.png"
replot
set term pop

# 
set title "Jacobi {/Symbol q}_2 function; q = 0.05, 0.5, 0.7, 0.9"
set xlabel "x"
plot [0.0:2.0][-6:6] \
  "test_jacobi_theta.txt" index 2 using 1:2 with lines title "{/Symbol q}_2(0.05,{/Symbol p}x)", \
                       "" index 2 using 1:3 with lines title "{/Symbol q}_2( 0.5,{/Symbol p}x)", \
                       "" index 2 using 1:4 with lines title "{/Symbol q}_2( 0.7,{/Symbol p}x)", \
                       "" index 2 using 1:5 with lines title "{/Symbol q}_2( 0.9,{/Symbol p}x)"

set term push
set term png
set output "jacobi_theta_c.png"
replot
set term pop

# 
set title "Jacobi {/Symbol q}_3 function; q = 0.05, 0.5, 0.7, 0.9"
set xlabel "x"
plot [0.0:2.0][0:6] \
  "test_jacobi_theta.txt" index 3 using 1:2 with lines title "{/Symbol q}_3(0.05,{/Symbol p}x)", \
                       "" index 3 using 1:3 with lines title "{/Symbol q}_3( 0.5,{/Symbol p}x)", \
                       "" index 3 using 1:4 with lines title "{/Symbol q}_3( 0.7,{/Symbol p}x)", \
                       "" index 3 using 1:5 with lines title "{/Symbol q}_3( 0.9,{/Symbol p}x)"

set term push
set term png
set output "jacobi_theta_d.png"
replot
set term pop

# 
set title "Jacobi {/Symbol q}_4 function; q = 0.05, 0.5, 0.7, 0.9"
set xlabel "x"
plot [0.0:2.0][0:6] \
  "test_jacobi_theta.txt" index 4 using 1:2 with lines title "{/Symbol q}_4(0.05,{/Symbol p}x)", \
                       "" index 4 using 1:3 with lines title "{/Symbol q}_4( 0.5,{/Symbol p}x)", \
                       "" index 4 using 1:4 with lines title "{/Symbol q}_4( 0.7,{/Symbol p}x)", \
                       "" index 4 using 1:5 with lines title "{/Symbol q}_4( 0.9,{/Symbol p}x)"

set term push
set term png
set output "jacobi_theta_e.png"
replot
set term pop

# 
set title "Jacobi {/Symbol q}_1 function; x = 0, 0.4, 5, 10, 40"
set xlabel "q"
plot [0.0:1.0][-3:1.5] \
  "test_jacobi_theta.txt" index 5 using 1:2 with lines title "{/Symbol q}_1(q,  0)", \
                       "" index 5 using 1:3 with lines title "{/Symbol q}_1(q,0.4)", \
                       "" index 5 using 1:4 with lines title "{/Symbol q}_1(q,  5)", \
                       "" index 5 using 1:5 with lines title "{/Symbol q}_1(q, 10)", \
                       "" index 5 using 1:6 with lines title "{/Symbol q}_1(q, 40)"

set term push
set term png
set output "jacobi_theta_f.png"
replot
set term pop

# 
set title "Jacobi {/Symbol q}_2 function; x = 0, 0.4, 5, 10, 40"
set xlabel "q"
plot [0.0:1.0][-2:8] \
  "test_jacobi_theta.txt" index 6 using 1:2 with lines title "{/Symbol q}_2(q,  0)", \
                       "" index 6 using 1:3 with lines title "{/Symbol q}_2(q,0.4)", \
                       "" index 6 using 1:4 with lines title "{/Symbol q}_2(q,  5)", \
                       "" index 6 using 1:5 with lines title "{/Symbol q}_2(q, 10)", \
                       "" index 6 using 1:6 with lines title "{/Symbol q}_2(q, 40)"

set term push
set term png
set output "jacobi_theta_g.png"
replot
set term pop

# 
set title "Jacobi {/Symbol q}_3 function; x = 0, 0.4, 5, 10, 40"
set xlabel "q"
plot [0.0:1.0][-0.5:7.5] \
  "test_jacobi_theta.txt" index 7 using 1:2 with lines title "{/Symbol q}_3(q,  0)", \
                       "" index 7 using 1:3 with lines title "{/Symbol q}_3(q,0.4)", \
                       "" index 7 using 1:4 with lines title "{/Symbol q}_3(q,  5)", \
                       "" index 7 using 1:5 with lines title "{/Symbol q}_3(q, 10)", \
                       "" index 7 using 1:6 with lines title "{/Symbol q}_3(q, 40)"

set term push
set term png
set output "jacobi_theta_h.png"
replot
set term pop

# 
set title "Jacobi {/Symbol q}_4 function; x = 0, 0.4, 5, 10, 40"
set xlabel "q"
plot [0.0:1.0][-0.5:3] \
  "test_jacobi_theta.txt" index 8 using 1:2 with lines title "{/Symbol q}_4(q,  0)", \
                       "" index 8 using 1:3 with lines title "{/Symbol q}_4(q,0.4)", \
                       "" index 8 using 1:4 with lines title "{/Symbol q}_4(q,  5)", \
                       "" index 8 using 1:5 with lines title "{/Symbol q}_4(q, 10)", \
                       "" index 8 using 1:6 with lines title "{/Symbol q}_4(q, 40)"

set term push
set term png
set output "jacobi_theta_i.png"
replot
set term pop

# These are my old plots

# 
set title "Jacobi theta functions q = 0.0"
set xlabel "x"
plot [0.0:2.0][-2:2] \
  "test_jacobi_theta.txt" index 9 using 1:2 with lines title "{/Symbol q}_1(q,{/Symbol p}x)", \
                       "" index 9 using 1:3 with lines title "{/Symbol q}_2(q,{/Symbol p}x)", \
                       "" index 9 using 1:4 with lines title "{/Symbol q}_3(q,{/Symbol p}x)", \
                       "" index 9 using 1:5 with lines title "{/Symbol q}_4(q,{/Symbol p}x)"

set term push
set term png
set output "jacobi_theta_j.png"
replot
set term pop

# 
set title "Jacobi theta functions q = 0.1"
set xlabel "x"
plot [0.0:2.0][-2:2] \
  "test_jacobi_theta.txt" index 10 using 1:2 with lines title "{/Symbol q}_1(q,{/Symbol p}x)", \
                       "" index 10 using 1:3 with lines title "{/Symbol q}_2(q,{/Symbol p}x)", \
                       "" index 10 using 1:4 with lines title "{/Symbol q}_3(q,{/Symbol p}x)", \
                       "" index 10 using 1:5 with lines title "{/Symbol q}_4(q,{/Symbol p}x)"

set term push
set term png
set output "jacobi_theta_k.png"
replot
set term pop

# 
set title "Jacobi theta functions q = 0.2"
set xlabel "x"
plot [0.0:2.0][-2:2] \
  "test_jacobi_theta.txt" index 11 using 1:2 with lines title "{/Symbol q}_1(q,{/Symbol p}x)", \
                       "" index 11 using 1:3 with lines title "{/Symbol q}_2(q,{/Symbol p}x)", \
                       "" index 11 using 1:4 with lines title "{/Symbol q}_3(q,{/Symbol p}x)", \
                       "" index 11 using 1:5 with lines title "{/Symbol q}_4(q,{/Symbol p}x)"

set term push
set term png
set output "jacobi_theta_l.png"
replot
set term pop

# 
set title "Jacobi theta functions q = 0.3"
set xlabel "x"
plot [0.0:2.0][-2:2] \
  "test_jacobi_theta.txt" index 12 using 1:2 with lines title "{/Symbol q}_1(q,{/Symbol p}x)", \
                       "" index 12 using 1:3 with lines title "{/Symbol q}_2(q,{/Symbol p}x)", \
                       "" index 12 using 1:4 with lines title "{/Symbol q}_3(q,{/Symbol p}x)", \
                       "" index 12 using 1:5 with lines title "{/Symbol q}_4(q,{/Symbol p}x)"

set term push
set term png
set output "jacobi_theta_m.png"
replot
set term pop

# 
set title "Jacobi theta functions q = 0.4"
set xlabel "x"
plot [0.0:2.0][-2:2] \
  "test_jacobi_theta.txt" index 13 using 1:2 with lines title "{/Symbol q}_1(q,{/Symbol p}x)", \
                       "" index 13 using 1:3 with lines title "{/Symbol q}_2(q,{/Symbol p}x)", \
                       "" index 13 using 1:4 with lines title "{/Symbol q}_3(q,{/Symbol p}x)", \
                       "" index 13 using 1:5 with lines title "{/Symbol q}_4(q,{/Symbol p}x)"

set term push
set term png
set output "jacobi_theta_n.png"
replot
set term pop

# 
set title "Jacobi theta functions q = 0.5"
set xlabel "x"
plot [0.0:2.0][-3:3] \
  "test_jacobi_theta.txt" index 14 using 1:2 with lines title "{/Symbol q}_1(q,{/Symbol p}x)", \
                       "" index 14 using 1:3 with lines title "{/Symbol q}_2(q,{/Symbol p}x)", \
                       "" index 14 using 1:4 with lines title "{/Symbol q}_3(q,{/Symbol p}x)", \
                       "" index 14 using 1:5 with lines title "{/Symbol q}_4(q,{/Symbol p}x)"

set term push
set term png
set output "jacobi_theta_o.png"
replot
set term pop

# 
set title "Jacobi theta functions q = 0.6"
set xlabel "x"
plot [0.0:2.0][-3:3] \
  "test_jacobi_theta.txt" index 15 using 1:2 with lines title "{/Symbol q}_1(q,{/Symbol p}x)", \
                       "" index 15 using 1:3 with lines title "{/Symbol q}_2(q,{/Symbol p}x)", \
                       "" index 15 using 1:4 with lines title "{/Symbol q}_3(q,{/Symbol p}x)", \
                       "" index 15 using 1:5 with lines title "{/Symbol q}_4(q,{/Symbol p}x)"

set term push
set term png
set output "jacobi_theta_p.png"
replot
set term pop

# 
set title "Jacobi theta functions q = 0.7"
set xlabel "x"
plot [0.0:2.0][-4:4] \
  "test_jacobi_theta.txt" index 16 using 1:2 with lines title "{/Symbol q}_1(q,{/Symbol p}x)", \
                       "" index 16 using 1:3 with lines title "{/Symbol q}_2(q,{/Symbol p}x)", \
                       "" index 16 using 1:4 with lines title "{/Symbol q}_3(q,{/Symbol p}x)", \
                       "" index 16 using 1:5 with lines title "{/Symbol q}_4(q,{/Symbol p}x)"

set term push
set term png
set output "jacobi_theta_q.png"
replot
set term pop

# 
set title "Jacobi theta functions q = 0.8"
set xlabel "x"
plot [0.0:2.0][-4:4] \
  "test_jacobi_theta.txt" index 17 using 1:2 with lines title "{/Symbol q}_1(q,{/Symbol p}x)", \
                       "" index 17 using 1:3 with lines title "{/Symbol q}_2(q,{/Symbol p}x)", \
                       "" index 17 using 1:4 with lines title "{/Symbol q}_3(q,{/Symbol p}x)", \
                       "" index 17 using 1:5 with lines title "{/Symbol q}_4(q,{/Symbol p}x)"

set term push
set term png
set output "jacobi_theta_r.png"
replot
set term pop

# 
set title "Jacobi theta functions q = 0.9"
set xlabel "x"
plot [0.0:2.0][-6:6] \
  "test_jacobi_theta.txt" index 18 using 1:2 with lines title "{/Symbol q}_1(q,{/Symbol p}x)", \
                       "" index 18 using 1:3 with lines title "{/Symbol q}_2(q,{/Symbol p}x)", \
                       "" index 18 using 1:4 with lines title "{/Symbol q}_3(q,{/Symbol p}x)", \
                       "" index 18 using 1:5 with lines title "{/Symbol q}_4(q,{/Symbol p}x)"

set term push
set term png
set output "jacobi_theta_s.png"
replot
set term pop

# 
set title "Jacobi {/Symbol q}_1 functions q = 0.0 - 0.9"
set xlabel "x"
plot [0.0:2.0][-6:6] \
  "test_jacobi_theta.txt" index  9 using 1:2 with lines title "{/Symbol q}_1(0.0,{/Symbol p}x)", \
                       "" index 10 using 1:2 with lines title "{/Symbol q}_1(0.1,{/Symbol p}x)", \
                       "" index 11 using 1:2 with lines title "{/Symbol q}_1(0.2,{/Symbol p}x)", \
                       "" index 12 using 1:2 with lines title "{/Symbol q}_1(0.3,{/Symbol p}x)", \
                       "" index 13 using 1:2 with lines title "{/Symbol q}_1(0.4,{/Symbol p}x)", \
                       "" index 14 using 1:2 with lines title "{/Symbol q}_1(0.5,{/Symbol p}x)", \
                       "" index 15 using 1:2 with lines title "{/Symbol q}_1(0.6,{/Symbol p}x)", \
                       "" index 16 using 1:2 with lines title "{/Symbol q}_1(0.7,{/Symbol p}x)", \
                       "" index 17 using 1:2 with lines title "{/Symbol q}_1(0.8,{/Symbol p}x)", \
                       "" index 18 using 1:2 with lines title "{/Symbol q}_1(0.9,{/Symbol p}x)"

set term push
set term png
set output "jacobi_theta_1.png"
replot
set term pop

# 
set title "Jacobi {/Symbol q}_2 functions q = 0.0 - 0.9"
set xlabel "x"
plot [0.0:2.0][-6:6] \
  "test_jacobi_theta.txt" index  9 using 1:3 with lines title "{/Symbol q}_2(0.0,{/Symbol p}x)", \
                       "" index 10 using 1:3 with lines title "{/Symbol q}_2(0.1,{/Symbol p}x)", \
                       "" index 11 using 1:3 with lines title "{/Symbol q}_2(0.2,{/Symbol p}x)", \
                       "" index 12 using 1:3 with lines title "{/Symbol q}_2(0.3,{/Symbol p}x)", \
                       "" index 13 using 1:3 with lines title "{/Symbol q}_2(0.4,{/Symbol p}x)", \
                       "" index 14 using 1:3 with lines title "{/Symbol q}_2(0.5,{/Symbol p}x)", \
                       "" index 15 using 1:3 with lines title "{/Symbol q}_2(0.6,{/Symbol p}x)", \
                       "" index 16 using 1:3 with lines title "{/Symbol q}_2(0.7,{/Symbol p}x)", \
                       "" index 17 using 1:3 with lines title "{/Symbol q}_2(0.8,{/Symbol p}x)", \
                       "" index 18 using 1:3 with lines title "{/Symbol q}_2(0.9,{/Symbol p}x)"

set term push
set term png
set output "jacobi_theta_2.png"
replot
set term pop

# 
set title "Jacobi {/Symbol q}_3 functions q = 0.0 - 0.9"
set xlabel "x"
plot [0.0:2.0][-0.5:6] \
  "test_jacobi_theta.txt" index  9 using 1:4 with lines title "{/Symbol q}_3(0.0,{/Symbol p}x)", \
                       "" index 10 using 1:4 with lines title "{/Symbol q}_3(0.1,{/Symbol p}x)", \
                       "" index 11 using 1:4 with lines title "{/Symbol q}_3(0.2,{/Symbol p}x)", \
                       "" index 12 using 1:4 with lines title "{/Symbol q}_3(0.3,{/Symbol p}x)", \
                       "" index 13 using 1:4 with lines title "{/Symbol q}_3(0.4,{/Symbol p}x)", \
                       "" index 14 using 1:4 with lines title "{/Symbol q}_3(0.5,{/Symbol p}x)", \
                       "" index 15 using 1:4 with lines title "{/Symbol q}_3(0.6,{/Symbol p}x)", \
                       "" index 16 using 1:4 with lines title "{/Symbol q}_3(0.7,{/Symbol p}x)", \
                       "" index 17 using 1:4 with lines title "{/Symbol q}_3(0.8,{/Symbol p}x)", \
                       "" index 18 using 1:4 with lines title "{/Symbol q}_3(0.9,{/Symbol p}x)"

set term push
set term png
set output "jacobi_theta_3.png"
replot
set term pop

# 
set title "Jacobi {/Symbol q}_4 functions q = 0.0 - 0.9"
set xlabel "x"
plot [0.0:2.0][-0.5:6] \
  "test_jacobi_theta.txt" index  9 using 1:5 with lines title "{/Symbol q}_4(0.0,{/Symbol p}x)", \
                       "" index 10 using 1:5 with lines title "{/Symbol q}_4(0.1,{/Symbol p}x)", \
                       "" index 11 using 1:5 with lines title "{/Symbol q}_4(0.2,{/Symbol p}x)", \
                       "" index 12 using 1:5 with lines title "{/Symbol q}_4(0.3,{/Symbol p}x)", \
                       "" index 13 using 1:5 with lines title "{/Symbol q}_4(0.4,{/Symbol p}x)", \
                       "" index 14 using 1:5 with lines title "{/Symbol q}_4(0.5,{/Symbol p}x)", \
                       "" index 15 using 1:5 with lines title "{/Symbol q}_4(0.6,{/Symbol p}x)", \
                       "" index 16 using 1:5 with lines title "{/Symbol q}_4(0.7,{/Symbol p}x)", \
                       "" index 17 using 1:5 with lines title "{/Symbol q}_4(0.8,{/Symbol p}x)", \
                       "" index 18 using 1:5 with lines title "{/Symbol q}_4(0.9,{/Symbol p}x)"

set term push
set term png
set output "jacobi_theta_4.png"
replot
set term pop
