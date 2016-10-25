
/usr/local/bin/gnuplot

set hidden3d

set palette model RGB defined (0.0 "blue", 0.3 "green", 0.7 "yellow", 1 "red")
set colorbox
set pm3d at sb corners2color geomean

splot [-5:5][-5:5][0:12] "plot/airy_complex_double_new.txt" index 0 with lines title "|Ai(z)|^{1/6}"
splot [-20:5][-5:5][-20:20] "plot/airy_complex_double_new.txt" index 1 with lines title "Re[Ai(z)]"
splot [-20:5][-5:5][-20:20] "plot/airy_complex_double_new.txt" index 2 with lines title "Im[Ai(z)]"
splot [-20:5][-5:5][-3.2:+3.2] "plot/airy_complex_double_new.txt" index 3 with lines title "Arg[Ai(z)]"
splot [-5:5][-5:5][0:12] "plot/airy_complex_double_new.txt" index 4 with lines title "|Bi(z)|^{1/6}"
splot [-20:5][-5:5][-2:2] "plot/airy_complex_double_new.txt" index 5 with lines title "Re[Bi(z)]"
splot [-20:5][-5:5][-2:2] "plot/airy_complex_double_new.txt" index 6 with lines title "Im[Bi(z)]"
splot [-20:5][-5:5][-3.2:+3.2] "plot/airy_complex_double_new.txt" index 7 with lines title "Arg[Bi(z)]"
#splot [-20:5][-5:5][-2:2] "plot/airy_complex_double_new.txt" index 8 with lines title "Wronski"
