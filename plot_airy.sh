
/usr/local/bin/gnuplot --persist -e "plot [-20.0:5.0][-1.5:+1.5] \
 'plot/airy_float.txt' using 1:2 with pm3d title 'Ai', \
        	    '' using 1:3 with pm3d title 'Ai''', \
        	    '' using 1:4 with pm3d title 'Bi', \
        	    '' using 1:5 with pm3d title 'Bi'''"

/usr/local/bin/gnuplot --persist -e "plot [-20.0:5.0][-1.5:+1.5] \
 'plot/airy_double.txt' using 1:2 with pm3d title 'Ai', \
        	     '' using 1:3 with pm3d title 'Ai''', \
        	     '' using 1:4 with pm3d title 'Bi', \
        	     '' using 1:5 with pm3d title 'Bi'''"

/usr/local/bin/gnuplot --persist -e "plot [-20.0:5.0][-1.0:+1.0] \
 'plot/airy_long_double.txt' using 1:2 with pm3d title 'Ai', \
                	  '' using 1:4 with pm3d title 'Bi'"

/usr/local/bin/gnuplot --persist -e "plot [-20.0:5.0][-1.5:+1.5] \
 'plot/airy_long_double.txt' using 1:2 with pm3d title 'Ai', \
                	  '' using 1:3 with pm3d title 'Ai'''"

/usr/local/bin/gnuplot --persist -e "plot [-20.0:5.0][-1.5:+1.5] \
 'plot/airy_long_double.txt' using 1:4 with pm3d title 'Bi', \
                	  '' using 1:5 with pm3d title 'Bi'''"


/usr/local/bin/gnuplot --persist -e "plot [-20.0:5.0][-1.0:+1.0] \
 'plot/scorer_double.txt' using 1:2 with pm3d title 'Gi', \
                	  '' using 1:4 with pm3d title 'Hi'"

/usr/local/bin/gnuplot --persist -e "plot [-20.0:5.0][-1.5:+1.5] \
 'plot/scorer_double.txt' using 1:2 with pm3d title 'Gi', \
                	  '' using 1:3 with pm3d title 'Gi'''"

/usr/local/bin/gnuplot --persist -e "plot [-20.0:5.0][-1.5:+1.5] \
 'plot/scorer_double.txt' using 1:4 with pm3d title 'Hi', \
                	  '' using 1:5 with pm3d title 'Hi'''"


/usr/local/bin/gnuplot --persist -e "splot [-20.0:5.0][-5.0:5.0][-1.5:+1.5] \
 'plot/airy_complex_double.txt' title 'Ai'"


/usr/local/bin/gnuplot
set xzeroaxis
set yzeroaxis
set grid

set palette model RGB defined (0.0 "blue", 0.3 "green", 0.7 "yellow", 1 "red")
set colorbox
set pm3d at sb corners2color geomean

plot [-5.0:2.0][-1.5:+1.5] "plot/fgh_double.txt" using 1:2 with pm3d title "fai", "" using 1:3 with pm3d title "gai", "" using 1:4 with pm3d title "hai"
plot [-5.0:2.0][-2.5:+2.5] "plot/fgh_double.txt" using 1:5 with pm3d title "fai'", "" using 1:6 with pm3d title "gai'", "" using 1:7 with pm3d title "hai'"
plot [-5.0:2.0][-2.5:+2.5] "plot/fgh_double.txt" using 1:($2+$3+$4) with pm3d title "Hi"
plot [-5.0:2.0][-2.5:+2.5] "plot/fgh_double.txt" using 1:($5+$6+$7) with pm3d title "Hi'"

plot [-20.0:5.0][-1.0:+1.0] "plot/scorer_double.txt" using 1:2 with pm3d title "Gi", "" using 1:4 with pm3d title "Hi"
plot [-20.0:5.0][-1.5:+1.5] "plot/scorer_double.txt" using 1:2 with pm3d title "Gi", "" using 1:3 with pm3d title "Gi'"
plot [-20.0:5.0][-1.5:+1.5] "plot/scorer_double.txt" using 1:4 with pm3d title "Hi", "" using 1:5 with pm3d title "Hi'"

plot [-20.0:5.0][-1.0:+1.0] "plot/scorer_float.txt" using 1:2 with pm3d title "Gi", "" using 1:3 with pm3d title "Gi'", "" using 1:4 with pm3d title "Hi", "" using 1:5 with pm3d title "Hi'"
plot [-20.0:5.0][-1.0:+1.0] "plot/scorer_double.txt" using 1:2 with pm3d title "Gi", "" using 1:3 with pm3d title "Gi'", "" using 1:4 with pm3d title "Hi", "" using 1:5 with pm3d title "Hi'"
plot [-20.0:5.0][-1.0:+1.0] "plot/scorer_long_double.txt" using 1:2 with pm3d title "Gi", "" using 1:3 with pm3d title "Gi'", "" using 1:4 with pm3d title "Hi", "" using 1:5 with pm3d title "Hi'"

plot [-20.0:5.0][-1.0:+1.0] "plot/airy_float.txt" using 1:2 with pm3d title "Ai(x)", "" using 1:4 with pm3d title "Bi(x)"
plot [-20.0:5.0][-1.5:+1.5] "plot/airy_float.txt" using 1:2 with pm3d title "Ai(x)", "" using 1:3 with pm3d title "Ai'(x)"
plot [-20.0:5.0][-1.5:+1.5] "plot/airy_float.txt" using 1:4 with pm3d title "Bi(x)", "" using 1:5 with pm3d title "Bi'(x)"

plot [-20.0:5.0][-1.0:+1.0] "plot/airy_double.txt" using 1:2 with pm3d title "Ai(x)", "" using 1:4 with pm3d title "Bi(x)"
plot [-20.0:5.0][-1.5:+1.5] "plot/airy_double.txt" using 1:2 with pm3d title "Ai(x)", "" using 1:3 with pm3d title "Ai'(x)"
plot [-20.0:5.0][-1.5:+1.5] "plot/airy_double.txt" using 1:4 with pm3d title "Bi(x)", "" using 1:5 with pm3d title "Bi'(x)"

plot [-20.0:5.0][-1.0:+1.0] "plot/airy_long_double.txt" using 1:2 with pm3d title "Ai(x)", "" using 1:4 with pm3d title "Bi(x)"
plot [-20.0:5.0][-1.5:+1.5] "plot/airy_long_double.txt" using 1:2 with pm3d title "Ai(x)", "" using 1:3 with pm3d title "Ai'(x)"
plot [-20.0:5.0][-1.5:+1.5] "plot/airy_long_double.txt" using 1:4 with pm3d title "Bi(x)", "" using 1:5 with pm3d title "Bi'(x)"

plot [0.0:30.0][-0.5:+6.0] "plot/struve_float.txt" using 1:2 with pm3d title "H_0", "" using 1:4 with pm3d title "H_1", "" using 1:6 with pm3d title "H_2", "" using 1:8 with pm3d title "H_3"
plot [0.0:5.0][-0.5:+4.5] "plot/struve_float.txt" using 1:3 with pm3d title "L_0", "" using 1:5 with pm3d title "L_1", "" using 1:7 with pm3d title "L_2", "" using 1:9 with pm3d title "L_3"

plot [0.0:30.0][-0.5:+6.0] "plot/struve_double.txt" using 1:2 with pm3d title "H_0", "" using 1:4 with pm3d title "H_1", "" using 1:6 with pm3d title "H_2", "" using 1:8 with pm3d title "H_3"
plot [0.0:5.0][-0.5:+4.5] "plot/struve_double.txt" using 1:3 with pm3d title "L_0", "" using 1:5 with pm3d title "L_1", "" using 1:7 with pm3d title "L_2", "" using 1:9 with pm3d title "L_3"

plot [0.0:30.0][-0.5:+6.0] "plot/struve_long_double.txt" using 1:2 with pm3d title "H_0", "" using 1:4 with pm3d title "H_1", "" using 1:6 with pm3d title "H_2", "" using 1:8 with pm3d title "H_3"
plot [0.0:5.0][-0.5:+4.5] "plot/struve_long_double.txt" using 1:3 with pm3d title "L_0", "" using 1:5 with pm3d title "L_1", "" using 1:7 with pm3d title "L_2", "" using 1:9 with pm3d title "L_3"

plot [0.0:40.0][-1.5:1.5] "plot/kelvin_float.txt" using 1:2 with pm3d title "ber(x)", "" using 1:3 with pm3d title "bei(x)", "" using 1:4 with pm3d title "ker(x)", "" using 1:5 with pm3d title "kei(x)"
plot [0.0:40.0][-1.5:1.5] "plot/kelvin_double.txt" using 1:2 with pm3d title "ber(x)", "" using 1:3 with pm3d title "bei(x)", "" using 1:4 with pm3d title "ker(x)", "" using 1:5 with pm3d title "kei(x)"
plot [0.0:40.0][-1.5:1.5] "plot/kelvin_long_double.txt" using 1:2 with pm3d title "ber(x)", "" using 1:3 with pm3d title "bei(x)", "" using 1:4 with pm3d title "ker(x)", "" using 1:5 with pm3d title "kei(x)"

plot [0.0:40.0][-1.5:1.5] "plot/kelvin_order_float.txt" using 1:2 with pm3d title "ber(x)", "" using 1:3 with pm3d title "bei(x)", "" using 1:4 with pm3d title "ker(x)", "" using 1:5 with pm3d title "kei(x)"
plot [0.0:40.0][-1.5:1.5] "plot/kelvin_order_double.txt" using 1:2 with pm3d title "ber(x)", "" using 1:3 with pm3d title "bei(x)", "" using 1:4 with pm3d title "ker(x)", "" using 1:5 with pm3d title "kei(x)"
plot [0.0:40.0][-1.5:1.5] "plot/kelvin_order_long_double.txt" using 1:2 with pm3d title "ber(x)", "" using 1:3 with pm3d title "bei(x)", "" using 1:4 with pm3d title "ker(x)", "" using 1:5 with pm3d title "kei(x)"

/usr/local/bin/gnuplot
set xzeroaxis
set yzeroaxis
set grid

set palette model RGB defined (0.0 "blue", 0.3 "green", 0.7 "yellow", 1 "red")
set colorbox
set pm3d at sb corners2color geomean

/usr/local/bin/gnuplot

set hidden3d

splot [-20.0:5.0][-5.0:5.0][-1.5:+1.5] "plot/airy_complex_double.txt" title "Ai"


splot [-5:5][-5:5][0:12] "plot/airy_complex_float.txt" index 0 with pm3d title "|Ai(z)|^{1/6}"
splot [-20:5][-5:5][-20:20] "plot/airy_complex_float.txt" index 1 with pm3d title "Re[Ai(z)]"
splot [-20:5][-5:5][-20:20] "plot/airy_complex_float.txt" index 2 with pm3d title "Im[Ai(z)]"
splot [-20:5][-5:5][-3.2:+3.2] "plot/airy_complex_float.txt" index 3 with pm3d title "Arg[Ai(z)]"
splot [-5:5][-5:5][0:12] "plot/airy_complex_float.txt" index 4 with pm3d title "|Bi(z)|^{1/6}"
splot [-20:5][-5:5][-20:20] "plot/airy_complex_float.txt" index 5 with pm3d title "Re[Bi(z)]"
splot [-20:5][-5:5][-20:20] "plot/airy_complex_float.txt" index 6 with pm3d title "Im[Bi(z)]"
splot [-20:5][-5:5][-3.2:+3.2] "plot/airy_complex_float.txt" index 7 with pm3d title "Arg[Bi(z)]"
splot [-20:5][-5:5][-2:2] "plot/airy_complex_float.txt" index 8 with pm3d title "Wronski"

splot [-5:5][-5:5][0:12] "plot/airy_complex_double.txt" index 0 with pm3d title "|Ai(z)|^{1/6}"
splot [-20:5][-5:5][-20:20] "plot/airy_complex_double.txt" index 1 with pm3d title "Re[Ai(z)]"
splot [-20:5][-5:5][-20:20] "plot/airy_complex_double.txt" index 2 with pm3d title "Im[Ai(z)]"
splot [-20:5][-5:5][-3.2:+3.2] "plot/airy_complex_double.txt" index 3 with pm3d title "Arg[Ai(z)]"
splot [-5:5][-5:5][0:12] "plot/airy_complex_double.txt" index 4 with pm3d title "|Bi(z)|^{1/6}"
splot [-20:5][-5:5][-20:20] "plot/airy_complex_double.txt" index 5 with pm3d title "Re[Bi(z)]"
splot [-20:5][-5:5][-20:20] "plot/airy_complex_double.txt" index 6 with pm3d title "Im[Bi(z)]"
splot [-20:5][-5:5][-3.2:+3.2] "plot/airy_complex_double.txt" index 7 with pm3d title "Arg[Bi(z)]"
splot [-20:5][-5:5][-2:2] "plot/airy_complex_double.txt" index 8 with pm3d title "Wronski"

splot [-5:5][-5:5][0:12] "plot/airy_complex_long_double.txt" index 0 with pm3d title "|Ai(z)|^{1/6}"
splot [-20:5][-5:5][-20:20] "plot/airy_complex_long_double.txt" index 1 with pm3d title "Re[Ai(z)]"
splot [-20:5][-5:5][-20:20] "plot/airy_complex_long_double.txt" index 2 with pm3d title "Im[Ai(z)]"
splot [-20:5][-5:5][-3.2:+3.2] "plot/airy_complex_long_double.txt" index 3 with pm3d title "Arg[Ai(z)]"
splot [-5:5][-5:5][0:12] "plot/airy_complex_long_double.txt" index 4 with pm3d title "|Bi(z)|^{1/6}"
splot [-20:5][-5:5][-20:20] "plot/airy_complex_long_double.txt" index 5 with pm3d title "Re[Bi(z)]"
splot [-20:5][-5:5][-20:20] "plot/airy_complex_long_double.txt" index 6 with pm3d title "Im[Bi(z)]"
splot [-20:5][-5:5][-3.2:+3.2] "plot/airy_complex_long_double.txt" index 7 with pm3d title "Arg[Bi(z)]"
splot [-20:5][-5:5][-2:2] "plot/airy_complex_long_double.txt" index 8 with pm3d title "Wronski"

