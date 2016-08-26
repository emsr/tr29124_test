
/usr/local/bin/gnuplot --persist -e "plot [-20.0:5.0][-1.5:+1.5] \
 'plot/airy_float.txt' using 1:2 with lines title 'Ai', \
        	    '' using 1:3 with lines title 'Ai''', \
        	    '' using 1:4 with lines title 'Bi', \
        	    '' using 1:5 with lines title 'Bi'''"

/usr/local/bin/gnuplot --persist -e "plot [-20.0:5.0][-1.5:+1.5] \
 'plot/airy_double.txt' using 1:2 with lines title 'Ai', \
        	     '' using 1:3 with lines title 'Ai''', \
        	     '' using 1:4 with lines title 'Bi', \
        	     '' using 1:5 with lines title 'Bi'''"

/usr/local/bin/gnuplot --persist -e "plot [-20.0:5.0][-1.0:+1.0] \
 'plot/airy_long_double.txt' using 1:2 with lines title 'Ai', \
                	  '' using 1:4 with lines title 'Bi'"

/usr/local/bin/gnuplot --persist -e "plot [-20.0:5.0][-1.5:+1.5] \
 'plot/airy_long_double.txt' using 1:2 with lines title 'Ai', \
                	  '' using 1:3 with lines title 'Ai'''"

/usr/local/bin/gnuplot --persist -e "plot [-20.0:5.0][-1.5:+1.5] \
 'plot/airy_long_double.txt' using 1:4 with lines title 'Bi', \
                	  '' using 1:5 with lines title 'Bi'''"


/usr/local/bin/gnuplot --persist -e "plot [-20.0:5.0][-1.0:+1.0] \
 'plot/scorer_double.txt' using 1:2 with lines title 'Gi', \
                	  '' using 1:4 with lines title 'Hi'"

/usr/local/bin/gnuplot --persist -e "plot [-20.0:5.0][-1.5:+1.5] \
 'plot/scorer_double.txt' using 1:2 with lines title 'Gi', \
                	  '' using 1:3 with lines title 'Gi'''"

/usr/local/bin/gnuplot --persist -e "plot [-20.0:5.0][-1.5:+1.5] \
 'plot/scorer_double.txt' using 1:4 with lines title 'Hi', \
                	  '' using 1:5 with lines title 'Hi'''"


/usr/local/bin/gnuplot --persist -e "splot [-20.0:5.0][-5.0:5.0][-1.5:+1.5] \
 'plot/airy_complex_double.txt' title 'Ai'"


/usr/local/bin/gnuplot
set xzeroaxis
set yzeroaxis
set grid

plot [-5.0:2.0][-1.5:+1.5] "plot/fgh_double.txt" using 1:2 with lines title "fai", "" using 1:3 with lines title "gai", "" using 1:4 with lines title "hai"
plot [-5.0:2.0][-2.5:+2.5] "plot/fgh_double.txt" using 1:5 with lines title "fai'", "" using 1:6 with lines title "gai'", "" using 1:7 with lines title "hai'"
plot [-5.0:2.0][-2.5:+2.5] "plot/fgh_double.txt" using 1:($2+$3+$4) with lines title "Hi"
plot [-5.0:2.0][-2.5:+2.5] "plot/fgh_double.txt" using 1:($5+$6+$7) with lines title "Hi'"

plot [-20.0:5.0][-1.0:+1.0] "plot/scorer_double.txt" using 1:2 with lines title "Gi", "" using 1:4 with lines title "Hi"
plot [-20.0:5.0][-1.5:+1.5] "plot/scorer_double.txt" using 1:2 with lines title "Gi", "" using 1:3 with lines title "Gi'"
plot [-20.0:5.0][-1.5:+1.5] "plot/scorer_double.txt" using 1:4 with lines title "Hi", "" using 1:5 with lines title "Hi'"

plot [-20.0:5.0][-1.0:+1.0] "plot/airy_float.txt" using 1:2 with lines title "Ai(x)", "" using 1:4 with lines title "Bi(x)"
plot [-20.0:5.0][-1.5:+1.5] "plot/airy_float.txt" using 1:2 with lines title "Ai(x)", "" using 1:3 with lines title "Ai'(x)"
plot [-20.0:5.0][-1.5:+1.5] "plot/airy_float.txt" using 1:4 with lines title "Bi(x)", "" using 1:5 with lines title "Bi'(x)"

plot [-20.0:5.0][-1.0:+1.0] "plot/airy_double.txt" using 1:2 with lines title "Ai(x)", "" using 1:4 with lines title "Bi(x)"
plot [-20.0:5.0][-1.5:+1.5] "plot/airy_double.txt" using 1:2 with lines title "Ai(x)", "" using 1:3 with lines title "Ai'(x)"
plot [-20.0:5.0][-1.5:+1.5] "plot/airy_double.txt" using 1:4 with lines title "Bi(x)", "" using 1:5 with lines title "Bi'(x)"

plot [-20.0:5.0][-1.0:+1.0] "plot/airy_long_double.txt" using 1:2 with lines title "Ai(x)", "" using 1:4 with lines title "Bi(x)"
plot [-20.0:5.0][-1.5:+1.5] "plot/airy_long_double.txt" using 1:2 with lines title "Ai(x)", "" using 1:3 with lines title "Ai'(x)"
plot [-20.0:5.0][-1.5:+1.5] "plot/airy_long_double.txt" using 1:4 with lines title "Bi(x)", "" using 1:5 with lines title "Bi'(x)"

plot [0.0:30.0][-0.5:+6.0] "plot/struve_float.txt" using 1:2 with lines title "H_0", "" using 1:4 with lines title "H_1", "" using 1:6 with lines title "H_2", "" using 1:8 with lines title "H_3"
plot [0.0:5.0][-0.5:+4.5] "plot/struve_float.txt" using 1:3 with lines title "L_0", "" using 1:5 with lines title "L_1", "" using 1:7 with lines title "L_2", "" using 1:9 with lines title "L_3"

plot [0.0:30.0][-0.5:+6.0] "plot/struve_double.txt" using 1:2 with lines title "H_0", "" using 1:4 with lines title "H_1", "" using 1:6 with lines title "H_2", "" using 1:8 with lines title "H_3"
plot [0.0:5.0][-0.5:+4.5] "plot/struve_double.txt" using 1:3 with lines title "L_0", "" using 1:5 with lines title "L_1", "" using 1:7 with lines title "L_2", "" using 1:9 with lines title "L_3"

plot [0.0:30.0][-0.5:+6.0] "plot/struve_long_double.txt" using 1:2 with lines title "H_0", "" using 1:4 with lines title "H_1", "" using 1:6 with lines title "H_2", "" using 1:8 with lines title "H_3"
plot [0.0:5.0][-0.5:+4.5] "plot/struve_long_double.txt" using 1:3 with lines title "L_0", "" using 1:5 with lines title "L_1", "" using 1:7 with lines title "L_2", "" using 1:9 with lines title "L_3"

plot [0.0:40.0][-1.5:1.5] "plot/kelvin_float.txt" using 1:2 with lines title "ber(x)", "" using 1:3 with lines title "bei(x)", "" using 1:4 with lines title "ker(x)", "" using 1:5 with lines title "kei(x)"
plot [0.0:40.0][-1.5:1.5] "plot/kelvin_order_float.txt" using 1:2 with lines title "ber(x)", "" using 1:3 with lines title "bei(x)", "" using 1:4 with lines title "ker(x)", "" using 1:5 with lines title "kei(x)"

plot [0.0:40.0][-1.5:1.5] "plot/kelvin_double.txt" using 1:2 with lines title "ber(x)", "" using 1:3 with lines title "bei(x)", "" using 1:4 with lines title "ker(x)", "" using 1:5 with lines title "kei(x)"
plot [0.0:40.0][-1.5:1.5] "plot/kelvin_order_double.txt" using 1:2 with lines title "ber(x)", "" using 1:3 with lines title "bei(x)", "" using 1:4 with lines title "ker(x)", "" using 1:5 with lines title "kei(x)"

plot [0.0:40.0][-1.5:1.5] "plot/kelvin_long_double.txt" using 1:2 with lines title "ber(x)", "" using 1:3 with lines title "bei(x)", "" using 1:4 with lines title "ker(x)", "" using 1:5 with lines title "kei(x)"
plot [0.0:40.0][-1.5:1.5] "plot/kelvin_order_long_double.txt" using 1:2 with lines title "ber(x)", "" using 1:3 with lines title "bei(x)", "" using 1:4 with lines title "ker(x)", "" using 1:5 with lines title "kei(x)"

/usr/local/bin/gnuplot
set xzeroaxis
set yzeroaxis
set grid

/usr/local/bin/gnuplot
set hidden3d
splot [-20.0:5.0][-5.0:5.0][-1.5:+1.5] "plot/airy_complex_double.txt" title "Ai"


/usr/local/bin/gnuplot
set hidden3d


splot [-5:5][-5:5][0:12] "plot/airy_complex_double.txt" index 0 with lines title "|Ai(z)|^{1/6}"
splot [-20:5][-5:5][-20:20] "plot/airy_complex_double.txt" index 1 with lines title "Re[Ai(z)]"
splot [-20:5][-5:5][-20:20] "plot/airy_complex_double.txt" index 2 with lines title "Im[Ai(z)]"
splot [-20:5][-5:5][-3.2:+3.2] "plot/airy_complex_double.txt" index 3 with lines title "Arg[Ai(z)]"
splot [-5:5][-5:5][0:12] "plot/airy_complex_double.txt" index 4 with lines title "|Bi(z)|^{1/6}"
splot [-20:5][-5:5][-20:20] "plot/airy_complex_double.txt" index 5 with lines title "Re[Bi(z)]"
splot [-20:5][-5:5][-20:20] "plot/airy_complex_double.txt" index 6 with lines title "Im[Bi(z)]"
splot [-20:5][-5:5][-3.2:+3.2] "plot/airy_complex_double.txt" index 7 with lines title "Arg[Bi(z)]"
splot [-20:5][-5:5][-2:2] "plot/airy_complex_double.txt" index 8 with lines title "Wronski"

