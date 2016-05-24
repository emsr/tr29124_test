
/usr/local/bin/gnuplot --persist -e "plot \
 'plot/airy_float.txt' using 1:2 with lines title 'Ai', \
        	    '' using 1:3 with lines title 'Ai''', \
        	    '' using 1:4 with lines title 'Bi', \
        	    '' using 1:5 with lines title 'Bi'''"

/usr/local/bin/gnuplot --persist -e "plot [-20.0:5.0][-1.5:+1.5] \
 'plot/airy_double.txt' using 1:2 with lines title 'Ai', \
        	     '' using 1:3 with lines title 'Ai''', \
        	     '' using 1:4 with lines title 'Bi', \
        	     '' using 1:5 with lines title 'Bi'''"

/usr/local/bin/gnuplot --persist -e "plot [-20.0:1.0][-1.0:+1.0] \
 'plot/airy_long_double.txt' using 1:2 with lines title 'Ai', \
                	  '' using 1:4 with lines title 'Bi'"

/usr/local/bin/gnuplot --persist -e "plot [-20.0:5.0][-1.5:+1.5] \
 'plot/airy_long_double.txt' using 1:2 with lines title 'Ai', \
                	  '' using 1:3 with lines title 'Ai'''"

/usr/local/bin/gnuplot --persist -e "plot [-20.0:5.0][-1.5:+1.5] \
 'plot/airy_long_double.txt' using 1:4 with lines title 'Bi', \
                	  '' using 1:5 with lines title 'Bi'''"


/usr/local/bin/gnuplot --persist -e "plot [-20.0:1.0][-1.0:+1.0] \
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




plot [-20.0:5.0][-1.0:+1.0] 'plot/airy_float.txt' using 1:2 with lines title 'Ai(x)', '' using 1:4 with lines title 'Bi(x)'
plot [-20.0:5.0][-1.5:+1.5] 'plot/airy_float.txt' using 1:2 with lines title 'Ai(x)',  '' using 1:3 with lines title 'Ai''(x)'
plot [-20.0:5.0][-1.5:+1.5] 'plot/airy_float.txt' using 1:4 with lines title 'Bi(x)',  '' using 1:5 with lines title 'Bi''(x)'

plot [-20.0:5.0][-1.0:+1.0] 'plot/airy_double.txt' using 1:2 with lines title 'Ai(x)', '' using 1:4 with lines title 'Bi(x)'
plot [-20.0:5.0][-1.5:+1.5] 'plot/airy_double.txt' using 1:2 with lines title 'Ai(x)',  '' using 1:3 with lines title 'Ai''(x)'
plot [-20.0:5.0][-1.5:+1.5] 'plot/airy_double.txt' using 1:4 with lines title 'Bi(x)',  '' using 1:5 with lines title 'Bi''(x)'

plot [-20.0:5.0][-1.0:+1.0] 'plot/airy_long_double.txt' using 1:2 with lines title 'Ai(x)', '' using 1:4 with lines title 'Bi(x)'
plot [-20.0:5.0][-1.5:+1.5] 'plot/airy_long_double.txt' using 1:2 with lines title 'Ai(x)',  '' using 1:3 with lines title 'Ai''(x)'
plot [-20.0:5.0][-1.5:+1.5] 'plot/airy_long_double.txt' using 1:4 with lines title 'Bi(x)',  '' using 1:5 with lines title 'Bi''(x)'


/usr/local/bin/gnuplot
set hidden3d
splot [-20.0:5.0][-5.0:5.0][-1.5:+1.5] "plot/airy_complex_double.txt" title "Ai"


/usr/local/bin/gnuplot
set hidden3d


splot [-20:5][-5:5][0:2] "plot/airy_complex_double.txt" index 0 with lines title "|Ai(z)|"
splot [-20:5][-5:5][-2:2] "plot/airy_complex_double.txt" index 1 with lines title "Re[Ai(z)]"
splot [-20:5][-5:5][-2:2] "plot/airy_complex_double.txt" index 2 with lines title "Im[Ai(z)]"
splot [-20:5][-5:5][-3.2:+3.2] "plot/airy_complex_double.txt" index 3 with lines title "Arg[Ai(z)]"
splot [-20:5][-5:5][0:2] "plot/airy_complex_double.txt" index 4 with lines title "|Bi(z)|"
splot [-20:5][-5:5][-2:2] "plot/airy_complex_double.txt" index 5 with lines title "Re[Bi(z)]"
splot [-20:5][-5:5][-2:2] "plot/airy_complex_double.txt" index 6 with lines title "Im[Bi(z)]"
splot [-20:5][-5:5][-3.2:+3.2] "plot/airy_complex_double.txt" index 7 with lines title "Arg[Bi(z)]"
splot [-20:5][-5:5][-2:2] "plot/airy_complex_double.txt" index 8 with lines title "Wronski"

