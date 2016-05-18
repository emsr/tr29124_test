
gnuplot --persist -e "plot \
 'plot/airy_float.txt' using 1:2 with lines title 'Ai', \
        	    '' using 1:3 with lines title 'Ai''', \
        	    '' using 1:4 with lines title 'Bi', \
        	    '' using 1:5 with lines title 'Bi'''"

gnuplot --persist -e "plot \
 'plot/airy_double.txt' using 1:2 with lines title 'Ai', \
        	     '' using 1:3 with lines title 'Ai''', \
        	     '' using 1:4 with lines title 'Bi', \
        	     '' using 1:5 with lines title 'Bi'''"

gnuplot --persist -e "plot [-20.0:1.0][-1.0:+1.0] \
 'plot/airy_long_double.txt' using 1:2 with lines title 'Ai', \
                	  '' using 1:4 with lines title 'Bi'"

gnuplot --persist -e "plot [-20.0:5.0][-1.5:+1.5] \
 'plot/airy_long_double.txt' [0:1200] using 1:2 with lines title 'Ai', \
                	  '' [0:1200] using 1:3 with lines title 'Ai'''"

gnuplot --persist -e "plot [-20.0:5.0][-1.5:+1.5] \
 'plot/airy_long_double.txt' [0:1200] using 1:4 with lines title 'Bi', \
                	  '' [0:1200] using 1:5 with lines title 'Bi'''"
