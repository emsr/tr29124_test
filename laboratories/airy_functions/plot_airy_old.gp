
# Output from airy_toy -DOLD.

gnuplot --persist -e \
    "plot [-20.0:5.0][-1.5:+1.5] \
 'plot_data/airy_float.txt' using 1:2 with lines title 'Ai', \
        	    '' using 1:3 with lines title 'Ai''', \
        	    '' using 1:4 with lines title 'Bi', \
        	    '' using 1:5 with lines title 'Bi'''"

gnuplot --persist -e \
    "plot [-20.0:5.0][-1.5:+1.5] \
 'plot_data/airy_double.txt' using 1:2 with lines title 'Ai', \
        	     '' using 1:3 with lines title 'Ai''', \
        	     '' using 1:4 with lines title 'Bi', \
        	     '' using 1:5 with lines title 'Bi'''"

gnuplot --persist -e \
    "plot [-20.0:5.0][-1.0:+1.0] \
 'plot_data/airy_long_double.txt' using 1:2 with lines title 'Ai', \
                	  '' using 1:4 with lines title 'Bi'"

gnuplot --persist -e \
    "plot [-20.0:5.0][-1.5:+1.5] \
 'plot_data/airy_long_double.txt' using 1:2 with lines title 'Ai', \
                	  '' using 1:3 with lines title 'Ai'''"

gnuplot --persist -e \
    "plot [-20.0:5.0][-1.5:+1.5] \
 'plot_data/airy_long_double.txt' using 1:4 with lines title 'Bi', \
                	  '' using 1:5 with lines title 'Bi'''"


gnuplot --persist -e \
    "plot [-20.0:5.0][-1.0:+1.0] \
 'plot_data/scorer_double.txt' using 1:2 with lines title 'Gi', \
                	  '' using 1:4 with lines title 'Hi'"

gnuplot --persist -e \
    "plot [-20.0:5.0][-1.5:+1.5] \
 'plot_data/scorer_double.txt' using 1:2 with lines title 'Gi', \
                	  '' using 1:3 with lines title 'Gi'''"

gnuplot --persist -e \
    "plot [-20.0:5.0][-1.5:+1.5] \
 'plot_data/scorer_double.txt' using 1:4 with lines title 'Hi', \
                	  '' using 1:5 with lines title 'Hi'''"


gnuplot --persist -e \
    "splot [-20.0:5.0][-5.0:5.0][-1.5:+1.5] \
 'plot_data/airy_complex_double.txt' title 'Ai'"


gnuplot
set termoption enhanced
set xzeroaxis
set yzeroaxis
set grid

set palette model RGB defined (0.0 "blue", 0.3 "green", 0.7 "yellow", 1 "red")
set colorbox
set pm3d at bs corners2color geomean


# FGH
plot [-5.0:2.0][-1.5:+1.5] \
    "../plot_data/fgh_double_old.txt" using 1:2 with lines title "fai", \
                       "" using 1:3 with lines title "gai", \
                       "" using 1:4 with lines title "hai"

plot [-5.0:2.0][-2.5:+2.5] \
    "../plot_data/fgh_double_old.txt" using 1:5 with lines title "fai'", \
                       "" using 1:6 with lines title "gai'", \
                       "" using 1:7 with lines title "hai'"

plot [-5.0:2.0][-2.5:+2.5] \
    "../plot_data/fgh_double_old.txt" using 1:($2+$3+$4) with lines title "Hi"

plot [-5.0:2.0][-2.5:+2.5] \
    "../plot_data/fgh_double_old.txt" using 1:($5+$6+$7) with lines title "Hi'"


# Scorer
plot [-20.0:5.0][-1.0:+1.0] \
    "../plot_data/scorer_double_old.txt" using 1:2 with lines title "Gi", \
                          "" using 1:4 with lines title "Hi"

plot [-20.0:5.0][-1.5:+1.5] \
    "../plot_data/scorer_double_old.txt" using 1:2 with lines title "Gi", \
                          "" using 1:3 with lines title "Gi'"

plot [-20.0:5.0][-1.5:+1.5] \
    "../plot_data/scorer_double_old.txt" using 1:4 with lines title "Hi", \
                          "" using 1:5 with lines title "Hi'"


# Scorer
plot [-20.0:5.0][-1.0:+1.0] \
    "../plot_data/scorer_float_old.txt" using 1:2 with lines title "Gi", \
                         "" using 1:3 with lines title "Gi'", \
                         "" using 1:4 with lines title "Hi", \
                         "" using 1:5 with lines title "Hi'"

plot [-20.0:5.0][-1.0:+1.0] \
    "../plot_data/scorer_double_old.txt" using 1:2 with lines title "Gi", \
                          "" using 1:3 with lines title "Gi'", \
                          "" using 1:4 with lines title "Hi", \
                          "" using 1:5 with lines title "Hi'"

plot [-20.0:5.0][-1.0:+1.0] \
    "../plot_data/scorer_long_double_old.txt" using 1:2 with lines title "Gi", \
                               "" using 1:3 with lines title "Gi'", \
                               "" using 1:4 with lines title "Hi", \
                               "" using 1:5 with lines title "Hi'"


# Airy
plot [-20.0:5.0][-1.0:+1.0] \
    "../plot_data/airy_float_old.txt" using 1:2 with lines title "Ai(x)", \
                       "" using 1:4 with lines title "Bi(x)"

plot [-20.0:5.0][-1.5:+1.5] \
    "../plot_data/airy_float_old.txt" using 1:2 with lines title "Ai(x)", \
                       "" using 1:3 with lines title "Ai'(x)"

plot [-20.0:5.0][-1.5:+1.5] \
    "../plot_data/airy_float_old.txt" using 1:4 with lines title "Bi(x)", \
                       "" using 1:5 with lines title "Bi'(x)"


plot [-20.0:5.0][-1.0:+1.0] \
    "../plot_data/airy_double_old.txt" using 1:2 with lines title "Ai(x)", \
                        "" using 1:4 with lines title "Bi(x)"

plot [-20.0:5.0][-1.5:+1.5] \
    "../plot_data/airy_double_old.txt" using 1:2 with lines title "Ai(x)", \
                        "" using 1:3 with lines title "Ai'(x)"

plot [-20.0:5.0][-1.5:+1.5] \
    "../plot_data/airy_double_old.txt" using 1:4 with lines title "Bi(x)", \
                        "" using 1:5 with lines title "Bi'(x)"


plot [-20.0:5.0][-1.0:+1.0] \
    "../plot_data/airy_long_double_old.txt" using 1:2 with lines title "Ai(x)", \
                             "" using 1:4 with lines title "Bi(x)"

plot [-20.0:5.0][-1.5:+1.5] \
    "../plot_data/airy_long_double_old.txt" using 1:2 with lines title "Ai(x)", \
                             "" using 1:3 with lines title "Ai'(x)"

plot [-20.0:5.0][-1.5:+1.5] \
    "../plot_data/airy_long_double_old.txt" using 1:4 with lines title "Bi(x)", \
                             "" using 1:5 with lines title "Bi'(x)"


# Struve
plot [0.0:30.0][-0.5:+6.0] \
    "../plot_data/struve_float_old.txt" using 1:2 with lines title "{/:Bold H}_0", \
                         "" using 1:4 with lines title "{/:Bold H}_1", \
                         "" using 1:6 with lines title "{/:Bold H}_2", \
                         "" using 1:8 with lines title "{/:Bold H}_3"

plot [0.0:5.0][-0.5:+4.5] \
    "../plot_data/struve_float_old.txt" using 1:3 with lines title "{/:Bold L}_0", \
                         "" using 1:5 with lines title "{/:Bold L}_1", \
                         "" using 1:7 with lines title "{/:Bold L}_2", \
                         "" using 1:9 with lines title "{/:Bold L}_3"

plot [0.0:30.0][-0.5:+6.0] \
    "../plot_data/struve_double_old.txt" using 1:2 with lines title "{/:Bold H}_0", \
                          "" using 1:4 with lines title "{/:Bold H}_1", \
                          "" using 1:6 with lines title "{/:Bold H}_2", \
                          "" using 1:8 with lines title "{/:Bold H}_3"

plot [0.0:5.0][-0.5:+4.5] \
    "../plot_data/struve_double_old.txt" using 1:3 with lines title "{/:Bold L}_0", \
                          "" using 1:5 with lines title "{/:Bold L}_1", \
                          "" using 1:7 with lines title "{/:Bold L}_2", \
                          "" using 1:9 with lines title "{/:Bold L}_3"

plot [0.0:30.0][-0.5:+6.0] \
    "../plot_data/struve_long_double_old.txt" using 1:2 with lines title "{/:Bold H}_0", \
                               "" using 1:4 with lines title "{/:Bold H}_1", \
                               "" using 1:6 with lines title "{/:Bold H}_2", \
                               "" using 1:8 with lines title "{/:Bold H}_3"

plot [0.0:5.0][-0.5:+4.5] \
    "../plot_data/struve_long_double_old.txt" using 1:3 with lines title "{/:Bold L}_0", \
                               "" using 1:5 with lines title "{/:Bold L}_1", \
                               "" using 1:7 with lines title "{/:Bold L}_2", \
                               "" using 1:9 with lines title "{/:Bold L}_3"


gnuplot
set xzeroaxis
set yzeroaxis
set grid

set palette model RGB defined (0.0 "blue", 0.3 "green", 0.7 "yellow", 1 "red")
set colorbox
set pm3d at sb corners2color geomean

set hidden3d


splot [-20.0:5.0][-5.0:5.0][-1.5:+1.5] "../plot_data/airy_complex_double_old.txt" title "Ai"


# Airy complex float
splot [-5:5][-5:5][0:12] "../plot_data/airy_complex_float_old.txt" index 0 with pm3d title "|Ai(z)|^{1/6}"
splot [-20:5][-5:5][-20:20] "../plot_data/airy_complex_float_old.txt" index 1 with pm3d title "Re[Ai(z)]"
splot [-20:5][-5:5][-20:20] "../plot_data/airy_complex_float_old.txt" index 2 with pm3d title "Im[Ai(z)]"
splot [-20:5][-5:5][-3.2:+3.2] "../plot_data/airy_complex_float_old.txt" index 3 with pm3d title "Arg[Ai(z)]"
splot [-5:5][-5:5][0:12] "../plot_data/airy_complex_float_old.txt" index 4 with pm3d title "|Bi(z)|^{1/6}"
splot [-20:5][-5:5][-20:20] "../plot_data/airy_complex_float_old.txt" index 5 with pm3d title "Re[Bi(z)]"
splot [-20:5][-5:5][-20:20] "../plot_data/airy_complex_float_old.txt" index 6 with pm3d title "Im[Bi(z)]"
splot [-20:5][-5:5][-3.2:+3.2] "../plot_data/airy_complex_float_old.txt" index 7 with pm3d title "Arg[Bi(z)]"
splot [-20:5][-5:5][-2:2] "../plot_data/airy_complex_float_old.txt" index 8 with pm3d title "Wronski"


# Airy complex double
splot [-5:5][-5:5][0:12] "../plot_data/airy_complex_double_old.txt" index 0 with pm3d title "|Ai(z)|^{1/6}"
splot [-20:5][-5:5][-20:20] "../plot_data/airy_complex_double_old.txt" index 1 with pm3d title "Re[Ai(z)]"
splot [-20:5][-5:5][-20:20] "../plot_data/airy_complex_double_old.txt" index 2 with pm3d title "Im[Ai(z)]"
splot [-20:5][-5:5][-3.2:+3.2] "../plot_data/airy_complex_double_old.txt" index 3 with pm3d title "Arg[Ai(z)]"
splot [-5:5][-5:5][0:12] "../plot_data/airy_complex_double_old.txt" index 4 with pm3d title "|Bi(z)|^{1/6}"
splot [-20:5][-5:5][-20:20] "../plot_data/airy_complex_double_old.txt" index 5 with pm3d title "Re[Bi(z)]"
splot [-20:5][-5:5][-20:20] "../plot_data/airy_complex_double_old.txt" index 6 with pm3d title "Im[Bi(z)]"
splot [-20:5][-5:5][-3.2:+3.2] "../plot_data/airy_complex_double_old.txt" index 7 with pm3d title "Arg[Bi(z)]"
splot [-20:5][-5:5][-2:2] "../plot_data/airy_complex_double_old.txt" index 8 with pm3d title "Wronski"


# Airy complex long double
splot [-5:5][-5:5][0:12] "../plot_data/airy_complex_long_double_old.txt" index 0 with pm3d title "|Ai(z)|^{1/6}"
splot [-20:5][-5:5][-20:20] "../plot_data/airy_complex_long_double_old.txt" index 1 with pm3d title "Re[Ai(z)]"
splot [-20:5][-5:5][-20:20] "../plot_data/airy_complex_long_double_old.txt" index 2 with pm3d title "Im[Ai(z)]"
splot [-20:5][-5:5][-3.2:+3.2] "../plot_data/airy_complex_long_double_old.txt" index 3 with pm3d title "Arg[Ai(z)]"
splot [-5:5][-5:5][0:12] "../plot_data/airy_complex_long_double_old.txt" index 4 with pm3d title "|Bi(z)|^{1/6}"
splot [-20:5][-5:5][-20:20] "../plot_data/airy_complex_long_double_old.txt" index 5 with pm3d title "Re[Bi(z)]"
splot [-20:5][-5:5][-20:20] "../plot_data/airy_complex_long_double_old.txt" index 6 with pm3d title "Im[Bi(z)]"
splot [-20:5][-5:5][-3.2:+3.2] "../plot_data/airy_complex_long_double_old.txt" index 7 with pm3d title "Arg[Bi(z)]"
splot [-20:5][-5:5][-2:2] "../plot_data/airy_complex_long_double_old.txt" index 8 with pm3d title "Wronski"

