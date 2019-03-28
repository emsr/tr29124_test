gnuplot
set termoption enhanced
set xzeroaxis
set yzeroaxis
set grid

set palette model RGB defined (0.0 "blue", 0.3 "green", 0.7 "yellow", 1 "red")
set colorbox
set pm3d at bs corners2color geomean



# Kelvin
set title "Kelvin functions ber(x), bei(x), ker(x), kei(x)"
set xlabel "x"
plot [0.0:40.0][-1.5:1.5] \
    "../plot_data/kelvin_float.txt" using 1:2 with lines title "ber(x)", \
                         "" using 1:3 with lines title "bei(x)", \
                         "" using 1:4 with lines title "ker(x)", \
                         "" using 1:5 with lines title "kei(x)"

set title "Kelvin functions ber(x), bei(x), ker(x), kei(x)"
set xlabel "x"
plot [0.0:40.0][-1.5:1.5] \
    "../plot_data/kelvin_double.txt" using 1:2 with lines title "ber(x)", \
                          "" using 1:3 with lines title "bei(x)", \
                          "" using 1:4 with lines title "ker(x)", \
                          "" using 1:5 with lines title "kei(x)"

set title "Kelvin functions ber(x), bei(x), ker(x), kei(x)"
set xlabel "x"
plot [0.0:40.0][-1.5:1.5] \
    "../plot_data/kelvin_long_double.txt" using 1:2 with lines title "ber(x)", \
                               "" using 1:3 with lines title "bei(x)", \
                               "" using 1:4 with lines title "ker(x)", \
                               "" using 1:5 with lines title "kei(x)"


set title "Kelvin functions ber_{/Symbol n}(x), bei_{/Symbol n}(x), ker_{/Symbol n}(x), kei_{/Symbol n}(x)"
set xlabel "x"
plot [0.0:40.0][-1.5:1.5] \
    "../plot_data/kelvin_order_float.txt" using 1:2 with lines title "ber_{/Symbol n}(x)", \
                               "" using 1:3 with lines title "bei_{/Symbol n}(x)", \
                               "" using 1:4 with lines title "ker_{/Symbol n}(x)", \
                               "" using 1:5 with lines title "kei_{/Symbol n}(x)"

set title "Kelvin functions ber_{/Symbol n}(x), bei_{/Symbol n}(x), ker_{/Symbol n}(x), kei_{/Symbol n}(x)"
set xlabel "x"
plot [0.0:40.0][-1.5:1.5] \
    "../plot_data/kelvin_order_double.txt" using 1:2 with lines title "ber_{/Symbol n}(x)", \
                        	"" using 1:3 with lines title "bei_{/Symbol n}(x)", \
                        	"" using 1:4 with lines title "ker_{/Symbol n}(x)", \
                        	"" using 1:5 with lines title "kei_{/Symbol n}(x)"

set title "Kelvin functions ber_{/Symbol n}(x), bei_{/Symbol n}(x), ker_{/Symbol n}(x), kei_{/Symbol n}(x)"
set xlabel "x"
plot [0.0:40.0][-1.5:1.5] \
    "../plot_data/kelvin_order_long_double.txt" using 1:2 with lines title "ber_{/Symbol n}(x)", \
                                     "" using 1:3 with lines title "bei_{/Symbol n}(x)", \
                                     "" using 1:4 with lines title "ker_{/Symbol n}(x)", \
                                     "" using 1:5 with lines title "kei_{/Symbol n}(x)"
