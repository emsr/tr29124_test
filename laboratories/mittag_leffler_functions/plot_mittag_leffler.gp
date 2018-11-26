
gnuplot

set termoption enhanced

set xzeroaxis
set yzeroaxis
set grid

asymp(x, c) = c

# Figure 1
set title "Mittag-Leffler function E_{1/4,1}(x)"
set xlabel "-x"
plot [0.0:10.0][-1.0:1.2] \
                    "test_mittag_leffler.txt" index 0 using 1:2 with lines title "E_{{/Symbol a},{/Symbol b}}(x)", \
                    "" index 0 using 1:3 with lines title "E'_{{/Symbol a},{/Symbol b}}(x)"

# Figure 2
set title "Mittag-Leffler function E_{7/4,1}(x)"
set xlabel "-x"
plot [0.0:50.0][-0.7:1.0] \
                    "test_mittag_leffler.txt" index 1 using 1:2 with lines title "E_{{/Symbol a},{/Symbol b}}(x)", \
                    "" index 1 using 1:3 with lines title "E'_{{/Symbol a},{/Symbol b}}(x)"

# Figure 3
set title "Mittag-Leffler function E_{9/4,1}(x)"
set xlabel "-x"
plot [0.0:100.0][-2.0:3.0] \
                    "test_mittag_leffler.txt" index 2 using 1:2 with lines title "E_{{/Symbol a},{/Symbol b}}(x)", \
                    "" index 2 using 1:3 with lines title "E'_{{/Symbol a},{/Symbol b}}(x)"

# Figure 4
set title "Mittag-Leffler function |E_{3/4,1}(z)|"
set xlabel "|z|, arg(z) = {/Symbol a}{/Symbol p}/4"
plot [0.0:5.0][0.0:400.0] \
                    "test_mittag_leffler.txt" index 3 using 1:2 with lines title "|E_{{/Symbol a},{/Symbol b}}(z)|"

# Figure 5
set title "Mittag-Leffler function |E_{3/4,1}(z)|"
set xlabel "|z|, arg(z) = {/Symbol a}{/Symbol p}/2"
plot [0.0:50.0][1.29:1.38] \
                    "test_mittag_leffler.txt" index 4 using 1:2 with lines title "|E_{{/Symbol a},{/Symbol b}}(z)|", \
                    asymp(x, 4./3.) with lines title "1/{/Symbol a} = 4/3"

# Figure 6
set title "Mittag-Leffler function |E_{3/4,1}(z)|"
set xlabel "|z|, arg(z) = 3{/Symbol a}{/Symbol p}/4"
plot [0.0:50.0][0.000:0.175] \
                    "test_mittag_leffler.txt" index 5 using 1:2 with lines title "|E_{{/Symbol a},{/Symbol b}}(z)|"

# Figure 7
set title "Mittag-Leffler function |E_{3/4,1}(z)|"
set xlabel "|z|, arg(z) = {/Symbol p}"
plot [0.0:20.0][0.00:1.00] \
                    "test_mittag_leffler.txt" index 6 using 1:2 with lines title "|E_{{/Symbol a},{/Symbol b}}(z)|"

# Figure 8
set title "Mittag-Leffler function |E_{5/4,1}(z)|"
set xlabel "|z|, arg(z) = {/Symbol a}{/Symbol p}/4"
plot [0.0:10.0][0.00:70.00] \
                    "test_mittag_leffler.txt" index 7 using 1:2 with lines title "|E_{{/Symbol a},{/Symbol b}}(x)|"

# Figure 9
set title "Mittag-Leffler function |E_{5/4,1}(z)|"
set xlabel "|z|, arg(z) = {/Symbol a}{/Symbol p}/2"
plot [0.0:50.0][0.72:0.85] \
                    "test_mittag_leffler.txt" index 8 using 1:2 with lines title "|E_{{/Symbol a},{/Symbol b}}(x)|", \
                    asymp(x, 4./5.) with lines title "1/{/Symbol a} = 4/5"

# Figure 10
set title "Mittag-Leffler function |E_{5/4,1}(z)|"
set xlabel "|z|, arg(z) = 3{/Symbol a}{/Symbol p}/2"
plot [0.0:50.0][0.00:0.27] \
                    "test_mittag_leffler.txt" index 9 using 1:2 with lines title "|E_{{/Symbol a},{/Symbol b}}(x)|"

# Figure 11
set title "Mittag-Leffler function |E_{5/4,1}(z)|"
set xlabel "|z|, arg(z) = {/Symbol p}"
plot [0.0:100.0][-0.11:0.15] \
                    "test_mittag_leffler.txt" index 10 using 1:3 with lines title "Re{E_{{/Symbol a},{/Symbol b}}(z)}"
