gnuplot

set termoption enhanced

set title "e^{+2x^{3/2}/3}Ai(x)"
plot [0:10][0:0.5] \
  "test_airy_scaled.txt" index 0 using 1:2 with lines title "gnu", \
                      "" index 0 using 1:3 with lines title "mul", \
                      "" index 0 using 1:4 with lines title "me"

set title "e^{-2x^{3/2}/3}Bi(x)"
plot [0:10][0:1] \
  "test_airy_scaled.txt" index 0 using 1:5 with lines title "gnu", \
                      "" index 0 using 1:6 with lines title "mul", \
                      "" index 0 using 1:7 with lines title "me"
