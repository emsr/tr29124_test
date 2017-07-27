
gnuplot

set termoption enhanced

set polar
set grid polar 15
set angles degrees
set rrange [0:1]
set trange [*:*]
plot  "test_weierstrass_ellint.txt" index 1 using 1:2 title "Re[g_2(q)]"
plot  "test_weierstrass_ellint.txt" index 1 using 1:3 title "Im[g_2(q)]"
plot  "test_weierstrass_ellint.txt" index 1 using 1:4 title "|g_2(q)|"
plot  "test_weierstrass_ellint.txt" index 1 using 1:5 title "Re[g_3(q)]"
plot  "test_weierstrass_ellint.txt" index 1 using 1:6 title "Im[g_3(q)]"
plot  "test_weierstrass_ellint.txt" index 1 using 1:7 title "|g_3(q)|"


set hidden3d

set palette model RGB defined (-1 "black", -0.1 "blue", 0 "white", +0.1 "red", 1 "white")
set colorbox
set pm3d at bs corners2color geomean

splot [-6.25:6.25][-6.25:6.25][-50.0:50.0] "test_weierstrass_ellint.txt" index 0 using 1:2:3 with pm3d title "Re[P(q,z)]"
splot [-6.25:6.25][-6.25:6.25][-50.0:50.0] "test_weierstrass_ellint.txt" index 0 using 1:2:4 with pm3d title "Im[P(q,z)]"

set palette model RGB defined (0 "blue", 0.1 "red", 1 "white")
splot [-6.25:6.25][-6.25:6.25][0.0:50.0] "test_weierstrass_ellint.txt" index 0 using 1:2:5 with pm3d title "|P(q,z)|"
