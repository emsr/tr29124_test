
gnuplot

set termoption enhanced

set hidden3d

set palette model RGB defined (-1 "black", -0.1 "blue", 0 "white", +0.1 "red", 1 "white")
set colorbox
set pm3d at bs corners2color geomean

splot [-6.25:6.25][-6.25:6.25][-0.06:0.06] "test_weierstrass_ellint.txt" index 0 using 1:2:3 with pm3d title "Re[P(q,z)]"
splot [-6.25:6.25][-6.25:6.25][-0.06:0.06] "test_weierstrass_ellint.txt" index 0 using 1:2:4 with pm3d title "Im[P(q,z)]"

set palette model RGB defined (0 "blue", 0.1 "red", 1 "white")
splot [-6.25:6.25][-6.25:6.25][-0.06:0.06] "test_weierstrass_ellint.txt" index 0 using 1:2:5 with pm3d title "|P(q,z)|"
