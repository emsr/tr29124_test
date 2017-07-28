
gnuplot

set termoption enhanced

#set polar
#set grid polar 15
#set angles radians
#set rrange [0.0:1.0]
#set trange [0.0:2.0*pi]
splot [-1.0:1.0][-1.0:1.0][-0.5:0.5] "test_weierstrass_ellint.txt" index 1 using 1:2:3 with pm3d title "Re[g_2(q)]"
splot [-1.0:1.0][-1.0:1.0][-0.5:0.5] "test_weierstrass_ellint.txt" index 1 using 1:2:4 with pm3d title "Im[g_2(q)]"
splot [-1.0:1.0][-1.0:1.0][0.0:0.5] "test_weierstrass_ellint.txt" index 1 using 1:2:5 with pm3d title "|g_2(q)|"
splot [-1.0:1.0][-1.0:1.0][-0.5:0.5] "test_weierstrass_ellint.txt" index 1 using 1:2:6 with pm3d title "Re[g_3(q)]"
splot [-1.0:1.0][-1.0:1.0][-0.5:0.5] "test_weierstrass_ellint.txt" index 1 using 1:2:7 with pm3d title "Im[g_3(q)]"
splot [-1.0:1.0][-1.0:1.0][0.0:0.5] "test_weierstrass_ellint.txt" index 1 using 1:2:8 with pm3d title "|g_3(q)|"


set hidden3d

set palette model RGB defined (-1 "black", -0.1 "blue", 0 "white", +0.1 "red", 1 "white")
set colorbox
set pm3d at bs corners2color geomean

splot [-6.25:6.25][-6.25:6.25][-50.0:50.0] "test_weierstrass_ellint.txt" index 0 using 1:2:3 with pm3d title "Re[P(q,z)]"
splot [-6.25:6.25][-6.25:6.25][-50.0:50.0] "test_weierstrass_ellint.txt" index 0 using 1:2:4 with pm3d title "Im[P(q,z)]"

set palette model RGB defined (0 "blue", 0.1 "red", 1 "white")
splot [-6.25:6.25][-6.25:6.25][0.0:50.0] "test_weierstrass_ellint.txt" index 0 using 1:2:5 with pm3d title "|P(q,z)|"
