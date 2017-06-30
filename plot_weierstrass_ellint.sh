
gnuplot

set termoption enhanced

set hidden3d

set palette model RGB defined (-1 "black", -0.1 "blue", 0 "white", +0.1 "red", 1 "white")
set colorbox
set pm3d at bs corners2color geomean

splot [-6.25:6.25][-6.25:6.25][-0.06:0.06] "test_weierstrass_ellint.txt" index 0 with pm3d title "P(q,z)"
