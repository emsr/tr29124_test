
# If the handcrafted tests get clobbered...

git checkout -- \
 check/check_airy_ai.cc \
 check/check_airy_bi.cc \
 check/check_beta.cc \
 check/check_chi.cc \
 check/check_dirichlet_beta.cc \
 check/check_dirichlet_lambda.cc \
 check/check_ellnome.cc  \
 check/check_shi.cc \
 check/check_theta_1.cc \
 check/check_theta_2.cc \
 check/check_theta_3.cc \
 check/check_theta_4.cc \
 check/check_gamma_reciprocal.cc

git checkout -- \
 check/complex_airy_ai.cc \
 check/complex_airy_bi.cc \
 check/complex_ellint_rc.cc \
 check/complex_ellint_rd.cc \
 check/complex_ellint_rf.cc \
 check/complex_ellint_rg.cc \
 check/complex_ellint_rj.cc

git checkout -- \
 check/origin_cyl_bessel_j.cc \
 check/origin_cyl_neumann.cc

git checkout -- \
 check/pr56216_cyl_bessel_i.cc \
 check/pr56216_cyl_hankel_1.cc \
 check/pr56216_cyl_hankel_2.cc
