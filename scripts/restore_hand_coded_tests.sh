git checkout -- \
  testsuite/tr1/5_numerical_facilities/special_functions/02_assoc_legendre/pr83140.cc \
  testsuite/special_functions/02_assoc_legendre/pr83140.cc \
  testsuite/tr1/5_numerical_facilities/special_functions/15_expint/pr68397.cc \
  testsuite/special_functions/14_expint/pr68397.cc \
  testsuite/special_functions/07_cyl_bessel_i/pr56216.cc \
  testsuite/special_functions/08_cyl_bessel_j/check_origin.cc \
  testsuite/special_functions/10_cyl_neumann/check_origin.cc \
  testsuite/ext/special_functions/theta_1/check_value.cc \
  testsuite/ext/special_functions/theta_2/check_value.cc \
  testsuite/ext/special_functions/theta_3/check_value.cc \
  testsuite/ext/special_functions/theta_4/check_value.cc \
  testsuite/ext/special_functions/sinhint/check_shi.cc \
  testsuite/ext/special_functions/ellnome/check_value.cc  \
  testsuite/ext/special_functions/dirichlet_lambda/check_value.cc  \
  testsuite/ext/special_functions/dirichlet_beta/check_value.cc \
  testsuite/ext/special_functions/cyl_hankel_1/pr56216.cc \
  testsuite/ext/special_functions/cyl_hankel_2/pr56216.cc \
  testsuite/ext/special_functions/coshint/check_chi.cc \
  testsuite/ext/special_functions/gamma_reciprocal/check_value.cc \
  testsuite/ext/special_functions/eulerian_2/check_value.cc

cp  testsuite/tr1/5_numerical_facilities/special_functions/02_assoc_legendre/pr83140.cc   check/tr1_pr83140.cc
cp  testsuite/special_functions/02_assoc_legendre/pr83140.cc   check/pr83140.cc
cp  testsuite/tr1/5_numerical_facilities/special_functions/15_expint/pr68397.cc   check/tr1_pr68397.cc
cp  testsuite/special_functions/14_expint/pr68397.cc                  check/pr68397.cc
cp  testsuite/special_functions/07_cyl_bessel_i/pr56216.cc            check/pr56216_cyl_bessel_i.cc
cp  testsuite/special_functions/08_cyl_bessel_j/check_origin.cc       check/origin_cyl_bessel_j.cc
cp  testsuite/special_functions/10_cyl_neumann/check_origin.cc        check/origin_cyl_neumann.cc
cp  testsuite/ext/special_functions/theta_1/check_value.cc            check/check_theta_1.cc
cp  testsuite/ext/special_functions/theta_2/check_value.cc            check/check_theta_2.cc
cp  testsuite/ext/special_functions/theta_3/check_value.cc            check/check_theta_3.cc
cp  testsuite/ext/special_functions/theta_4/check_value.cc            check/check_theta_4.cc
cp  testsuite/ext/special_functions/sinhint/check_shi.cc              check/check_shi.cc
cp  testsuite/ext/special_functions/ellnome/check_value.cc            check/check_ellnome.cc
cp  testsuite/ext/special_functions/dirichlet_lambda/check_value.cc   check/check_dirichlet_lambda.cc
cp  testsuite/ext/special_functions/dirichlet_beta/check_value.cc     check/check_dirichlet_beta.cc
cp  testsuite/ext/special_functions/cyl_hankel_1/pr56216.cc           check/pr56216_cyl_hankel_1.cc
cp  testsuite/ext/special_functions/cyl_hankel_2/pr56216.cc           check/pr56216_cyl_hankel_2.cc
cp  testsuite/ext/special_functions/coshint/check_chi.cc              check/check_chi.cc
cp  testsuite/ext/special_functions/gamma_reciprocal/check_value.cc   check/check_gamma_reciprocal.cc
cp  testsuite/ext/special_functions/eulerian_2/check_value.cc         check/check_eulerian_2.cc

git add \
  check/tr1_pr83140.cc \
  check/pr83140.cc \
  check/tr1_pr68397.cc \
  check/pr68397.cc \
  check/pr56216_cyl_bessel_i.cc \
  check/origin_cyl_bessel_j.cc \
  check/origin_cyl_neumann.cc \
  check/check_theta_1.cc \
  check/check_theta_2.cc \
  check/check_theta_3.cc \
  check/check_theta_4.cc \
  check/check_shi.cc \
  check/check_ellnome.cc \
  check/check_dirichlet_lambda.cc \
  check/check_dirichlet_beta.cc \
  check/pr56216_cyl_hankel_1.cc \
  check/pr56216_cyl_hankel_2.cc \
  check/check_chi.cc \
  check/check_gamma_reciprocal.cc \
  check/check_eulerian_2.cc