git checkout -- \
  testsuite/special_functions/10_cyl_neumann/check_origin.cc \
  testsuite/special_functions/07_cyl_bessel_i/pr56216.cc \
  testsuite/special_functions/08_cyl_bessel_j/check_origin.cc \
  testsuite/ext/special_functions/theta_4/check_value.cc \
  testsuite/ext/special_functions/theta_3/check_value.cc \
  testsuite/ext/special_functions/theta_2/check_value.cc \
  testsuite/ext/special_functions/theta_1/check_value.cc \
  testsuite/ext/special_functions/sinhint/check_shi.cc \
  testsuite/ext/special_functions/ellnome/check_value.cc  \
  testsuite/ext/special_functions/dirichlet_lambda/check_value.cc  \
  testsuite/ext/special_functions/dirichlet_beta/check_value.cc \
  testsuite/ext/special_functions/cyl_hankel_2/pr56216.cc \
  testsuite/ext/special_functions/coshint/check_chi.cc \
  testsuite/ext/special_functions/cyl_hankel_1/pr56216.cc

cp  testsuite/special_functions/10_cyl_neumann/check_origin.cc        check/origin_cyl_neumann.cc
cp  testsuite/special_functions/07_cyl_bessel_i/pr56216.cc            check/pr56216_cyl_bessel_i.cc
cp  testsuite/special_functions/08_cyl_bessel_j/check_origin.cc       check/origin_cyl_bessel_j.cc
cp  testsuite/ext/special_functions/theta_4/check_value.cc            check/check_theta_4.cc
cp  testsuite/ext/special_functions/theta_3/check_value.cc            check/check_theta_3.cc
cp  testsuite/ext/special_functions/theta_2/check_value.cc            check/check_theta_2.cc
cp  testsuite/ext/special_functions/theta_1/check_value.cc            check/check_theta_1.cc
cp  testsuite/ext/special_functions/sinhint/check_shi.cc              check/check_shi.cc
cp  testsuite/ext/special_functions/ellnome/check_value.cc            check/check_ellnome.cc
cp  testsuite/ext/special_functions/dirichlet_lambda/check_value.cc   check/check_dirichlet_lambda.cc
cp  testsuite/ext/special_functions/dirichlet_beta/check_value.cc     check/check_dirichlet_beta.cc
cp  testsuite/ext/special_functions/cyl_hankel_2/pr56216.cc           check/pr56216_cyl_hankel_2.cc
cp  testsuite/ext/special_functions/coshint/check_chi.cc              check/check_chi.cc
cp  testsuite/ext/special_functions/cyl_hankel_1/pr56216.cc           check/pr56216_cyl_hankel_1.cc

git add \
  check/origin_cyl_neumann.cc \
  check/pr56216_cyl_bessel_i.cc \
  check/origin_cyl_bessel_j.cc \
  check/check_theta_4.cc \
  check/check_theta_3.cc \
  check/check_theta_2.cc \
  check/check_theta_1.cc \
  check/check_shi.cc \
  check/check_ellnome.cc \
  check/check_dirichlet_lambda.cc \
  check/check_dirichlet_beta.cc \
  check/pr56216_cyl_hankel_2.cc \
  check/check_chi.cc \
  check/pr56216_cyl_hankel_1.cc
