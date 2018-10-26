#!  /bin/bash

copy="cp -f"

suffix="_tr29124"
if [ $# -ge 1 ]; then
  suffix="$1"
fi

base_dir="$HOME/gcc${suffix}"
test_dir="${base_dir}/libstdc++-v3/testsuite"
util_dir="${test_dir}/util"
gcc_dir="${test_dir}/special_functions"
ext_dir="${test_dir}/ext/special_functions"
tr1_dir="${test_dir}/tr1/5_numerical_facilities/special_functions"

${copy} check/check_tr1_assoc_laguerre.cc ${tr1_dir}/01_assoc_laguerre/check_value.cc
${copy} check/check_tr1_assoc_legendre.cc ${tr1_dir}/02_assoc_legendre/check_value.cc
${copy} check/check_tr1_beta.cc 	  ${tr1_dir}/03_beta/check_value.cc
${copy} check/check_tr1_comp_ellint_1.cc  ${tr1_dir}/04_comp_ellint_1/check_value.cc
${copy} check/check_tr1_comp_ellint_2.cc  ${tr1_dir}/05_comp_ellint_2/check_value.cc
${copy} check/check_tr1_comp_ellint_3.cc  ${tr1_dir}/06_comp_ellint_3/check_value.cc
${copy} check/check_tr1_conf_hyperg.cc    ${tr1_dir}/07_conf_hyperg/check_value.cc
${copy} check/check_tr1_cyl_bessel_i.cc   ${tr1_dir}/08_cyl_bessel_i/check_value.cc
${copy} check/check_tr1_cyl_bessel_j.cc   ${tr1_dir}/09_cyl_bessel_j/check_value.cc
${copy} check/check_tr1_cyl_bessel_k.cc   ${tr1_dir}/10_cyl_bessel_k/check_value.cc
${copy} check/check_tr1_cyl_neumann.cc    ${tr1_dir}/11_cyl_neumann/check_value.cc
${copy} check/check_tr1_ellint_1.cc	  ${tr1_dir}/12_ellint_1/check_value.cc
${copy} check/check_tr1_ellint_2.cc	  ${tr1_dir}/13_ellint_2/check_value.cc
${copy} check/check_tr1_ellint_3.cc	  ${tr1_dir}/14_ellint_3/check_value.cc
${copy} check/check_tr1_expint.cc	  ${tr1_dir}/15_expint/check_value_neg.cc
${copy} check/check_tr1_hermite.cc	  ${tr1_dir}/16_hermite/check_value.cc
${copy} check/check_tr1_hyperg.cc         ${tr1_dir}/17_hyperg/check_value.cc
${copy} check/check_tr1_laguerre.cc	  ${tr1_dir}/18_laguerre/check_value.cc
${copy} check/check_tr1_legendre.cc	  ${tr1_dir}/19_legendre/check_value.cc
${copy} check/check_tr1_riemann_zeta.cc   ${tr1_dir}/20_riemann_zeta/check_value_neg.cc
${copy} check/check_tr1_sph_bessel.cc	  ${tr1_dir}/21_sph_bessel/check_value.cc
${copy} check/check_tr1_sph_legendre.cc   ${tr1_dir}/22_sph_legendre/check_value.cc
${copy} check/check_tr1_sph_neumann.cc    ${tr1_dir}/23_sph_neumann/check_value.cc

${copy} check/check_assoc_laguerre.cc ${gcc_dir}/01_assoc_laguerre/check_value.cc
${copy} check/check_assoc_legendre.cc ${gcc_dir}/02_assoc_legendre/check_value.cc
${copy} check/check_beta.cc           ${gcc_dir}/03_beta/check_value.cc
${copy} check/check_comp_ellint_1.cc  ${gcc_dir}/04_comp_ellint_1/check_value.cc
${copy} check/check_comp_ellint_2.cc  ${gcc_dir}/05_comp_ellint_2/check_value.cc
${copy} check/check_comp_ellint_3.cc  ${gcc_dir}/06_comp_ellint_3/check_value.cc
${copy} check/check_cyl_bessel_i.cc   ${gcc_dir}/07_cyl_bessel_i/check_value.cc
${copy} check/check_cyl_bessel_j.cc   ${gcc_dir}/08_cyl_bessel_j/check_value.cc
${copy} check/check_cyl_bessel_k.cc   ${gcc_dir}/09_cyl_bessel_k/check_value.cc
${copy} check/check_cyl_neumann.cc    ${gcc_dir}/10_cyl_neumann/check_value.cc
${copy} check/check_ellint_1.cc       ${gcc_dir}/11_ellint_1/check_value.cc
${copy} check/check_ellint_2.cc       ${gcc_dir}/12_ellint_2/check_value.cc
${copy} check/check_ellint_3.cc       ${gcc_dir}/13_ellint_3/check_value.cc
${copy} check/check_expint.cc         ${gcc_dir}/14_expint/check_value.cc
${copy} check/check_hermite.cc        ${gcc_dir}/15_hermite/check_value.cc
${copy} check/check_laguerre.cc       ${gcc_dir}/16_laguerre/check_value.cc
${copy} check/check_legendre.cc       ${gcc_dir}/17_legendre/check_value.cc
${copy} check/check_riemann_zeta.cc   ${gcc_dir}/18_riemann_zeta/check_value.cc
${copy} check/check_sph_bessel.cc     ${gcc_dir}/19_sph_bessel/check_value.cc
${copy} check/check_sph_legendre.cc   ${gcc_dir}/20_sph_legendre/check_value.cc
${copy} check/check_sph_neumann.cc    ${gcc_dir}/21_sph_neumann/check_value.cc

${copy} check/check_conf_hyperg.cc    ${ext_dir}/conf_hyperg/check_value.cc
${copy} check/check_hyperg.cc         ${ext_dir}/hyperg/check_value.cc
