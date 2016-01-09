#!  /bin/bash

copy="cp -f"

base_dir="$HOME/gcc_tr29124"
test_dir="${base_dir}/libstdc++-v3/testsuite"
util_dir="${test_dir}/util"
gcc_dir="${test_dir}/special_functions"
ext_dir="${test_dir}/ext/special_functions"

${copy} check_assoc_laguerre.cc ${gcc_dir}/01_assoc_laguerre/check_value.cc
${copy} check_assoc_legendre.cc ${gcc_dir}/02_assoc_legendre/check_value.cc
${copy} check_beta.cc           ${gcc_dir}/03_beta/check_value.cc
${copy} check_comp_ellint_1.cc  ${gcc_dir}/04_comp_ellint_1/check_value.cc
${copy} check_comp_ellint_2.cc  ${gcc_dir}/05_comp_ellint_2/check_value.cc
${copy} check_comp_ellint_3.cc  ${gcc_dir}/06_comp_ellint_3/check_value.cc
${copy} check_cyl_bessel_i.cc   ${gcc_dir}/07_cyl_bessel_i/check_value.cc
${copy} check_cyl_bessel_j.cc   ${gcc_dir}/08_cyl_bessel_j/check_value.cc
${copy} check_cyl_bessel_k.cc   ${gcc_dir}/09_cyl_bessel_k/check_value.cc
${copy} check_cyl_neumann.cc    ${gcc_dir}/10_cyl_neumann/check_value.cc
${copy} check_ellint_1.cc       ${gcc_dir}/11_ellint_1/check_value.cc
${copy} check_ellint_2.cc       ${gcc_dir}/12_ellint_2/check_value.cc
${copy} check_ellint_3.cc       ${gcc_dir}/13_ellint_3/check_value.cc
${copy} check_expint.cc         ${gcc_dir}/14_expint/check_value.cc
${copy} check_hermite.cc        ${gcc_dir}/15_hermite/check_value.cc
${copy} check_laguerre.cc       ${gcc_dir}/16_laguerre/check_value.cc
${copy} check_legendre.cc       ${gcc_dir}/17_legendre/check_value.cc
${copy} check_riemann_zeta.cc   ${gcc_dir}/18_riemann_zeta/check_value.cc
${copy} check_sph_bessel.cc     ${gcc_dir}/19_sph_bessel/check_value.cc
${copy} check_sph_legendre.cc   ${gcc_dir}/20_sph_legendre/check_value.cc
${copy} check_sph_neumann.cc    ${gcc_dir}/21_sph_neumann/check_value.cc


${copy} check_conf_hyperg.cc    ${ext_dir}/conf_hyperg/check_value.cc
${copy} check_hyperg.cc         ${ext_dir}/hyperg/check_value.cc
