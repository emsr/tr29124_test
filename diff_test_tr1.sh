#!/bin/bash

test_dir="$HOME/tr29124_test/testsuite/special_functions"
tr1_dir="$HOME/gcc_tr29124/libstdc++-v3/testsuite/tr1/5_numerical_facilities/special_functions"

tool="kdiff3"

${tool} ${test_dir}/01_assoc_laguerre/ ${tr1_dir}/01_assoc_laguerre/
${tool} ${test_dir}/02_assoc_legendre/ ${tr1_dir}/02_assoc_legendre/
${tool} ${test_dir}/03_beta/           ${tr1_dir}/03_beta/
${tool} ${test_dir}/04_comp_ellint_1/  ${tr1_dir}/04_comp_ellint_1/
${tool} ${test_dir}/05_comp_ellint_2/  ${tr1_dir}/05_comp_ellint_2/
${tool} ${test_dir}/06_comp_ellint_3/  ${tr1_dir}/06_comp_ellint_3/
${tool} ${test_dir}/07_cyl_bessel_i/   ${tr1_dir}/08_cyl_bessel_i/
${tool} ${test_dir}/08_cyl_bessel_j/   ${tr1_dir}/09_cyl_bessel_j/
${tool} ${test_dir}/09_cyl_bessel_k/   ${tr1_dir}/10_cyl_bessel_k/
${tool} ${test_dir}/10_cyl_neumann/    ${tr1_dir}/11_cyl_neumann/
${tool} ${test_dir}/11_ellint_1/       ${tr1_dir}/12_ellint_1/
${tool} ${test_dir}/12_ellint_2/       ${tr1_dir}/13_ellint_2/
${tool} ${test_dir}/13_ellint_3/       ${tr1_dir}/14_ellint_3/
${tool} ${test_dir}/14_expint/         ${tr1_dir}/15_expint/
${tool} ${test_dir}/15_hermite/        ${tr1_dir}/16_hermite/
${tool} ${test_dir}/16_laguerre/       ${tr1_dir}/18_laguerre/
${tool} ${test_dir}/17_legendre/       ${tr1_dir}/19_legendre/
${tool} ${test_dir}/18_riemann_zeta/   ${tr1_dir}/20_riemann_zeta/
${tool} ${test_dir}/19_sph_bessel/     ${tr1_dir}/21_sph_bessel/
${tool} ${test_dir}/20_sph_legendre/   ${tr1_dir}/22_sph_legendre/
${tool} ${test_dir}/21_sph_neumann/    ${tr1_dir}/23_sph_neumann/
