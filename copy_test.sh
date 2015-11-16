#!/bin/bash

test_dir="$HOME/tr29124_test/testsuite/special_functions"
gcc_dir="$HOME/gcc_specfun/libstdc++-v3/testsuite/special_functions"

tool="cp -f"

${tool} ${test_dir}/01_assoc_laguerre/* ${gcc_dir}/01_assoc_laguerre/
${tool} ${test_dir}/02_assoc_legendre/* ${gcc_dir}/02_assoc_legendre/
${tool} ${test_dir}/03_beta/*           ${gcc_dir}/03_beta/
${tool} ${test_dir}/04_comp_ellint_1/*  ${gcc_dir}/04_comp_ellint_1/
${tool} ${test_dir}/05_comp_ellint_2/*  ${gcc_dir}/05_comp_ellint_2/
${tool} ${test_dir}/06_comp_ellint_3/*  ${gcc_dir}/06_comp_ellint_3/
${tool} ${test_dir}/07_cyl_bessel_i/*   ${gcc_dir}/07_cyl_bessel_i/
${tool} ${test_dir}/08_cyl_bessel_j/*   ${gcc_dir}/08_cyl_bessel_j/
${tool} ${test_dir}/09_cyl_bessel_k/*   ${gcc_dir}/09_cyl_bessel_k/
${tool} ${test_dir}/10_cyl_neumann/*    ${gcc_dir}/10_cyl_neumann/
${tool} ${test_dir}/11_ellint_1/*       ${gcc_dir}/11_ellint_1/
${tool} ${test_dir}/12_ellint_2/*       ${gcc_dir}/12_ellint_2/
${tool} ${test_dir}/13_ellint_3/*       ${gcc_dir}/13_ellint_3/
${tool} ${test_dir}/14_expint/*	      ${gcc_dir}/14_expint/
${tool} ${test_dir}/15_hermite/*        ${gcc_dir}/15_hermite/
${tool} ${test_dir}/16_laguerre/*       ${gcc_dir}/16_laguerre/
${tool} ${test_dir}/17_legendre/*       ${gcc_dir}/17_legendre/
${tool} ${test_dir}/18_riemann_zeta/*   ${gcc_dir}/18_riemann_zeta/
${tool} ${test_dir}/19_sph_bessel/*     ${gcc_dir}/19_sph_bessel/
${tool} ${test_dir}/20_sph_legendre/*   ${gcc_dir}/20_sph_legendre/
${tool} ${test_dir}/21_sph_neumann/*    ${gcc_dir}/21_sph_neumann/

text_dir="$HOME/tr29124_test/testsuite/ext/special_functions"
ext_dir="$HOME/gcc_specfun/libstdc++-v3/testsuite/ext/special_functions"

${tool} ${text_dir}/airy/*              ${ext_dir}/airy/
${tool} ${text_dir}/conf_hyperg/*       ${ext_dir}/conf_hyperg/
${tool} ${text_dir}/hyperg/*            ${ext_dir}/hyperg/
