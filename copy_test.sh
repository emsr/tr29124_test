#!/bin/bash

test_dir="$HOME/tr29124_test/testsuite/special_functions"
gcc_dir="$HOME/gcc_tr29124/libstdc++-v3/testsuite/special_functions"

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
${tool} ${test_dir}/14_expint/*	        ${gcc_dir}/14_expint/
${tool} ${test_dir}/15_hermite/*        ${gcc_dir}/15_hermite/
${tool} ${test_dir}/16_laguerre/*       ${gcc_dir}/16_laguerre/
${tool} ${test_dir}/17_legendre/*       ${gcc_dir}/17_legendre/
${tool} ${test_dir}/18_riemann_zeta/*   ${gcc_dir}/18_riemann_zeta/
${tool} ${test_dir}/19_sph_bessel/*     ${gcc_dir}/19_sph_bessel/
${tool} ${test_dir}/20_sph_legendre/*   ${gcc_dir}/20_sph_legendre/
${tool} ${test_dir}/21_sph_neumann/*    ${gcc_dir}/21_sph_neumann/

text_dir="$HOME/tr29124_test/testsuite/ext/special_functions"
ext_dir="$HOME/gcc_tr29124/libstdc++-v3/testsuite/ext/special_functions"

${tool} ${text_dir}/airy_ai/*            ${ext_dir}/airy_ai
${tool} ${text_dir}/airy_bi/*            ${ext_dir}/airy_bi
${tool} ${text_dir}/bincoef/*            ${ext_dir}/bincoef
${tool} ${text_dir}/comp_ellint_d/*      ${ext_dir}/comp_ellint_d
${tool} ${text_dir}/conf_hyperg/*        ${ext_dir}/conf_hyperg
${tool} ${text_dir}/conf_hyperg_lim/*    ${ext_dir}/conf_hyperg_lim
${tool} ${text_dir}/coshint/*            ${ext_dir}/coshint
${tool} ${text_dir}/cosint/*             ${ext_dir}/cosint
${tool} ${text_dir}/cyl_hankel_1/*       ${ext_dir}/cyl_hankel_1
${tool} ${text_dir}/cyl_hankel_2/*       ${ext_dir}/cyl_hankel_2
${tool} ${text_dir}/dawson/*             ${ext_dir}/dawson
${tool} ${text_dir}/dirichlet_eta/*      ${ext_dir}/dirichlet_eta
${tool} ${text_dir}/dilog/*              ${ext_dir}/dilog
${tool} ${text_dir}/double_factorial/*   ${ext_dir}/double_factorial
${tool} ${text_dir}/ellint_d/*           ${ext_dir}/ellint_d
${tool} ${text_dir}/ellint_rc/*          ${ext_dir}/ellint_rc
${tool} ${text_dir}/ellint_rd/*          ${ext_dir}/ellint_rd
${tool} ${text_dir}/ellint_rf/*          ${ext_dir}/ellint_rf
${tool} ${text_dir}/ellint_rg/*          ${ext_dir}/ellint_rg
${tool} ${text_dir}/ellint_rj/*          ${ext_dir}/ellint_rj
${tool} ${text_dir}/ellnome/*            ${ext_dir}/ellnome
${tool} ${text_dir}/expint_e1/*          ${ext_dir}/expint_e1
${tool} ${text_dir}/factorial/*          ${ext_dir}/factorial
${tool} ${text_dir}/fresnel_c/*          ${ext_dir}/fresnel_c
${tool} ${text_dir}/fresnel_s/*          ${ext_dir}/fresnel_s
${tool} ${text_dir}/gamma_l/*            ${ext_dir}/gamma_l
${tool} ${text_dir}/gamma_u/*            ${ext_dir}/gamma_u
${tool} ${text_dir}/gegenbauer/*         ${ext_dir}/gegenbauer
${tool} ${text_dir}/hurwitz_zeta/*       ${ext_dir}/hurwitz_zeta
${tool} ${text_dir}/hyperg/*             ${ext_dir}/hyperg
${tool} ${text_dir}/ibeta/*              ${ext_dir}/ibeta
${tool} ${text_dir}/jacobi_sn/*          ${ext_dir}/jacobi_sn
${tool} ${text_dir}/jacobi_cn/*          ${ext_dir}/jacobi_cn
${tool} ${text_dir}/jacobi_dn/*          ${ext_dir}/jacobi_dn
${tool} ${text_dir}/lbincoef/*           ${ext_dir}/lbincoef
${tool} ${text_dir}/ldouble_factorial/*  ${ext_dir}/ldouble_factorial
${tool} ${text_dir}/legendre_q/*         ${ext_dir}/legendre_q
${tool} ${text_dir}/lfactorial/*         ${ext_dir}/lfactorial
${tool} ${text_dir}/lpochhammer_l/*      ${ext_dir}/lpochhammer_l
${tool} ${text_dir}/lpochhammer_u/*      ${ext_dir}/lpochhammer_u
${tool} ${text_dir}/pochhammer_l/*       ${ext_dir}/pochhammer_l
${tool} ${text_dir}/pochhammer_u/*       ${ext_dir}/pochhammer_u
${tool} ${text_dir}/psi/*                ${ext_dir}/psi
${tool} ${text_dir}/sinc/*               ${ext_dir}/sinc
${tool} ${text_dir}/sinc_pi/*            ${ext_dir}/sinc_pi
${tool} ${text_dir}/sinhint/*            ${ext_dir}/sinhint
${tool} ${text_dir}/sinint/*             ${ext_dir}/sinint
${tool} ${text_dir}/sph_bessel_i/*       ${ext_dir}/sph_bessel_i
${tool} ${text_dir}/sph_bessel_k/*       ${ext_dir}/sph_bessel_k
${tool} ${text_dir}/sph_hankel_1/*       ${ext_dir}/sph_hankel_1
${tool} ${text_dir}/sph_hankel_2/*       ${ext_dir}/sph_hankel_2
${tool} ${text_dir}/sph_harmonic/*       ${ext_dir}/sph_harmonic
${tool} ${text_dir}/theta_1/*            ${ext_dir}/theta_1
${tool} ${text_dir}/theta_2/*            ${ext_dir}/theta_2
${tool} ${text_dir}/theta_3/*            ${ext_dir}/theta_3
${tool} ${text_dir}/theta_4/*            ${ext_dir}/theta_4
