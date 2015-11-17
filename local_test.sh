#!  /bin/bash

copy="cp -f"
makedir="mkdir -p"

rm -rf "testsuite"

utildir="testsuite/util"
testdir="testsuite/special_functions"
textdir="testsuite/ext/special_functions"

${makedir} ${testdir}/01_assoc_laguerre
${makedir} ${testdir}/02_assoc_legendre
${makedir} ${testdir}/03_beta
${makedir} ${testdir}/04_comp_ellint_1
${makedir} ${testdir}/05_comp_ellint_2
${makedir} ${testdir}/06_comp_ellint_3
${makedir} ${testdir}/07_cyl_bessel_i
${makedir} ${testdir}/08_cyl_bessel_j
${makedir} ${testdir}/09_cyl_bessel_k
${makedir} ${testdir}/10_cyl_neumann
${makedir} ${testdir}/11_ellint_1
${makedir} ${testdir}/12_ellint_2
${makedir} ${testdir}/13_ellint_3
${makedir} ${testdir}/14_expint
${makedir} ${testdir}/15_hermite
${makedir} ${testdir}/16_laguerre
${makedir} ${testdir}/17_legendre
${makedir} ${testdir}/18_riemann_zeta
${makedir} ${testdir}/19_sph_bessel
${makedir} ${testdir}/20_sph_legendre
${makedir} ${testdir}/21_sph_neumann

${makedir} ${utildir}
${copy} specfun_testcase.h        ${utildir}

${copy} check_assoc_laguerre.cc   ${testdir}/01_assoc_laguerre/check_value.cc
${copy} check_assoc_legendre.cc   ${testdir}/02_assoc_legendre/check_value.cc
${copy} check_beta.cc             ${testdir}/03_beta/check_value.cc
${copy} check_comp_ellint_1.cc    ${testdir}/04_comp_ellint_1/check_value.cc
${copy} check_comp_ellint_2.cc    ${testdir}/05_comp_ellint_2/check_value.cc
${copy} check_comp_ellint_3.cc    ${testdir}/06_comp_ellint_3/check_value.cc
${copy} check_cyl_bessel_i.cc     ${testdir}/07_cyl_bessel_i/check_value.cc
${copy} check_cyl_bessel_j.cc     ${testdir}/08_cyl_bessel_j/check_value.cc
${copy} check_cyl_bessel_k.cc     ${testdir}/09_cyl_bessel_k/check_value.cc
${copy} check_cyl_neumann.cc      ${testdir}/10_cyl_neumann/check_value.cc
${copy} check_ellint_1.cc         ${testdir}/11_ellint_1/check_value.cc
${copy} check_ellint_2.cc         ${testdir}/12_ellint_2/check_value.cc
${copy} check_ellint_3.cc         ${testdir}/13_ellint_3/check_value.cc
${copy} check_expint.cc           ${testdir}/14_expint/check_value.cc

${copy} check_laguerre.cc         ${testdir}/16_laguerre/check_value.cc
${copy} check_legendre.cc         ${testdir}/17_legendre/check_value.cc
${copy} check_riemann_zeta.cc     ${testdir}/18_riemann_zeta/check_value.cc
${copy} check_sph_bessel.cc       ${testdir}/19_sph_bessel/check_value.cc
${copy} check_sph_legendre.cc     ${testdir}/20_sph_legendre/check_value.cc
${copy} check_sph_neumann.cc      ${testdir}/21_sph_neumann/check_value.cc

${makedir} ${textdir}/airy
${makedir} ${textdir}/conf_hyperg
${makedir} ${textdir}/hurwitz_zeta
${makedir} ${textdir}/hyperg

${copy} check_airy.cc             ${textdir}/airy/check_value.cc
${copy} check_conf_hyperg.cc      ${textdir}/conf_hyperg/check_value.cc
${copy} check_hurwitz_zeta.cc     ${textdir}/hurwitz_zeta/check_value.cc
${copy} check_hyperg.cc           ${textdir}/hyperg/check_value.cc
