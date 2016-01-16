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
${copy} check_hermite.cc          ${testdir}/15_hermite/check_value.cc
${copy} check_laguerre.cc         ${testdir}/16_laguerre/check_value.cc
${copy} check_legendre.cc         ${testdir}/17_legendre/check_value.cc
${copy} check_riemann_zeta.cc     ${testdir}/18_riemann_zeta/check_value.cc
${copy} check_sph_bessel.cc       ${testdir}/19_sph_bessel/check_value.cc
${copy} check_sph_legendre.cc     ${testdir}/20_sph_legendre/check_value.cc
${copy} check_sph_neumann.cc      ${testdir}/21_sph_neumann/check_value.cc


${makedir} ${textdir}/airy_ai
${makedir} ${textdir}/airy_bi
${makedir} ${textdir}/bincoef
${makedir} ${textdir}/conf_hyperg
${makedir} ${textdir}/conf_hyperg_lim
${makedir} ${textdir}/coshint
${makedir} ${textdir}/cosint
${makedir} ${textdir}/dawson
${makedir} ${textdir}/dilog
${makedir} ${textdir}/double_factorial
${makedir} ${textdir}/ellint_rc
${makedir} ${textdir}/ellint_rd
${makedir} ${textdir}/ellint_rf
${makedir} ${textdir}/ellint_rj
${makedir} ${textdir}/expint_e1
${makedir} ${textdir}/factorial
${makedir} ${textdir}/fresnel_c
${makedir} ${textdir}/fresnel_s
${makedir} ${textdir}/gamma_l
${makedir} ${textdir}/gamma_u
${makedir} ${textdir}/gegenbauer
${makedir} ${textdir}/hurwitz_zeta
${makedir} ${textdir}/hyperg
${makedir} ${textdir}/ibeta
${makedir} ${textdir}/jacobi_sn
${makedir} ${textdir}/jacobi_cn
${makedir} ${textdir}/jacobi_dn
${makedir} ${textdir}/lbincoef
${makedir} ${textdir}/ldouble_factorial
${makedir} ${textdir}/legendre_q
${makedir} ${textdir}/lfactorial
${makedir} ${textdir}/lpochhammer_l
${makedir} ${textdir}/lpochhammer_u
${makedir} ${textdir}/pochhammer_l
${makedir} ${textdir}/pochhammer_u
${makedir} ${textdir}/psi
${makedir} ${textdir}/sinc
${makedir} ${textdir}/sinc_pi
${makedir} ${textdir}/sinhint
${makedir} ${textdir}/sinint
${makedir} ${textdir}/sph_bessel_i
${makedir} ${textdir}/sph_bessel_k

${copy} check_airy_ai.cc            ${textdir}/airy_ai/check_value.cc
${copy} check_airy_bi.cc            ${textdir}/airy_bi/check_value.cc
${copy} check_bincoef.cc            ${textdir}/bincoef/check_value.cc
${copy} check_conf_hyperg.cc        ${textdir}/conf_hyperg/check_value.cc
#${copy} check_conf_hyperg_lim.cc    ${textdir}/conf_hyperg_lim/check_value.cc
${copy} check_coshint.cc            ${textdir}/coshint/check_value.cc
${copy} check_cosint.cc             ${textdir}/cosint/check_value.cc
${copy} check_dawson.cc             ${textdir}/dawson/check_value.cc
${copy} check_dilog.cc              ${textdir}/dilog/check_value.cc
${copy} check_double_factorial.cc   ${textdir}/double_factorial/check_value.cc
${copy} check_ellint_rc.cc          ${textdir}/ellint_rc/check_value.cc
${copy} check_ellint_rd.cc          ${textdir}/ellint_rd/check_value.cc
${copy} check_ellint_rf.cc          ${textdir}/ellint_rf/check_value.cc
${copy} check_ellint_rj.cc          ${textdir}/ellint_rj/check_value.cc
${copy} check_expint_e1.cc          ${textdir}/expint_e1/check_value.cc
${copy} check_factorial.cc          ${textdir}/factorial/check_value.cc
${copy} check_fresnel_c.cc          ${textdir}/fresnel_c/check_value.cc
${copy} check_fresnel_s.cc          ${textdir}/fresnel_s/check_value.cc
#${copy} check_gamma_l.cc            ${textdir}/gamma_l/check_value.cc
${copy} check_gamma_u.cc            ${textdir}/gamma_u/check_value.cc
${copy} check_gegenbauer.cc         ${textdir}/gegenbauer/check_value.cc
${copy} check_hurwitz_zeta.cc       ${textdir}/hurwitz_zeta/check_value.cc
${copy} check_hyperg.cc             ${textdir}/hyperg/check_value.cc
${copy} check_ibeta.cc              ${textdir}/ibeta/check_value.cc
${copy} check_jacobi_sn.cc          ${textdir}/jacobi_sn/check_value.cc
${copy} check_jacobi_cn.cc          ${textdir}/jacobi_cn/check_value.cc
${copy} check_jacobi_dn.cc          ${textdir}/jacobi_dn/check_value.cc
${copy} check_lbincoef.cc           ${textdir}/lbincoef/check_value.cc
${copy} check_ldouble_factorial.cc  ${textdir}/ldouble_factorial/check_value.cc
${copy} check_legendre_q.cc         ${textdir}/legendre_q/check_value.cc
${copy} check_lfactorial.cc         ${textdir}/lfactorial/check_value.cc
#${copy} check_lpochhammer_l.cc      ${textdir}/lpochhammer_l/check_value.cc
${copy} check_lpochhammer_u.cc      ${textdir}/lpochhammer_u/check_value.cc
#${copy} check_pochhammer_l.cc       ${textdir}/pochhammer_l/check_value.cc
${copy} check_pochhammer_u.cc       ${textdir}/pochhammer_u/check_value.cc
${copy} check_psi.cc                ${textdir}/psi/check_value.cc
${copy} check_sinc.cc               ${textdir}/sinc/check_value.cc
${copy} check_sinc_pi.cc            ${textdir}/sinc_pi/check_value.cc
${copy} check_sinhint.cc            ${textdir}/sinhint/check_value.cc
${copy} check_sinint.cc             ${textdir}/sinint/check_value.cc
${copy} check_sph_bessel_i.cc       ${textdir}/sph_bessel_i/check_value.cc
${copy} check_sph_bessel_k.cc       ${textdir}/sph_bessel_k/check_value.cc
