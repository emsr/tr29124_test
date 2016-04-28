#!  /bin/bash

copy="cp -f"
makedir="mkdir -p"

rm -rf "testsuite"

utildir="testsuite/util"
test_dir="testsuite/special_functions"
text_dir="testsuite/ext/special_functions"


${makedir} ${test_dir}/01_assoc_laguerre
${makedir} ${test_dir}/02_assoc_legendre
${makedir} ${test_dir}/03_beta
${makedir} ${test_dir}/04_comp_ellint_1
${makedir} ${test_dir}/05_comp_ellint_2
${makedir} ${test_dir}/06_comp_ellint_3
${makedir} ${test_dir}/07_cyl_bessel_i
${makedir} ${test_dir}/08_cyl_bessel_j
${makedir} ${test_dir}/09_cyl_bessel_k
${makedir} ${test_dir}/10_cyl_neumann
${makedir} ${test_dir}/11_ellint_1
${makedir} ${test_dir}/12_ellint_2
${makedir} ${test_dir}/13_ellint_3
${makedir} ${test_dir}/14_expint
${makedir} ${test_dir}/15_hermite
${makedir} ${test_dir}/16_laguerre
${makedir} ${test_dir}/17_legendre
${makedir} ${test_dir}/18_riemann_zeta
${makedir} ${test_dir}/19_sph_bessel
${makedir} ${test_dir}/20_sph_legendre
${makedir} ${test_dir}/21_sph_neumann


${makedir} ${utildir}

${copy} specfun_testcase.h        ${utildir}

${copy} check_assoc_laguerre.cc   ${test_dir}/01_assoc_laguerre/check_value.cc
${copy} check_assoc_legendre.cc   ${test_dir}/02_assoc_legendre/check_value.cc
${copy} check_beta.cc             ${test_dir}/03_beta/check_value.cc
${copy} check_comp_ellint_1.cc    ${test_dir}/04_comp_ellint_1/check_value.cc
${copy} check_comp_ellint_2.cc    ${test_dir}/05_comp_ellint_2/check_value.cc
${copy} check_comp_ellint_3.cc    ${test_dir}/06_comp_ellint_3/check_value.cc
${copy} check_cyl_bessel_i.cc     ${test_dir}/07_cyl_bessel_i/check_value.cc
${copy} pr56216_cyl_bessel_i.cc   ${test_dir}/07_cyl_bessel_i/pr56216.cc
${copy} check_cyl_bessel_j.cc     ${test_dir}/08_cyl_bessel_j/check_value.cc
${copy} check_cyl_bessel_k.cc     ${test_dir}/09_cyl_bessel_k/check_value.cc
${copy} check_cyl_neumann.cc      ${test_dir}/10_cyl_neumann/check_value.cc
${copy} check_ellint_1.cc         ${test_dir}/11_ellint_1/check_value.cc
${copy} check_ellint_2.cc         ${test_dir}/12_ellint_2/check_value.cc
${copy} check_ellint_3.cc         ${test_dir}/13_ellint_3/check_value.cc
${copy} check_expint.cc           ${test_dir}/14_expint/check_value.cc
${copy} check_hermite.cc          ${test_dir}/15_hermite/check_value.cc
${copy} check_laguerre.cc         ${test_dir}/16_laguerre/check_value.cc
${copy} check_legendre.cc         ${test_dir}/17_legendre/check_value.cc
${copy} check_riemann_zeta.cc     ${test_dir}/18_riemann_zeta/check_value.cc
${copy} check_sph_bessel.cc       ${test_dir}/19_sph_bessel/check_value.cc
${copy} check_sph_legendre.cc     ${test_dir}/20_sph_legendre/check_value.cc
${copy} check_sph_neumann.cc      ${test_dir}/21_sph_neumann/check_value.cc


${makedir} ${text_dir}/airy_ai
${makedir} ${text_dir}/airy_bi
${makedir} ${text_dir}/bincoef
${makedir} ${text_dir}/clausen
${makedir} ${text_dir}/comp_ellint_d
${makedir} ${text_dir}/conf_hyperg
${makedir} ${text_dir}/conf_hyperg_lim
${makedir} ${text_dir}/coshint
${makedir} ${text_dir}/cosint
${makedir} ${text_dir}/cyl_hankel_1
${makedir} ${text_dir}/cyl_hankel_2
${makedir} ${text_dir}/dawson
${makedir} ${text_dir}/dilog
${makedir} ${text_dir}/dirichlet_eta
${makedir} ${text_dir}/double_factorial
${makedir} ${text_dir}/ellint_d
${makedir} ${text_dir}/ellint_rc
${makedir} ${text_dir}/ellint_rd
${makedir} ${text_dir}/ellint_rf
${makedir} ${text_dir}/ellint_rg
${makedir} ${text_dir}/ellint_rj
${makedir} ${text_dir}/expint_e1
${makedir} ${text_dir}/ellnome
${makedir} ${text_dir}/factorial
${makedir} ${text_dir}/fresnel_c
${makedir} ${text_dir}/fresnel_s
${makedir} ${text_dir}/gamma_l
${makedir} ${text_dir}/gamma_u
${makedir} ${text_dir}/gegenbauer
${makedir} ${text_dir}/heuman_lambda
${makedir} ${text_dir}/hurwitz_zeta
${makedir} ${text_dir}/hyperg
${makedir} ${text_dir}/ibeta
${makedir} ${text_dir}/jacobi_sn
${makedir} ${text_dir}/jacobi_cn
${makedir} ${text_dir}/jacobi_dn
${makedir} ${text_dir}/lbincoef
${makedir} ${text_dir}/ldouble_factorial
${makedir} ${text_dir}/legendre_q
${makedir} ${text_dir}/lfactorial
${makedir} ${text_dir}/lpochhammer_l
${makedir} ${text_dir}/lpochhammer_u
${makedir} ${text_dir}/owens_t
${makedir} ${text_dir}/pochhammer_l
${makedir} ${text_dir}/pochhammer_u
${makedir} ${text_dir}/psi
${makedir} ${text_dir}/radpoly
${makedir} ${text_dir}/sinc
${makedir} ${text_dir}/sinc_pi
${makedir} ${text_dir}/sinhint
${makedir} ${text_dir}/sinint
${makedir} ${text_dir}/sph_bessel_i
${makedir} ${text_dir}/sph_bessel_k
${makedir} ${text_dir}/sph_hankel_1
${makedir} ${text_dir}/sph_hankel_2
${makedir} ${text_dir}/sph_harmonic
${makedir} ${text_dir}/theta_1
${makedir} ${text_dir}/theta_2
${makedir} ${text_dir}/theta_3
${makedir} ${text_dir}/theta_4
${makedir} ${text_dir}/zernike

${copy} check_airy_ai.cc            ${text_dir}/airy_ai/check_value.cc
${copy} check_airy_bi.cc            ${text_dir}/airy_bi/check_value.cc
${copy} check_bincoef.cc            ${text_dir}/bincoef/check_value.cc
${copy} check_chi.cc                ${text_dir}/coshint/check_chi.cc
${copy} check_clausen_c.cc          ${text_dir}/clausen/check_value.cc
${copy} check_comp_ellint_d.cc      ${text_dir}/comp_ellint_d/check_value.cc
${copy} check_conf_hyperg.cc        ${text_dir}/conf_hyperg/check_value.cc
${copy} check_conf_hyperg_lim.cc    ${text_dir}/conf_hyperg_lim/check_value.cc
${copy} check_coshint.cc            ${text_dir}/coshint/check_value.cc
${copy} check_cosint.cc             ${text_dir}/cosint/check_value.cc
${copy} check_cyl_hankel_1.cc       ${text_dir}/cyl_hankel_1/check_value.cc
${copy} pr56216_cyl_hankel_1.cc     ${test_dir}/cyl_hankel_1/pr56216.cc
${copy} check_cyl_hankel_2.cc       ${text_dir}/cyl_hankel_2/check_value.cc
${copy} pr56216_cyl_hankel_2.cc     ${test_dir}/cyl_hankel_2/pr56216.cc
${copy} check_dawson.cc             ${text_dir}/dawson/check_value.cc
${copy} check_dilog.cc              ${text_dir}/dilog/check_value.cc
${copy} check_dirichlet_eta.cc      ${text_dir}/dirichlet_eta/check_value.cc
${copy} check_double_factorial.cc   ${text_dir}/double_factorial/check_value.cc
${copy} check_ellint_d.cc           ${text_dir}/ellint_d/check_value.cc
${copy} check_ellint_rc.cc          ${text_dir}/ellint_rc/check_value.cc
${copy} check_ellint_rd.cc          ${text_dir}/ellint_rd/check_value.cc
${copy} check_ellint_rf.cc          ${text_dir}/ellint_rf/check_value.cc
${copy} check_ellint_rg.cc          ${text_dir}/ellint_rg/check_value.cc
${copy} check_ellint_rj.cc          ${text_dir}/ellint_rj/check_value.cc
${copy} check_ellnome.cc            ${text_dir}/ellnome/check_value.cc
${copy} check_expint_e1.cc          ${text_dir}/expint_e1/check_value.cc
${copy} check_factorial.cc          ${text_dir}/factorial/check_value.cc
${copy} check_fresnel_c.cc          ${text_dir}/fresnel_c/check_value.cc
${copy} check_fresnel_s.cc          ${text_dir}/fresnel_s/check_value.cc
${copy} check_gamma_l.cc            ${text_dir}/gamma_l/check_value.cc
${copy} check_gamma_u.cc            ${text_dir}/gamma_u/check_value.cc
${copy} check_gegenbauer.cc         ${text_dir}/gegenbauer/check_value.cc
${copy} check_heuman_lambda.cc      ${text_dir}/heuman_lambda/check_value.cc
${copy} check_hurwitz_zeta.cc       ${text_dir}/hurwitz_zeta/check_value.cc
${copy} check_hyperg.cc             ${text_dir}/hyperg/check_value.cc
${copy} check_ibeta.cc              ${text_dir}/ibeta/check_value.cc
${copy} check_jacobi.cc             ${text_dir}/jacobi/check_value.cc
${copy} check_jacobi_sn.cc          ${text_dir}/jacobi_sn/check_value.cc
${copy} check_jacobi_cn.cc          ${text_dir}/jacobi_cn/check_value.cc
${copy} check_jacobi_dn.cc          ${text_dir}/jacobi_dn/check_value.cc
${copy} check_lbincoef.cc           ${text_dir}/lbincoef/check_value.cc
${copy} check_ldouble_factorial.cc  ${text_dir}/ldouble_factorial/check_value.cc
${copy} check_legendre_q.cc         ${text_dir}/legendre_q/check_value.cc
${copy} check_lfactorial.cc         ${text_dir}/lfactorial/check_value.cc
${copy} check_lpochhammer_l.cc      ${text_dir}/lpochhammer_l/check_value.cc
${copy} check_lpochhammer_u.cc      ${text_dir}/lpochhammer_u/check_value.cc
${copy} check_owens_t.cc            ${text_dir}/owens_t/check_value.cc
${copy} check_pochhammer_l.cc       ${text_dir}/pochhammer_l/check_value.cc
${copy} check_pochhammer_u.cc       ${text_dir}/pochhammer_u/check_value.cc
${copy} check_radpoly.cc            ${text_dir}/radpoly/check_value.cc
${copy} check_psi.cc                ${text_dir}/psi/check_value.cc
${copy} check_shi.cc                ${text_dir}/sinhint/check_shi.cc
${copy} check_sinc.cc               ${text_dir}/sinc/check_value.cc
${copy} check_sinc_pi.cc            ${text_dir}/sinc_pi/check_value.cc
${copy} check_sinhint.cc            ${text_dir}/sinhint/check_value.cc
${copy} check_sinint.cc             ${text_dir}/sinint/check_value.cc
${copy} check_sph_bessel_i.cc       ${text_dir}/sph_bessel_i/check_value.cc
${copy} check_sph_bessel_k.cc       ${text_dir}/sph_bessel_k/check_value.cc
${copy} check_sph_hankel_1.cc       ${text_dir}/sph_hankel_1/check_value.cc
${copy} check_sph_hankel_2.cc       ${text_dir}/sph_hankel_2/check_value.cc
${copy} check_sph_harmonic.cc       ${text_dir}/sph_harmonic/check_value.cc
${copy} check_theta_1.cc            ${text_dir}/theta_1/check_value.cc
${copy} check_theta_2.cc            ${text_dir}/theta_2/check_value.cc
${copy} check_theta_3.cc            ${text_dir}/theta_3/check_value.cc
${copy} check_theta_4.cc            ${text_dir}/theta_4/check_value.cc
${copy} check_zernike.cc            ${text_dir}/zernike/check_value.cc

