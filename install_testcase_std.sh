#!  /bin/bash

tool="cp -f"
makedir="mkdir -p"

suffix="_tr29124"
if [ $# -ge 1 ]; then
  suffix="$1"
fi

base_dir="$HOME/gcc${suffix}"
test_dir="${base_dir}/libstdc++-v3/testsuite"
util_dir="${test_dir}/util"
gcc_dir="${test_dir}/special_functions"
ext_dir="${test_dir}/ext/special_functions"


${makedir} ${utildir}

${tool} specfun_testcase.h              ${util_dir}


${makedir} ${gcc_dir}/01_assoc_laguerre
${makedir} ${gcc_dir}/02_assoc_legendre
${makedir} ${gcc_dir}/03_beta
${makedir} ${gcc_dir}/04_comp_ellint_1
${makedir} ${gcc_dir}/05_comp_ellint_2
${makedir} ${gcc_dir}/06_comp_ellint_3
${makedir} ${gcc_dir}/07_cyl_bessel_i
${makedir} ${gcc_dir}/08_cyl_bessel_j
${makedir} ${gcc_dir}/09_cyl_bessel_k
${makedir} ${gcc_dir}/10_cyl_neumann
${makedir} ${gcc_dir}/11_ellint_1
${makedir} ${gcc_dir}/12_ellint_2
${makedir} ${gcc_dir}/13_ellint_3
${makedir} ${gcc_dir}/14_expint
${makedir} ${gcc_dir}/15_hermite
${makedir} ${gcc_dir}/16_laguerre
${makedir} ${gcc_dir}/17_legendre
${makedir} ${gcc_dir}/18_riemann_zeta
${makedir} ${gcc_dir}/19_sph_bessel
${makedir} ${gcc_dir}/20_sph_legendre
${makedir} ${gcc_dir}/21_sph_neumann

${tool} check/check_assoc_laguerre.cc   ${gcc_dir}/01_assoc_laguerre/check_value.cc
${tool} check/check_assoc_legendre.cc   ${gcc_dir}/02_assoc_legendre/check_value.cc
${tool} check/check_beta.cc             ${gcc_dir}/03_beta/check_value.cc
${tool} check/check_comp_ellint_1.cc    ${gcc_dir}/04_comp_ellint_1/check_value.cc
${tool} check/check_comp_ellint_2.cc    ${gcc_dir}/05_comp_ellint_2/check_value.cc
${tool} check/check_comp_ellint_3.cc    ${gcc_dir}/06_comp_ellint_3/check_value.cc
${tool} check/check_cyl_bessel_i.cc     ${gcc_dir}/07_cyl_bessel_i/check_value.cc
${tool} check/pr56216_cyl_bessel_i.cc   ${gcc_dir}/07_cyl_bessel_i/pr56216.cc
${tool} check/check_cyl_bessel_j.cc     ${gcc_dir}/08_cyl_bessel_j/check_value.cc
${tool} check/origin_bessel_j.cc        ${gcc_dir}/08_cyl_bessel_j/check_origin.cc
${tool} check/check_cyl_bessel_k.cc     ${gcc_dir}/09_cyl_bessel_k/check_value.cc
${tool} check/check_cyl_neumann.cc      ${gcc_dir}/10_cyl_neumann/check_value.cc
${tool} check/origin_cyl_neumann.cc     ${gcc_dir}/10_cyl_neumann/check_origin.cc
${tool} check/check_ellint_1.cc         ${gcc_dir}/11_ellint_1/check_value.cc
${tool} check/check_ellint_2.cc         ${gcc_dir}/12_ellint_2/check_value.cc
${tool} check/check_ellint_3.cc         ${gcc_dir}/13_ellint_3/check_value.cc
${tool} check/check_expint.cc           ${gcc_dir}/14_expint/check_value.cc
${tool} check/check_hermite.cc          ${gcc_dir}/15_hermite/check_value.cc
${tool} check/check_laguerre.cc         ${gcc_dir}/16_laguerre/check_value.cc
${tool} check/check_legendre.cc         ${gcc_dir}/17_legendre/check_value.cc
${tool} check/check_riemann_zeta.cc     ${gcc_dir}/18_riemann_zeta/check_value.cc
${tool} check/check_sph_bessel.cc       ${gcc_dir}/19_sph_bessel/check_value.cc
${tool} check/check_sph_legendre.cc     ${gcc_dir}/20_sph_legendre/check_value.cc
${tool} check/check_sph_neumann.cc      ${gcc_dir}/21_sph_neumann/check_value.cc


${makedir} ${ext_dir}/airy_ai
${makedir} ${ext_dir}/airy_bi
${makedir} ${ext_dir}/bernoulli
${makedir} ${ext_dir}/bincoef
${makedir} ${ext_dir}/clausen
${makedir} ${ext_dir}/comp_ellint_d
${makedir} ${ext_dir}/conf_hyperg
${makedir} ${ext_dir}/conf_hyperg_lim
${makedir} ${ext_dir}/coshint
${makedir} ${ext_dir}/cosint
${makedir} ${ext_dir}/cos_pi
${makedir} ${ext_dir}/cyl_hankel_1
${makedir} ${ext_dir}/cyl_hankel_2
${makedir} ${ext_dir}/dawson
${makedir} ${ext_dir}/dilog
${makedir} ${ext_dir}/dirichlet_beta
${makedir} ${ext_dir}/dirichlet_eta
${makedir} ${ext_dir}/dirichlet_lambda
${makedir} ${ext_dir}/double_factorial
${makedir} ${ext_dir}/ellint_d
${makedir} ${ext_dir}/ellint_rc
${makedir} ${ext_dir}/ellint_rd
${makedir} ${ext_dir}/ellint_rf
${makedir} ${ext_dir}/ellint_rg
${makedir} ${ext_dir}/ellint_rj
${makedir} ${ext_dir}/expint
${makedir} ${ext_dir}/ellnome
${makedir} ${ext_dir}/factorial
${makedir} ${ext_dir}/fresnel_c
${makedir} ${ext_dir}/fresnel_s
${makedir} ${ext_dir}/gegenbauer
${makedir} ${ext_dir}/heuman_lambda
${makedir} ${ext_dir}/hurwitz_zeta
${makedir} ${ext_dir}/hyperg
${makedir} ${ext_dir}/ibeta
${makedir} ${ext_dir}/jacobi
${makedir} ${ext_dir}/jacobi_sn
${makedir} ${ext_dir}/jacobi_cn
${makedir} ${ext_dir}/jacobi_dn
${makedir} ${ext_dir}/lbincoef
${makedir} ${ext_dir}/ldouble_factorial
${makedir} ${ext_dir}/legendre_q
${makedir} ${ext_dir}/lfactorial
${makedir} ${ext_dir}/lgamma
${makedir} ${ext_dir}/lpochhammer_lower
${makedir} ${ext_dir}/lpochhammer
${makedir} ${ext_dir}/owens_t
${makedir} ${ext_dir}/pgamma
${makedir} ${ext_dir}/pochhammer_lower
${makedir} ${ext_dir}/pochhammer
${makedir} ${ext_dir}/psi
${makedir} ${ext_dir}/qgamma
${makedir} ${ext_dir}/radpoly
${makedir} ${ext_dir}/sinc
${makedir} ${ext_dir}/sinc_pi
${makedir} ${ext_dir}/sinhint
${makedir} ${ext_dir}/sinint
${makedir} ${ext_dir}/sin_pi
${makedir} ${ext_dir}/sph_bessel_i
${makedir} ${ext_dir}/sph_bessel_k
${makedir} ${ext_dir}/sph_hankel_1
${makedir} ${ext_dir}/sph_hankel_2
${makedir} ${ext_dir}/sph_harmonic
${makedir} ${ext_dir}/tgamma
${makedir} ${ext_dir}/tgamma_lower
${makedir} ${ext_dir}/theta_1
${makedir} ${ext_dir}/theta_2
${makedir} ${ext_dir}/theta_3
${makedir} ${ext_dir}/theta_4
${makedir} ${ext_dir}/zernike

${tool} check/check_airy_ai.cc            ${ext_dir}/airy_ai/check_value.cc
${tool} check/check_airy_bi.cc            ${ext_dir}/airy_bi/check_value.cc
${tool} check/check_bernoulli.cc          ${ext_dir}/bernoulli/check_value.cc
${tool} check/check_bincoef.cc            ${ext_dir}/bincoef/check_value.cc
${tool} check/check_chi.cc                ${ext_dir}/coshint/check_chi.cc
${tool} check/check_clausen_c.cc          ${ext_dir}/clausen/check_value.cc
${tool} check/check_comp_ellint_d.cc      ${ext_dir}/comp_ellint_d/check_value.cc
${tool} check/check_conf_hyperg.cc        ${ext_dir}/conf_hyperg/check_value.cc
${tool} check/check_conf_hyperg_lim.cc    ${ext_dir}/conf_hyperg_lim/check_value.cc
${tool} check/check_coshint.cc            ${ext_dir}/coshint/check_value.cc
${tool} check/check_cosint.cc             ${ext_dir}/cosint/check_value.cc
${tool} check/check_cyl_hankel_1.cc       ${ext_dir}/cyl_hankel_1/check_value.cc
${tool} check/pr56216_cyl_hankel_1.cc     ${ext_dir}/cyl_hankel_1/pr56216.cc
${tool} check/check_cyl_hankel_2.cc       ${ext_dir}/cyl_hankel_2/check_value.cc
${tool} check/pr56216_cyl_hankel_2.cc     ${ext_dir}/cyl_hankel_2/pr56216.cc
${tool} check/check_dawson.cc             ${ext_dir}/dawson/check_value.cc
${tool} check/check_dilog.cc              ${ext_dir}/dilog/check_value.cc
${tool} check/check_dirichlet_beta.cc     ${ext_dir}/dirichlet_beta/check_value.cc
${tool} check/check_dirichlet_eta.cc      ${ext_dir}/dirichlet_eta/check_value.cc
${tool} check/check_dirichlet_lambda.cc   ${ext_dir}/dirichlet_lambda/check_value.cc
${tool} check/check_double_factorial.cc   ${ext_dir}/double_factorial/check_value.cc
${tool} check/check_ellint_d.cc           ${ext_dir}/ellint_d/check_value.cc
${tool} check/check_ellint_rc.cc          ${ext_dir}/ellint_rc/check_value.cc
${tool} check/check_ellint_rd.cc          ${ext_dir}/ellint_rd/check_value.cc
${tool} check/check_ellint_rf.cc          ${ext_dir}/ellint_rf/check_value.cc
${tool} check/check_ellint_rg.cc          ${ext_dir}/ellint_rg/check_value.cc
${tool} check/check_ellint_rj.cc          ${ext_dir}/ellint_rj/check_value.cc
${tool} check/check_ellnome.cc            ${ext_dir}/ellnome/check_value.cc
${tool} check/check_expint_en.cc          ${ext_dir}/expint/check_value.cc
${tool} check/check_factorial.cc          ${ext_dir}/factorial/check_value.cc
${tool} check/check_fresnel_c.cc          ${ext_dir}/fresnel_c/check_value.cc
${tool} check/check_fresnel_s.cc          ${ext_dir}/fresnel_s/check_value.cc
${tool} check/check_gegenbauer.cc         ${ext_dir}/gegenbauer/check_value.cc
${tool} check/check_heuman_lambda.cc      ${ext_dir}/heuman_lambda/check_value.cc
${tool} check/check_hurwitz_zeta.cc       ${ext_dir}/hurwitz_zeta/check_value.cc
${tool} check/check_hyperg.cc             ${ext_dir}/hyperg/check_value.cc
${tool} check/check_ibeta.cc              ${ext_dir}/ibeta/check_value.cc
${tool} check/check_jacobi.cc             ${ext_dir}/jacobi/check_value.cc
${tool} check/check_jacobi_sn.cc          ${ext_dir}/jacobi_sn/check_value.cc
${tool} check/check_jacobi_cn.cc          ${ext_dir}/jacobi_cn/check_value.cc
${tool} check/check_jacobi_dn.cc          ${ext_dir}/jacobi_dn/check_value.cc
${tool} check/check_lbincoef.cc           ${ext_dir}/lbincoef/check_value.cc
${tool} check/check_ldouble_factorial.cc  ${ext_dir}/ldouble_factorial/check_value.cc
${tool} check/check_legendre_q.cc         ${ext_dir}/legendre_q/check_value.cc
${tool} check/check_lfactorial.cc         ${ext_dir}/lfactorial/check_value.cc
${tool} check/check_lgamma.cc             ${ext_dir}/lgamma/check_value.cc
${tool} check/check_lpochhammer_lower.cc  ${ext_dir}/lpochhammer_lower/check_value.cc
${tool} check/check_lpochhammer.cc        ${ext_dir}/lpochhammer/check_value.cc
${tool} check/check_owens_t.cc            ${ext_dir}/owens_t/check_value.cc
${tool} check/check_pgamma.cc             ${ext_dir}/pgamma/check_value.cc
${tool} check/check_pochhammer_lower.cc   ${ext_dir}/pochhammer_lower/check_value.cc
${tool} check/check_pochhammer.cc         ${ext_dir}/pochhammer/check_value.cc
${tool} check/check_psi.cc                ${ext_dir}/psi/check_value.cc
${tool} check/check_qgamma.cc             ${ext_dir}/qgamma/check_value.cc
${tool} check/check_radpoly.cc            ${ext_dir}/radpoly/check_value.cc
${tool} check/check_shi.cc                ${ext_dir}/sinhint/check_shi.cc
${tool} check/check_sinc.cc               ${ext_dir}/sinc/check_value.cc
${tool} check/check_sinc_pi.cc            ${ext_dir}/sinc_pi/check_value.cc
${tool} check/check_sinhint.cc            ${ext_dir}/sinhint/check_value.cc
${tool} check/check_sinint.cc             ${ext_dir}/sinint/check_value.cc
${tool} check/check_sph_bessel_i.cc       ${ext_dir}/sph_bessel_i/check_value.cc
${tool} check/check_sph_bessel_k.cc       ${ext_dir}/sph_bessel_k/check_value.cc
${tool} check/check_sph_hankel_1.cc       ${ext_dir}/sph_hankel_1/check_value.cc
${tool} check/check_sph_hankel_2.cc       ${ext_dir}/sph_hankel_2/check_value.cc
${tool} check/check_sph_harmonic.cc       ${ext_dir}/sph_harmonic/check_value.cc
${tool} check/check_tgamma.cc             ${ext_dir}/tgamma/check_value.cc
${tool} check/check_tgamma_lower.cc       ${ext_dir}/tgamma_lower/check_value.cc
${tool} check/check_theta_1.cc            ${ext_dir}/theta_1/check_value.cc
${tool} check/check_theta_2.cc            ${ext_dir}/theta_2/check_value.cc
${tool} check/check_theta_3.cc            ${ext_dir}/theta_3/check_value.cc
${tool} check/check_theta_4.cc            ${ext_dir}/theta_4/check_value.cc
${tool} check/check_zernike.cc            ${ext_dir}/zernike/check_value.cc
