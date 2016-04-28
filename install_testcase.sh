#!  /bin/bash

copy="cp -f"

base_dir="$HOME/gcc_specfun"
test_dir="${base_dir}/libstdc++-v3/testsuite"
util_dir="${test_dir}/util"
gcc_dir="${test_dir}/special_functions"
ext_dir="${test_dir}/ext/special_functions"

${copy} check_assoc_laguerre.cc   ${gcc_dir}/01_assoc_laguerre/check_value.cc
${copy} check_assoc_legendre.cc   ${gcc_dir}/02_assoc_legendre/check_value.cc
${copy} check_beta.cc             ${gcc_dir}/03_beta/check_value.cc
${copy} check_comp_ellint_1.cc    ${gcc_dir}/04_comp_ellint_1/check_value.cc
${copy} check_comp_ellint_2.cc    ${gcc_dir}/05_comp_ellint_2/check_value.cc
${copy} check_comp_ellint_3.cc    ${gcc_dir}/06_comp_ellint_3/check_value.cc
${copy} check_cyl_bessel_i.cc     ${gcc_dir}/07_cyl_bessel_i/check_value.cc
${copy} pr56216_cyl_bessel_i.cc   ${gcc_dir}/07_cyl_bessel_i/pr56216.cc
${copy} check_cyl_bessel_j.cc     ${gcc_dir}/08_cyl_bessel_j/check_value.cc
${copy} origin_bessel_j.cc        ${gcc_dir}/08_cyl_bessel_j/check_origin.cc
${copy} check_cyl_bessel_k.cc     ${gcc_dir}/09_cyl_bessel_k/check_value.cc
${copy} check_cyl_neumann.cc      ${gcc_dir}/10_cyl_neumann/check_value.cc
${copy} origin_neumann.cc         ${gcc_dir}/10_cyl_neumann/check_origin.cc
${copy} check_ellint_1.cc         ${gcc_dir}/11_ellint_1/check_value.cc
${copy} check_ellint_2.cc         ${gcc_dir}/12_ellint_2/check_value.cc
${copy} check_ellint_3.cc         ${gcc_dir}/13_ellint_3/check_value.cc
${copy} check_expint.cc           ${gcc_dir}/14_expint/check_value.cc
${copy} check_hermite.cc          ${gcc_dir}/15_hermite/check_value.cc
${copy} check_laguerre.cc         ${gcc_dir}/16_laguerre/check_value.cc
${copy} check_legendre.cc         ${gcc_dir}/17_legendre/check_value.cc
${copy} check_riemann_zeta.cc     ${gcc_dir}/18_riemann_zeta/check_value.cc
${copy} check_sph_bessel.cc       ${gcc_dir}/19_sph_bessel/check_value.cc
${copy} check_sph_legendre.cc     ${gcc_dir}/20_sph_legendre/check_value.cc
${copy} check_sph_neumann.cc      ${gcc_dir}/21_sph_neumann/check_value.cc


${copy} check_airy_ai.cc            ${ext_dir}/airy_ai/check_value.cc
${copy} check_airy_bi.cc            ${ext_dir}/airy_bi/check_value.cc
${copy} check_bincoef.cc            ${ext_dir}/bincoef/check_value.cc
${copy} check_chi.cc                ${ext_dir}/coshint/check_chi.cc
${copy} check_clausen_c.cc          ${ext_dir}/clausen/check_value.cc
${copy} check_comp_ellint_d.cc      ${ext_dir}/comp_ellint_d/check_value.cc
${copy} check_conf_hyperg.cc        ${ext_dir}/conf_hyperg/check_value.cc
${copy} check_conf_hyperg_lim.cc    ${ext_dir}/conf_hyperg_lim/check_value.cc
${copy} check_coshint.cc            ${ext_dir}/coshint/check_value.cc
${copy} check_cosint.cc             ${ext_dir}/cosint/check_value.cc
${copy} check_cyl_hankel_1.cc       ${ext_dir}/cyl_hankel_1/check_value.cc
${copy} pr56216_cyl_hankel_1.cc     ${ext_dir}/cyl_hankel_1/pr56216.cc
${copy} check_cyl_hankel_2.cc       ${ext_dir}/cyl_hankel_2/check_value.cc
${copy} pr56216_cyl_hankel_2.cc     ${ext_dir}/cyl_hankel_2/pr56216.cc
${copy} check_dawson.cc             ${ext_dir}/dawson/check_value.cc
${copy} check_dilog.cc              ${ext_dir}/dilog/check_value.cc
#${copy} check_dirichlet_beta.cc     ${ext_dir}/dirichlet_beta/check_value.cc
${copy} check_dirichlet_eta.cc      ${ext_dir}/dirichlet_eta/check_value.cc
#${copy} check_dirichlet_lambda.cc   ${ext_dir}/dirichlet_lambda/check_value.cc
${copy} check_double_factorial.cc   ${ext_dir}/double_factorial/check_value.cc
${copy} check_ellint_d.cc           ${ext_dir}/ellint_d/check_value.cc
${copy} check_ellint_rc.cc          ${ext_dir}/ellint_rc/check_value.cc
${copy} check_ellint_rd.cc          ${ext_dir}/ellint_rd/check_value.cc
${copy} check_ellint_rf.cc          ${ext_dir}/ellint_rf/check_value.cc
${copy} check_ellint_rg.cc          ${ext_dir}/ellint_rg/check_value.cc
${copy} check_ellint_rj.cc          ${ext_dir}/ellint_rj/check_value.cc
${copy} check_ellnome.cc            ${ext_dir}/ellnome/check_value.cc
${copy} check_expint_e1.cc          ${ext_dir}/expint_e1/check_value.cc
${copy} check_factorial.cc          ${ext_dir}/factorial/check_value.cc
${copy} check_fresnel_c.cc          ${ext_dir}/fresnel_c/check_value.cc
${copy} check_fresnel_s.cc          ${ext_dir}/fresnel_s/check_value.cc
${copy} check_gamma_l.cc            ${ext_dir}/gamma_l/check_value.cc
${copy} check_gamma_u.cc            ${ext_dir}/gamma_u/check_value.cc
${copy} check_gegenbauer.cc         ${ext_dir}/gegenbauer/check_value.cc
${copy} check_heuman_lambda.cc      ${ext_dir}/heuman_lambda/check_value.cc
${copy} check_hurwitz_zeta.cc       ${ext_dir}/hurwitz_zeta/check_value.cc
${copy} check_hyperg.cc             ${ext_dir}/hyperg/check_value.cc
${copy} check_ibeta.cc              ${ext_dir}/ibeta/check_value.cc
${copy} check_jacobi.cc             ${ext_dir}/jacobi/check_value.cc
${copy} check_jacobi_sn.cc          ${ext_dir}/jacobi_sn/check_value.cc
${copy} check_jacobi_cn.cc          ${ext_dir}/jacobi_cn/check_value.cc
${copy} check_jacobi_dn.cc          ${ext_dir}/jacobi_dn/check_value.cc
${copy} check_lbincoef.cc           ${ext_dir}/lbincoef/check_value.cc
${copy} check_ldouble_factorial.cc  ${ext_dir}/ldouble_factorial/check_value.cc
${copy} check_legendre_q.cc         ${ext_dir}/legendre_q/check_value.cc
${copy} check_lfactorial.cc         ${ext_dir}/lfactorial/check_value.cc
${copy} check_lpochhammer_l.cc      ${ext_dir}/lpochhammer_l/check_value.cc
${copy} check_lpochhammer_u.cc      ${ext_dir}/lpochhammer_u/check_value.cc
${copy} check_owens_t.cc            ${ext_dir}/owens_t/check_value.cc
${copy} check_pochhammer_l.cc       ${ext_dir}/pochhammer_l/check_value.cc
${copy} check_pochhammer_u.cc       ${ext_dir}/pochhammer_u/check_value.cc
${copy} check_psi.cc                ${ext_dir}/psi/check_value.cc
${copy} check_radpoly.cc            ${ext_dir}/radpoly/check_value.cc
${copy} check_shi.cc                ${ext_dir}/sinhint/check_shi.cc
${copy} check_sinc.cc               ${ext_dir}/sinc/check_value.cc
${copy} check_sinc_pi.cc            ${ext_dir}/sinc_pi/check_value.cc
${copy} check_sinhint.cc            ${ext_dir}/sinhint/check_value.cc
${copy} check_sinint.cc             ${ext_dir}/sinint/check_value.cc
${copy} check_sph_bessel_i.cc       ${ext_dir}/sph_bessel_i/check_value.cc
${copy} check_sph_bessel_k.cc       ${ext_dir}/sph_bessel_k/check_value.cc
${copy} check_sph_hankel_1.cc       ${ext_dir}/sph_hankel_1/check_value.cc
${copy} check_sph_hankel_2.cc       ${ext_dir}/sph_hankel_2/check_value.cc
${copy} check_sph_harmonic.cc       ${ext_dir}/sph_harmonic/check_value.cc
${copy} check_theta_1.cc            ${ext_dir}/theta_1/check_value.cc
${copy} check_theta_2.cc            ${ext_dir}/theta_2/check_value.cc
${copy} check_theta_3.cc            ${ext_dir}/theta_3/check_value.cc
${copy} check_theta_4.cc            ${ext_dir}/theta_4/check_value.cc
${copy} check_zernike.cc            ${ext_dir}/zernike/check_value.cc
