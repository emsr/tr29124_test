#!/bin/bash

tool="kdiff3"

suffix="_tr29124"
if [ $# -ge 1 ]; then
  suffix="$1"
fi

tr29124_dir="$HOME/gcc${suffix}/libstdc++-v3/testsuite/special_functions"
ext_dir="$HOME/gcc${suffix}/libstdc++-v3/testsuite/ext/special_functions"

#tr29124_dir="$HOME/gcc${suffix}/libstdc++-v3/testsuite/special_functions"
#ext_dir="$HOME/gcc${suffix}/libstdc++-v3/testsuite/ext/special_functions"

tr1_dir="$HOME/gcc${suffix}/libstdc++-v3/testsuite/tr1/5_numerical_facilities/special_functions"

${tool} ${tr29124_dir}/01_assoc_laguerre/	${tr1_dir}/01_assoc_laguerre/
${tool} ${tr29124_dir}/02_assoc_legendre/	${tr1_dir}/02_assoc_legendre/
${tool} ${tr29124_dir}/03_beta/			${tr1_dir}/03_beta/
${tool} ${tr29124_dir}/04_comp_ellint_1/	${tr1_dir}/04_comp_ellint_1/
${tool} ${tr29124_dir}/05_comp_ellint_2/	${tr1_dir}/05_comp_ellint_2/
${tool} ${tr29124_dir}/06_comp_ellint_3/	${tr1_dir}/06_comp_ellint_3/
${tool} ${ext_dir}/conf_hyperg/			${tr1_dir}/07_conf_hyperg/
${tool} ${tr29124_dir}/07_cyl_bessel_i/		${tr1_dir}/08_cyl_bessel_i/
${tool} ${tr29124_dir}/08_cyl_bessel_j/		${tr1_dir}/09_cyl_bessel_j/
${tool} ${tr29124_dir}/09_cyl_bessel_k/		${tr1_dir}/10_cyl_bessel_k/
${tool} ${tr29124_dir}/10_cyl_neumann/		${tr1_dir}/11_cyl_neumann/
${tool} ${tr29124_dir}/11_ellint_1/		${tr1_dir}/12_ellint_1/
${tool} ${tr29124_dir}/12_ellint_2/		${tr1_dir}/13_ellint_2/
${tool} ${tr29124_dir}/13_ellint_3/		${tr1_dir}/14_ellint_3/
${tool} ${tr29124_dir}/14_expint/		${tr1_dir}/15_expint/
${tool} ${tr29124_dir}/15_hermite/		${tr1_dir}/16_hermite/
${tool} ${ext_dir}/hyperg/			${tr1_dir}/17_hyperg/
${tool} ${tr29124_dir}/16_laguerre/		${tr1_dir}/18_laguerre/
${tool} ${tr29124_dir}/17_legendre/		${tr1_dir}/19_legendre/
${tool} ${tr29124_dir}/18_riemann_zeta/		${tr1_dir}/20_riemann_zeta/
${tool} ${tr29124_dir}/19_sph_bessel/		${tr1_dir}/21_sph_bessel/
${tool} ${tr29124_dir}/20_sph_legendre/		${tr1_dir}/22_sph_legendre/
${tool} ${tr29124_dir}/21_sph_neumann/		${tr1_dir}/23_sph_neumann/     
