#!/bin/bash

tool="cp -f"
gcc_dir="$HOME/gcc_specfun/libstdc++-v3/include/bits"
ext_dir="$HOME/gcc_specfun/libstdc++-v3/include/ext"
bin_dir="$HOME/bin_specfun/include/c++/6.0.0/bits"

${tool} ${gcc_dir}/specfun.h         ${bin_dir}/specfun.h
${tool} ${gcc_dir}/specfun_util.h    ${bin_dir}/specfun_util.h
${tool} ${gcc_dir}/complex_util.h    ${bin_dir}/complex_util.h
${tool} ${gcc_dir}/complex_util.tcc  ${bin_dir}/complex_util.tcc
${tool} ${gcc_dir}/float128.h        ${bin_dir}/float128.h
${tool} ${gcc_dir}/float128.tcc      ${bin_dir}/float128.tcc
${tool} ${gcc_dir}/numeric_limits.h  ${bin_dir}/numeric_limits.h
${tool} ${ext_dir}/math_const.h      ${bin_dir}/math_const.h
${tool} ${gcc_dir}/sf_bessel.tcc     ${bin_dir}/sf_bessel.tcc
${tool} ${gcc_dir}/sf_beta.tcc       ${bin_dir}/sf_beta.tcc
${tool} ${gcc_dir}/sf_cardinal.tcc   ${bin_dir}/sf_cardinal.tcc
${tool} ${gcc_dir}/sf_chebyshev.tcc  ${bin_dir}/sf_chebyshev.tcc
${tool} ${gcc_dir}/sf_dawson.tcc     ${bin_dir}/sf_dawson.tcc
${tool} ${gcc_dir}/sf_ellint.tcc     ${bin_dir}/sf_ellint.tcc
${tool} ${gcc_dir}/sf_expint.tcc     ${bin_dir}/sf_expint.tcc
${tool} ${gcc_dir}/sf_fresnel.tcc    ${bin_dir}/sf_fresnel.tcc
${tool} ${gcc_dir}/sf_gamma.tcc      ${bin_dir}/sf_gamma.tcc
${tool} ${gcc_dir}/sf_gegenbauer.tcc ${bin_dir}/sf_gegenbauer.tcc
${tool} ${gcc_dir}/sf_hermite.tcc    ${bin_dir}/sf_hermite.tcc
${tool} ${gcc_dir}/sf_hydrogen.tcc   ${bin_dir}/sf_hydrogen.tcc
${tool} ${gcc_dir}/sf_hyperg.tcc     ${bin_dir}/sf_hyperg.tcc
${tool} ${gcc_dir}/sf_hypint.tcc     ${bin_dir}/sf_hypint.tcc
${tool} ${gcc_dir}/sf_jacobi.tcc     ${bin_dir}/sf_jacobi.tcc
${tool} ${gcc_dir}/sf_laguerre.tcc   ${bin_dir}/sf_laguerre.tcc
${tool} ${gcc_dir}/sf_legendre.tcc   ${bin_dir}/sf_legendre.tcc
${tool} ${gcc_dir}/sf_mod_bessel.tcc ${bin_dir}/sf_mod_bessel.tcc
${tool} ${gcc_dir}/sf_theta.tcc      ${bin_dir}/sf_theta.tcc
${tool} ${gcc_dir}/sf_trigint.tcc    ${bin_dir}/sf_trigint.tcc
${tool} ${gcc_dir}/sf_zeta.tcc       ${bin_dir}/sf_zeta.tcc
${tool} ${gcc_dir}/sf_owens_t.tcc    ${bin_dir}/sf_owens_t.tcc
${tool} ${gcc_dir}/sf_polylog.tcc    ${bin_dir}/sf_polylog.tcc
${tool} ${gcc_dir}/sf_airy.tcc       ${bin_dir}/sf_airy.tcc
${tool} ${gcc_dir}/sf_hankel.tcc     ${bin_dir}/sf_hankel.tcc
