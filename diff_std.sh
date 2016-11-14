#!/bin/bash

tool="kdiff3"

suffix="_tr29124"
if [ $# -ge 1 ]; then
  suffix="$1"
fi

inc_dir="$HOME/gcc${suffix}/libstdc++-v3/include"
gcc_dir="$HOME/gcc${suffix}/libstdc++-v3/include/bits"
ext_dir="$HOME/gcc${suffix}/libstdc++-v3/include/ext"

${tool} ext/   ${ext_dir}/
${tool} bits/  ${gcc_dir}/

${tool} cmath                      ${inc_dir}/c_global/cmath
${tool} limits                     ${inc_dir}/std/limits
${tool} complex                    ${inc_dir}/std/complex
${tool} ext/math_const.h           ${ext_dir}/math_const.h
${tool} ext/math_util.h            ${ext_dir}/math_util.h
${tool} ext/polynomial.h           ${ext_dir}/polynomial.h
${tool} ext/polynomial.tcc         ${ext_dir}/polynomial.tcc
${tool} ext/cmath                  ${gcc_dir}/cmath
${tool} bits/specfun_util.h        ${gcc_dir}/specfun_util.h
${tool} bits/specfun_state.h       ${gcc_dir}/specfun_state.h
${tool} bits/specfun.h             ${gcc_dir}/specfun.h
${tool} bits/complex_util.h        ${gcc_dir}/complex_util.h
${tool} bits/complex_util.tcc      ${gcc_dir}/complex_util.tcc
${tool} bits/float128_io.h         ${gcc_dir}/float128_io.h
${tool} bits/float128_io.tcc       ${gcc_dir}/float128_io.tcc
${tool} bits/float128_limits.h     ${gcc_dir}/float128_limits.h
${tool} bits/float128_math.h       ${gcc_dir}/float128_math.h
${tool} bits/numeric_limits_float128.h     ${gcc_dir}/numeric_limits_float128.h
${tool} bits/summation.h           ${gcc_dir}/summation.h
${tool} bits/summation.tcc         ${gcc_dir}/summation.tcc
${tool} bits/sf_airy.tcc           ${gcc_dir}/sf_airy.tcc
${tool} bits/sf_bessel.tcc         ${gcc_dir}/sf_bessel.tcc
${tool} bits/sf_beta.tcc           ${gcc_dir}/sf_beta.tcc
${tool} bits/sf_cardinal.tcc       ${gcc_dir}/sf_cardinal.tcc
${tool} bits/sf_chebyshev.tcc      ${gcc_dir}/sf_chebyshev.tcc
${tool} bits/sf_dawson.tcc         ${gcc_dir}/sf_dawson.tcc
${tool} bits/sf_distributions.tcc  ${gcc_dir}/sf_distributions.tcc
${tool} bits/sf_ellint.tcc         ${gcc_dir}/sf_ellint.tcc
${tool} bits/sf_expint.tcc         ${gcc_dir}/sf_expint.tcc
${tool} bits/sf_fresnel.tcc        ${gcc_dir}/sf_fresnel.tcc
${tool} bits/sf_gamma.tcc          ${gcc_dir}/sf_gamma.tcc
${tool} bits/sf_gegenbauer.tcc     ${gcc_dir}/sf_gegenbauer.tcc
${tool} bits/sf_hankel.tcc         ${gcc_dir}/sf_hankel.tcc
${tool} bits/sf_hermite.tcc        ${gcc_dir}/sf_hermite.tcc
${tool} bits/sf_hydrogen.tcc       ${gcc_dir}/sf_hydrogen.tcc
${tool} bits/sf_hyperg.tcc         ${gcc_dir}/sf_hyperg.tcc
${tool} bits/sf_hypint.tcc         ${gcc_dir}/sf_hypint.tcc
${tool} bits/sf_jacobi.tcc         ${gcc_dir}/sf_jacobi.tcc
${tool} bits/sf_laguerre.tcc       ${gcc_dir}/sf_laguerre.tcc
${tool} bits/sf_legendre.tcc       ${gcc_dir}/sf_legendre.tcc
${tool} bits/sf_mod_bessel.tcc     ${gcc_dir}/sf_mod_bessel.tcc
${tool} bits/sf_owens_t.tcc        ${gcc_dir}/sf_owens_t.tcc
${tool} bits/sf_polylog.tcc        ${gcc_dir}/sf_polylog.tcc 
${tool} bits/sf_theta.tcc          ${gcc_dir}/sf_theta.tcc
${tool} bits/sf_trig.tcc           ${gcc_dir}/sf_trig.tcc
${tool} bits/sf_trigint.tcc        ${gcc_dir}/sf_trigint.tcc
${tool} bits/sf_zeta.tcc           ${gcc_dir}/sf_zeta.tcc
