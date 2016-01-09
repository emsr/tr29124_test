#!/bin/bash

tool="cp -f"
gcc_dir="$HOME/gcc_specfun/libstdc++-v3/include/bits"
ext_dir="$HOME/gcc_specfun/libstdc++-v3/include/ext"

${tool} ${gcc_dir}/specfun.h         ./specfun.h
${tool} ${gcc_dir}/specfun_util.h    ./specfun_util.h
${tool} ${gcc_dir}/complex_util.h    ./complex_util.h
${tool} ${gcc_dir}/complex_util.tcc  ./complex_util.tcc
${tool} ${gcc_dir}/float128.h        ./float128.h
${tool} ${gcc_dir}/float128.tcc      ./float128.tcc
${tool} ${ext_dir}/math_const.h      ./math_const.h
${tool} ${gcc_dir}/sf_bessel.tcc     ./sf_bessel.tcc
${tool} ${gcc_dir}/sf_beta.tcc       ./sf_beta.tcc
${tool} ${gcc_dir}/sf_chebyshev.tcc  ./sf_chebyshev.tcc
${tool} ${gcc_dir}/sf_dawson.tcc     ./sf_dawson.tcc
${tool} ${gcc_dir}/sf_ellint.tcc     ./sf_ellint.tcc
${tool} ${gcc_dir}/sf_expint.tcc     ./sf_expint.tcc
${tool} ${gcc_dir}/sf_fresnel.tcc    ./sf_fresnel.tcc
${tool} ${gcc_dir}/sf_gamma.tcc      ./sf_gamma.tcc
${tool} ${gcc_dir}/sf_gegenbauer.tcc ./sf_gegenbauer.tcc
${tool} ${gcc_dir}/sf_hermite.tcc    ./sf_hermite.tcc
${tool} ${gcc_dir}/sf_hydrogen.tcc   ./sf_hydrogen.tcc
${tool} ${gcc_dir}/sf_hyperg.tcc     ./sf_hyperg.tcc
${tool} ${gcc_dir}/sf_hypint.tcc     ./sf_hypint.tcc
${tool} ${gcc_dir}/sf_jacobi.tcc     ./sf_jacobi.tcc
${tool} ${gcc_dir}/sf_laguerre.tcc   ./sf_laguerre.tcc
${tool} ${gcc_dir}/sf_legendre.tcc   ./sf_legendre.tcc
${tool} ${gcc_dir}/sf_mod_bessel.tcc ./sf_mod_bessel.tcc
${tool} ${gcc_dir}/sf_theta.tcc      ./sf_theta.tcc
${tool} ${gcc_dir}/sf_trigint.tcc    ./sf_trigint.tcc
${tool} ${gcc_dir}/sf_zeta.tcc       ./sf_zeta.tcc
${tool} ${gcc_dir}/sf_airy.tcc       ./sf_airy.tcc
${tool} ${gcc_dir}/sf_hankel.tcc     ./sf_hankel.tcc
