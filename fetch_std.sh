#!/bin/bash

tool="cp -f"
gcc_dir="$HOME/gcc_specfun/libstdc++-v3/include/bits"

${tool} ${gcc_dir}/specfun_util.h    ./specfun_util.h
${tool} ${gcc_dir}/specfun.h         ./specfun.h
${tool} ${gcc_dir}/sf_bessel.tcc     ./sf_bessel.tcc
${tool} ${gcc_dir}/sf_beta.tcc       ./sf_beta.tcc
${tool} ${gcc_dir}/sf_ellint.tcc     ./sf_ellint.tcc
${tool} ${gcc_dir}/sf_expint.tcc     ./sf_expint.tcc
${tool} ${gcc_dir}/sf_gamma.tcc      ./sf_gamma.tcc
${tool} ${gcc_dir}/sf_hyperg.tcc     ./sf_hyperg.tcc
${tool} ${gcc_dir}/sf_legendre.tcc   ./sf_legendre.tcc
${tool} ${gcc_dir}/sf_mod_bessel.tcc ./sf_mod_bessel.tcc
${tool} ${gcc_dir}/sf_hermite.tcc    ./sf_hermite.tcc
${tool} ${gcc_dir}/sf_laguerre.tcc   ./sf_laguerre.tcc
${tool} ${gcc_dir}/sf_zeta.tcc       ./sf_zeta.tcc
