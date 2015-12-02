#!/bin/bash

tool="kdiff3"
gcc_dir="$HOME/gcc_specfun/libstdc++-v3/include/bits"
ext_dir="$HOME/gcc_specfun/libstdc++-v3/include/ext"

${tool} ./specfun_util.h    ${gcc_dir}/specfun_util.h
${tool} ./specfun.h         ${gcc_dir}/specfun.h
${tool} ./math_const.h      ${ext_dir}/math_const.h
${tool} ./sf_bessel.tcc     ${gcc_dir}/sf_bessel.tcc
${tool} ./sf_beta.tcc       ${gcc_dir}/sf_beta.tcc
${tool} ./sf_ellint.tcc     ${gcc_dir}/sf_ellint.tcc
${tool} ./sf_expint.tcc     ${gcc_dir}/sf_expint.tcc
${tool} ./sf_gamma.tcc      ${gcc_dir}/sf_gamma.tcc
${tool} ./sf_hyperg.tcc     ${gcc_dir}/sf_hyperg.tcc
${tool} ./sf_legendre.tcc   ${gcc_dir}/sf_legendre.tcc
${tool} ./sf_mod_bessel.tcc ${gcc_dir}/sf_mod_bessel.tcc
${tool} ./sf_hermite.tcc    ${gcc_dir}/sf_hermite.tcc
${tool} ./sf_laguerre.tcc   ${gcc_dir}/sf_laguerre.tcc
${tool} ./sf_zeta.tcc       ${gcc_dir}/sf_zeta.tcc
