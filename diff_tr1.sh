#!/bin/bash

tool="kdiff3"

suffix="_tr29124"
if [ $# -ge 1 ]; then
  suffix="$1"
fi

gcc_dir="$HOME/gcc${suffix}/libstdc++-v3/include/tr1"

${tool} ./special_function_util.h  ${gcc_dir}/special_function_util.h
${tool} ./special_function.h	   ${gcc_dir}/special_function.h
${tool} ./bessel_function.tcc	   ${gcc_dir}/bessel_function.tcc
${tool} ./beta_function.tcc	   ${gcc_dir}/beta_function.tcc
${tool} ./ell_integral.tcc	   ${gcc_dir}/ell_integral.tcc
${tool} ./exp_integral.tcc	   ${gcc_dir}/exp_integral.tcc
${tool} ./gamma.tcc		   ${gcc_dir}/gamma.tcc
${tool} ./hypergeometric.tcc	   ${gcc_dir}/hypergeometric.tcc
${tool} ./legendre_function.tcc    ${gcc_dir}/legendre_function.tcc
${tool} ./modified_bessel_func.tcc ${gcc_dir}/modified_bessel_func.tcc
${tool} ./poly_hermite.tcc	   ${gcc_dir}/poly_hermite.tcc
${tool} ./poly_laguerre.tcc	   ${gcc_dir}/poly_laguerre.tcc
${tool} ./riemann_zeta.tcc	   ${gcc_dir}/riemann_zeta.tcc
