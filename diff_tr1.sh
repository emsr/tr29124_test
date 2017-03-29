#!/bin/bash

tool="kdiff3"

suffix="_tr29124"
if [ $# -ge 1 ]; then
  suffix="$1"
fi

gcc_dir="$HOME/gcc${suffix}/libstdc++-v3/include"
gcc_tr1_dir="${gcc_dir}/tr1"
gcc_bits_dir="${gcc_dir}/bits"


${tool} ./tr1/special_function_util.h  ${gcc_tr1_dir}/special_function_util.h
${tool} ./bits/specfun.h               ${gcc_bits_dir}/specfun.h
${tool} ./tr1/bessel_function.tcc      ${gcc_tr1_dir}/bessel_function.tcc
${tool} ./tr1/beta_function.tcc	       ${gcc_tr1_dir}/beta_function.tcc
${tool} ./tr1/ell_integral.tcc	       ${gcc_tr1_dir}/ell_integral.tcc
${tool} ./tr1/exp_integral.tcc	       ${gcc_tr1_dir}/exp_integral.tcc
${tool} ./tr1/gamma.tcc		       ${gcc_tr1_dir}/gamma.tcc
${tool} ./tr1/hypergeometric.tcc       ${gcc_tr1_dir}/hypergeometric.tcc
${tool} ./tr1/legendre_function.tcc    ${gcc_tr1_dir}/legendre_function.tcc
${tool} ./tr1/modified_bessel_func.tcc ${gcc_tr1_dir}/modified_bessel_func.tcc
${tool} ./tr1/poly_hermite.tcc	       ${gcc_tr1_dir}/poly_hermite.tcc
${tool} ./tr1/poly_laguerre.tcc	       ${gcc_tr1_dir}/poly_laguerre.tcc
${tool} ./tr1/riemann_zeta.tcc	       ${gcc_tr1_dir}/riemann_zeta.tcc
