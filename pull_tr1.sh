#!/bin/bash

tool="cp -f"

suffix="_tr29124"
if [ $# -ge 1 ]; then
  suffix="$1"
fi

src_dir="$HOME/gcc${suffix}/libstdc++-v3/include"
if [ -d "$1" ]; then 
  src_dir="$1"
fi
src_tr1_dir="${src_dir}/tr1"

dst_dir="."
if [ -d "$2" ]; then 
  dst_dir="$2"
fi
dst_tr1_dir="${dst_dir}/tr1"

${tool} "${src_tr1_dir}/cmath"                    "${dst_tr1_dir}"
${tool} "${src_tr1_dir}/bessel_function.tcc"      "${dst_tr1_dir}"
${tool} "${src_tr1_dir}/beta_function.tcc"        "${dst_tr1_dir}"
${tool} "${src_tr1_dir}/ell_integral.tcc"         "${dst_tr1_dir}"
${tool} "${src_tr1_dir}/exp_integral.tcc"         "${dst_tr1_dir}"
${tool} "${src_tr1_dir}/gamma.tcc"                "${dst_tr1_dir}"
${tool} "${src_tr1_dir}/hypergeometric.tcc"       "${dst_tr1_dir}"
${tool} "${src_tr1_dir}/legendre_function.tcc"    "${dst_tr1_dir}"
${tool} "${src_tr1_dir}/modified_bessel_func.tcc" "${dst_tr1_dir}"
${tool} "${src_tr1_dir}/poly_hermite.tcc"         "${dst_tr1_dir}"
${tool} "${src_tr1_dir}/poly_laguerre.tcc"        "${dst_tr1_dir}"
${tool} "${src_tr1_dir}/riemann_zeta.tcc"         "${dst_tr1_dir}"
${tool} "${src_tr1_dir}/special_function_util.h"  "${dst_tr1_dir}"
