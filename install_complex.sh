#!/bin/bash

tool="cp -f"

suffix="_tr29124"
if [ $# -ge 1 ]; then
  suffix="$1"
fi

${tool} check/complex_ellint_rc.cc ../gcc${suffix}/libstdc++-v3/testsuite/ext/special_functions/ellint_rc/complex.cc
${tool} check/complex_ellint_rd.cc ../gcc${suffix}/libstdc++-v3/testsuite/ext/special_functions/ellint_rd/complex.cc
${tool} check/complex_ellint_rf.cc ../gcc${suffix}/libstdc++-v3/testsuite/ext/special_functions/ellint_rf/complex.cc
${tool} check/complex_ellint_rg.cc ../gcc${suffix}/libstdc++-v3/testsuite/ext/special_functions/ellint_rg/complex.cc
${tool} check/complex_ellint_rj.cc ../gcc${suffix}/libstdc++-v3/testsuite/ext/special_functions/ellint_rj/complex.cc
${tool} check/complex_airy_ai.cc ../gcc${suffix}/libstdc++-v3/testsuite/ext/special_functions/airy_ai/complex.cc
${tool} check/complex_airy_bi.cc ../gcc${suffix}/libstdc++-v3/testsuite/ext/special_functions/airy_bi/complex.cc
