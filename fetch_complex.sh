#!  /bin/bash

tool="cp -f"

suffix="_tr29124"
if [ $# -ge 1 ]; then
  suffix="$1"
fi

${tool} ../gcc${suffix}/libstdc++-v3/testsuite/ext/special_functions/ellint_rc/complex.cc check/complex_ellint_rc.cc
${tool} ../gcc${suffix}/libstdc++-v3/testsuite/ext/special_functions/ellint_rd/complex.cc check/complex_ellint_rd.cc
${tool} ../gcc${suffix}/libstdc++-v3/testsuite/ext/special_functions/ellint_rf/complex.cc check/complex_ellint_rf.cc
${tool} ../gcc${suffix}/libstdc++-v3/testsuite/ext/special_functions/ellint_rg/complex.cc check/complex_ellint_rg.cc
${tool} ../gcc${suffix}/libstdc++-v3/testsuite/ext/special_functions/ellint_rj/complex.cc check/complex_ellint_rj.cc
${tool} ../gcc${suffix}/libstdc++-v3/testsuite/ext/special_functions/airy_ai/complex.cc check/complex_airy_ai.cc
${tool} ../gcc${suffix}/libstdc++-v3/testsuite/ext/special_functions/airy_bi/complex.cc check/complex_airy_bi.cc
