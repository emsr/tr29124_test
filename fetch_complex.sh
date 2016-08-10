#!  /bin/bash

tool="cp -f"

suffix="_tr29124"
if [ $# -ge 1 ]; then
  suffix="$1"
fi

${tool} ../gcc${suffix}/libstdc++-v3/testsuite/ext/special_functions/ellint_rc/complex.cc complex_ellint_rc.cc 
${tool} ../gcc${suffix}/libstdc++-v3/testsuite/ext/special_functions/ellint_rd/complex.cc complex_ellint_rd.cc 
${tool} ../gcc${suffix}/libstdc++-v3/testsuite/ext/special_functions/ellint_rf/complex.cc complex_ellint_rf.cc 
${tool} ../gcc${suffix}/libstdc++-v3/testsuite/ext/special_functions/ellint_rg/complex.cc complex_ellint_rg.cc 
${tool} ../gcc${suffix}/libstdc++-v3/testsuite/ext/special_functions/ellint_rj/complex.cc complex_ellint_rj.cc 
