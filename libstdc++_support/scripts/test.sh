#!/bin/bash

suffix="_tr29124"
if [ $# -ge 1 ]; then
  suffix="$1"
fi

#  This speeds up things like crazy!!!
cd ~/obj${suffix}/x86_64-pc-linux-gnu/libstdc++-v3/testsuite
make -k check RUNTESTFLAGS="conformance.exp=check_value.cc"

#  Directories alsowork!!!
make -k check RUNTESTFLAGS="conformance.exp=special_functions/*/check_value.cc"
make -k check RUNTESTFLAGS="conformance.exp=ext/special_functions/*/check_value.cc"
make -k check RUNTESTFLAGS="conformance.exp=tr1/5_numerical_facilities/special_functions/*/check_value.cc"
make -k check RUNTESTFLAGS="conformance.exp=tr1/5_numerical_facilities/special_functions/*/*.cc"

