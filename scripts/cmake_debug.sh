#! /bin/bash

mkdir -p $HOME/builds/cxx_math/debug
cd $HOME/builds/cxx_math/debug
cmake -DCMAKE_BUILD_TYPE=debug $HOME/cxx_math
make -j$(nproc)
