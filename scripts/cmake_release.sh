#! /bin/bash

mkdir -p $HOME/builds/cxx_math/release
cd $HOME/builds/cxx_math/release
cmake -DCMAKE_BUILD_TYPE=release $HOME/cxx_math
make -j$(nproc)
