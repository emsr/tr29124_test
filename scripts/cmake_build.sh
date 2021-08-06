#! /bin/bash

src_dir=$(pwd)
echo "Source directory: $src_dir"

mkdir -p $HOME/builds/cxx_math/release
cd $HOME/builds/cxx_math/release
cmake -DCMAKE_BUILD_TYPE=release $src_dir
make -j$(nproc)

mkdir -p $HOME/builds/cxx_math/debug
cd $HOME/builds/cxx_math/debug
cmake -DCMAKE_BUILD_TYPE=debug $src_dir
make -j$(nproc)
