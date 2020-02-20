#! /bin/bash

mkdir -p $HOME/builds/tr29124_test/release
cd $HOME/builds/tr29124_test/release
cmake -DCMAKE_BUILD_TYPE=release $HOME/tr29124_test
make -j$(nproc)

mkdir -p $HOME/builds/tr29124_test/debug
cd $HOME/builds/tr29124_test/debug
cmake -DCMAKE_BUILD_TYPE=debug $HOME/tr29124_test
make -j$(nproc)
