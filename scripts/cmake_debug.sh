#! /bin/bash

dir=build/debug
mkdir -p $dir
cd $dir
cmake ../.. -DCMAKE_BUILD_TYPE=Debug
make
make test
