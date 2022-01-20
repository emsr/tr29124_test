#! /bin/bash

dir=build/release
mkdir -p $dir
cd $dir
cmake ../.. -DCMAKE_BUILD_TYPE=Release
make
make test
