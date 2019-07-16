#! /bin/bash

rm -rf $(find . -name CMakeCache.txt)
rm -rf $(find . -name CMakeFiles)
# Can't do this for the whole tree:
#rm -rf $(find . -name Makefile)
rm -rf $(find . -name cmake_install.cmake)

rm -rf $(find . -name CTestTestfile.cmake)
rm -rf $(find . -name Testing)

rm -rf $(find . -name cmake_install.cmake)

rm -rf $(find . -name CMakeDoxyfile.in)
rm -rf $(find . -name MakeDoxygenDefaults.cmake)

