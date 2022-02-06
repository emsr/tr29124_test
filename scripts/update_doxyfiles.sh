#! /bin/bash

# Run around regenerating Doxyfiles. Both subdirectories and proper submodules are done here.
# Subsubmodules are ignored.
# You'll need to do the commits by hand depending on subdirectory or submodule situation.

for file in cxx_differentiation cxx_math_constants cxx_chebyshev cxx_integration cxx_continued_fractions cxx_polynomial 
do
  cd $file
  doxygen -u
  cd ..
done
