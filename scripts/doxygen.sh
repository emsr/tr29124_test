#! /bin/bash

# Run around regenerating doxygen. Both subdirectories and proper submodules are done here.
# Subsubmodules are ignored.
# You'll need to do the commits by hand depending on subdirectory or submodule situation.

for file in cxx_differentiation cxx_math_constants cxx_chebyshev cxx_integration cxx_continued_fractions cxx_polynomial 
do
  cd $file
  rm -rf docs/html
  rm -rf docs/latex
  doxygen
  git add docs
  cd ..
done

# Top-level
rm -rf docs/html
rm -rf docs/latex
doxygen
git add docs
