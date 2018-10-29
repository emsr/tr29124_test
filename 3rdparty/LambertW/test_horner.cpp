/*
$HOME/bin/bin/g++ -std=c++17 -o test_horner test_horner.cpp

clang++ -std=c++1z -stdlib=libc++ -o test_horner test_horner.cpp
*/

#include "horner.h"

#include <iostream>
#include <iomanip>

int
main()
{
  std::cout << horner_big_end(3.0, 1.0, 2.0, 1.0) << '\n';
  std::cout << horner_big_end(3.0, 1.0, 2.0, 3.0) << '\n';
  std::cout << horner(3.0, 3.0, 2.0, 1.0) << '\n';
}
