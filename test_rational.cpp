/*
$HOME/bin/bin/g++ -g -Wall -Wextra -Wno-psabi -I. -o test_rational test_rational.cpp -lquadmath
./test_rational > test_rational.txt
*/

#include "rational.h"
#include <ext/polynomial.h>

int
main()
{
  using Rat = __gnu_cxx::_Rational<int>;

  std::cout << Rat(1, 2) + Rat(1, 3) << '\n';
  std::cout << Rat(1, 2) - Rat(1, 3) << '\n';
  std::cout << Rat(1, 2) * Rat(1, 3) << '\n';
  std::cout << Rat(1, 2) / Rat(1, 3) << '\n';

  __gnu_cxx::_Polynomial<Rat> rpoly{{1,2}, {3, 4}, {5, 6}, {7, 8}};
  std::cout << "P = " << rpoly << '\n';
  for (int i = 0; i <= 10; ++i)
    std::cout << "P(" << i << ") = " << rpoly(i*1.0) << '\n';
}
