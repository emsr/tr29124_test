// $HOME/bin/bin/g++ -g -I. -o test_rational test_rational.cpp

// LD_LIBRARY_PATH=$HOME/bin/lib64:$LD_LIBRARY_PATH ./test_rational

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
  std::cout << rpoly << '\n';
}
