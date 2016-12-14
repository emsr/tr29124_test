/*
$HOME/bin_tr29124/bin/g++ -std=gnu++17 -g -Wall -Wextra -Wno-psabi -I. -o test_maxint test_maxint.cpp -lquadmath
./test_maxint > test_maxint.txt

$HOME/bin/bin/g++ -std=gnu++17 -g -Wall -Wextra -Wno-psabi -I. -o test_maxint test_maxint.cpp -lquadmath
./test_maxint > test_maxint.txt
*/

#include <mpreal.h>

#include <limits>
#include <iostream>
#include <ext/cmath>
#include <bits/numeric_limits_mpreal.h>

template<typename _Tp>
  void
  test_maxint(_Tp proto = _Tp{})
  {
    std::cout.precision(__gnu_cxx::__max_digits10(proto));
    auto width = std::cout.precision() + 8;
    std::cout << std::showpoint << std::scientific;

    // Try 1/epsilon.
    auto maxint = _Tp{1} / std::numeric_limits<_Tp>::epsilon();
    std::cout << "\n\nTrying maxint = " << std::setw(width) << maxint << '\n';
    if (maxint + 1 == maxint)
      std::cout << "\nmaxint FAIL\n";
    for (int i = 1; i < 100; ++i, maxint += _Tp{1})
      if (maxint + 1 == maxint)
	{
	  std::cout << "\nmaxint FAIL at " << maxint << "\n";
	  break;
	}

    // Try ldexp(1, std::numeric_limis<_Tp>::digits);
    auto maxint2 = std::ldexp(_Tp{1}, std::numeric_limits<_Tp>::digits);
    std::cout << "\n\nTrying maxint2 = " << std::setw(width) << maxint2 << '\n';
    if (maxint2 + 1 == maxint2)
      std::cout << "\nmaxint2 FAIL\n";
    for (int i = 1; i < 100; ++i, maxint2 += _Tp{1})
      if (maxint2 + 1 == maxint2)
	{
	  std::cout << "\nmaxint2 FAIL at " << maxint2 << "\n";
	  break;
	}
  }

int
main()
{
  test_maxint<float>();
  test_maxint<double>();
  test_maxint<long double>();
  test_maxint<__float128>();
  //test_maxint<mpfr::mpreal>((1,  256));
}
