/*
$HOME/bin/bin/g++ -std=c++2a -o test_ulp test_ulp.cpp
*/

#include <cmath>
#include <limits>
#include <iostream>

template<typename _Tp>
  _Tp
  ulp(_Tp __x)
  {
    const int digs = std::numeric_limits<_Tp>::digits;
    const int min_exp = std::numeric_limits<_Tp>::min_exponent;
    int exp;
    /*_Tp __frac = */std::frexp(std::abs(__x), &exp);
    return ldexp(_Tp{1}, std::max(exp, min_exp) - digs + 1);
  }

int
main()
{
  auto a = ulp(-123.457);
  std::cout << "ulp(-123.457) = " << a << '\n';
}
