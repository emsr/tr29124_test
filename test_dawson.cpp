
// g++ -std=c++14 -o test_dawson test_dawson.cpp

// ./test_dawson

#include <iostream>
#include <limits>
#include <cmath>

int
main()
{
  using _Tp = __float128;
  std::cout.precision(std::numeric_limits<_Tp>::digits10);
  constexpr auto _S_1_sqrtpi{0.5641895835477562869480794515607726L};
  constexpr auto _S_H(0.2L);
  /// @todo this needs some compile-time construction!
  constexpr auto __n_max = 100;
  static _Tp __c[__n_max + 1];
  static auto __init = false;
  for (unsigned int __i = 1; __i <= __n_max; ++__i)
   {
     _Tp __y = _Tp(2 * __i - 1) * _S_H;
     __c[__i] = std::exp(-__y * __y);
   }
  for (auto __CC : __c)
    std::cout << ' ' << __CC << ',' << '\n';
}
