
/*
$HOME/bin/bin/g++ -std=gnu++17 -g -Wall -Wextra -I. -o test_dawson test_dawson.cpp -L$HOME/bin/lib64 -lquadmath
./test_dawson > test_dawson.txt
*/

#include <limits>
#include <iostream>
#include <iomanip>
#include <cmath>
#include <bits/float128_io.h>

int
main()
{
  using _Tp = __float128;
  std::cout.precision(__gnu_cxx::__digits10<_Tp>());
  constexpr auto _S_1_sqrtpi{0.5641895835477562869480794515607726L};
  constexpr auto _S_H{0.2L};
  /// @todo this needs some compile-time construction!
  constexpr auto _S_n_max = 100;
  static _Tp _S_c[_S_n_max];
  static auto __init = false;
  _S_c[0] = _Tp{1};
  for (unsigned int __i = 0; __i < _S_n_max; ++__i)
   {
     _Tp __y = _Tp(2 * __i + 1) * _S_H;
     _S_c[__i] = std::exp(-__y * __y);
   }
  std::cout << std::scientific;
  for (auto __CC : _S_c)
    std::cout << '\t' << __CC << "L,\n";
}
