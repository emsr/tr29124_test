/**
 *
 */

#include <limits>
#include <iostream>
#include <iomanip>
#include <cmath>

#include <emsr/float128_io.h>

int
main()
{
  //using _Tp = __float128;
  using _Tp = long double;
  std::cout.precision(emsr::digits10<_Tp>());
  //constexpr auto s_1_sqrtpi{0.5641895835477562869480794515607726L};
  constexpr auto s_H{0.2L};
  /// @todo this needs some compile-time construction!
  constexpr auto s_n_max = 100;
  static _Tp s_c[s_n_max];
  //static auto init = false;
  s_c[0] = _Tp{1};
  for (unsigned int i = 0; i < s_n_max; ++i)
   {
     _Tp y = _Tp(2 * i + 1) * s_H;
     s_c[i] = std::exp(-y * y);
   }
  std::cout << std::scientific;
  for (auto CC : s_c)
    std::cout << '\t' << CC << "L,\n";
}
