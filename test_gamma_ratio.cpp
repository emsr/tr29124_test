/*
$HOME/bin_tr29124/bin/g++ -std=gnu++17 -g -Wall -Wextra -I. -o test_gamma_ratio test_gamma_ratio.cpp wrap_boost.cpp -lquadmath
./test_gamma_ratio > test_gamma_ratio.txt

*/

#include <ext/cmath>
#include <bits/float128_io.h>
#include <limits>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <vector>
#include <string>
#include <complex>

#include "wrap_boost.h"

/**
 * Return the Stirling numbers of the first kind.
 *
 */
template<typename _Tp>
  _Tp
  stirling_number_1(int __n, int __k)
  {
    return _Tp{0};
  }

/**
 * Return the Stirling numbers of the first kind.
 *
 */
template<typename _Tp>
  _Tp
  stirling_number_1(int __n, _Tp __kappa)
  {
    return _Tp{0};
  }

/**
 * Return the ratio
 * @f[
 *   \frac{\Gamma(a+\alpha)}{\Gamma(a)}
 * @f]
 * for large @f$ a @f$.
 */
template<typename _Tp>
  _Tp
  __gamma_ratio_asymp()
  {
  }

template<typename _Tp>
  void
  test_gamma_ratio(Tp proto = _Tp{})
  {
    using _Val = _Tp;
    using _Real = std::__detail::__num_traits_t<_Val>;

    std::cout.precision(std::numeric_limits<_Real>::digits10);
    std::cout << std::showpoint << std::scientific;
    auto width = 8 + std::cout.precision();
  }

int
main()
{
  test_gamma_ratio(1.0f);
  
  test_gamma_ratio(1.0);

  test_gamma_ratio(1.0l);

  test_gamma_ratio(1.0q);
}
