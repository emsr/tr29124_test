// $HOME/bin_specfun/bin/g++ -std=gnu++14 -g -Wall -Wextra -o test_continued_fraction test_continued_fraction.cpp -L$HOME/bin/lib64 -lquadmath

// LD_LIBRARY_PATH=$HOME/bin_specfun/lib64:$LD_LIBRARY_PATH ./test_continued_fraction > test_continued_fraction.txt

// g++ -std=gnu++14 -Wall -Wextra -DNO_LOGBQ -I. -o test_continued_fraction test_continued_fraction.cpp -lquadmath

// ./test_continued_fraction > test_continued_fraction.txt

#include <limits>
#include <iostream>
#include <iomanip>
#include <vector>
#include <complex>
#include <ext/math_const.h>
#include <bits/float128.h>
#include <bits/numeric_limits.h>
#include <bits/specfun_util.h>
#include <bits/complex_util.h>
#include "polynomial.h"
//#include "rational.h"

/**
 * Convert a series to a continued fraction.
 * An example using the confluent hypergeometric limit function.
 */
template<typename _Tp>
  void
  run_conf_hyperg_lim(_Tp __c, _Tp __z, int __n = 40)
  {
    auto __s = __b;
    for (int __k = __n; __k >= 2; --__k)
      {
        auto __r = __z / (_Tp(__k) * (__c + __k));
        auto __a = -__r;
	auto __b = _Tp{1} + __r;
        __s = __b + __a / __s;
      }
      __s += __z / __c;
      __s += _Tp{1};
    return __s;
  }


int
main()
{
  std::cout << "\nfloat\n=====\n\n";
  run_conf_hyperg_lim<float>(0.5f, 0.25f);

  std::cout << "\ndouble\n======\n";
  run_conf_hyperg_lim<double>(0.5, 0.25);

  std::cout << "\nlong double\n===========\n";
  run_conf_hyperg_lim<long double>(0.5l, 0.25l);

  std::cout << "\n__float128\n==========\n";
  run_conf_hyperg_lim<__float128>(0.5q, 0.25q);
}
