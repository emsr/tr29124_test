/*
$HOME/bin_tr29124/bin/g++ -std=gnu++17 -g -Wall -Wextra -I. -o test_continued_fraction test_continued_fraction.cpp -L$HOME/bin/lib64 -lquadmath
./test_continued_fraction > test_continued_fraction.txt

g++ -std=gnu++14 -g -Wall -Wextra -DNO_LOGBQ -I. -o test_continued_fraction test_continued_fraction.cpp -lquadmath
./test_continued_fraction > test_continued_fraction.txt
*/

#include <limits>
#include <iostream>
#include <iomanip>
#include <vector>
#include <complex>
#include <ext/math_const.h>
#include <bits/float128_io.h>
#include <bits/numeric_limits.h>
#include <bits/specfun_util.h>
#include <bits/complex_util.h>
#include <ext/polynomial.h>
//#include "rational.h"

/**
 * Convert a series to a continued fraction.
 * An example using the confluent hypergeometric limit function.
 */
template<typename _Tp>
  void
  run_conf_hyperg_lim(_Tp __c, int __n = 20)
  {
    std::cout.precision(std::numeric_limits<_Tp>::digits10);
    auto width = 8 + std::cout.precision();

    for (int i = 0; i <= 100; ++i)
      {
        auto __z = _Tp{0.01Q} * i;
	auto __f = _Tp{1} + __z / (_Tp(__n) * (__c + __n)); // __b;
	for (int __k = __n - 1; __k >= 2; --__k)
	  {
            auto __r = __z / (_Tp(__k) * (__c + __k - 1));
            auto __a = -__r;
	    auto __b = _Tp{1} + __r;
            __f = __b + __a / __f;
	  }
	__f += __z / __c;
	__f += _Tp{1};

	auto __t = _Tp{1};
	auto __s = __t;
	for (int __k = 1; __k <= __n; ++__k)
	  {
	    __t *= __z / (_Tp(__k) * (__c + __k - 1));
	    __s += __t;
	  }

	std::cout << ' ' << std::setw(width) << __z
		  << ' ' << std::setw(width) << __s
		  << ' ' << std::setw(width) << __f
		  << '\n';
      }
  }


int
main()
{
  std::cout << "\nfloat\n=====\n\n";
  run_conf_hyperg_lim<float>(0.5f);

  std::cout << "\ndouble\n======\n";
  run_conf_hyperg_lim<double>(0.5);

  std::cout << "\nlong double\n===========\n";
  run_conf_hyperg_lim<long double>(0.5l);

  std::cout << "\n__float128\n==========\n";
  run_conf_hyperg_lim<__float128>(0.5q);
}
