
#include <limits>
#include <iostream>
#include <iomanip>
#include <vector>
#include <complex>
#include <emsr/math_constants.h>
#include <ext/float128_io.h>
#include <emsr/numeric_limits.h>
#include <emsr/fp_type_util.h>
#include <emsr/complex_util.h>
//#include "rational.h"
#include <emsr/polynomial.h>

/**
 * Convert a series to a continued fraction.
 * An example using the confluent hypergeometric limit function.
 */
template<typename Tp>
  void
  run_conf_hyperg_lim(Tp c, int n = 20)
  {
    std::cout.precision(std::numeric_limits<Tp>::digits10);
    auto width = 8 + std::cout.precision();

    for (int i = 0; i <= 100; ++i)
      {
        auto z = Tp{0.01Q} * i;
	auto f = Tp{1} + z / (Tp(n) * (c + n)); // b;
	for (int k = n - 1; k >= 2; --k)
	  {
            auto r = z / (Tp(k) * (c + k - 1));
            auto a = -r;
	    auto b = Tp{1} + r;
            f = b + a / f;
	  }
	f += z / c;
	f += Tp{1};

	auto t = Tp{1};
	auto s = t;
	for (int k = 1; k <= n; ++k)
	  {
	    t *= z / (Tp(k) * (c + k - 1));
	    s += t;
	  }

	std::cout << ' ' << std::setw(width) << z
		  << ' ' << std::setw(width) << s
		  << ' ' << std::setw(width) << f
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

  //std::cout << "\n__float128\n==========\n";
  //run_conf_hyperg_lim<__float128>(0.5q);
}
