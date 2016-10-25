/*
$HOME/bin_tr29124/bin/g++ -std=gnu++17 -g -Wall -Wextra -I. -o test_expint test_expint.cpp wrap_boost.cpp -lquadmath
./test_expint > test_expint.txt

$HOME/bin/bin/g++ -std=gnu++17 -g -Wall -Wextra -I. -o test_expint test_expint.cpp wrap_boost.cpp -lquadmath
./test_expint > test_expint.txt

g++ -std=gnu++17 -DNO_LOGBQ -g -Wall -Wextra -I. -o test_expint test_expint.cpp wrap_boost.cpp -lquadmath
./test_expint > test_expint.txt
*/

#include <cmath>
#include <limits>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <vector>
#include <string>
#include <complex>
#include <bits/float128_io.h>

#include "wrap_boost.h"

template<typename _Tp>
  void
  test_expint()
  {
    using _Val = _Tp;
    using _Real = std::__detail::__num_traits_t<_Val>;

    std::cout.precision(std::numeric_limits<_Real>::digits10);
    std::cout << std::showpoint << std::scientific;
    auto width = 8 + std::cout.precision();
    std::vector<unsigned int> order({0, 1, 2, 5, 10, 20, 50, 100});

    for (auto n : order)
      {
	std::cout << '\n'
		  << ' ' << std::setw(width) << "n"
		  << ' ' << std::setw(width) << "x"
		  << ' ' << std::setw(width) << "log_gamma"
		  << ' ' << std::setw(width) << "std::lgamma"
		  << ' ' << std::setw(width) << "__log_gamma"
		  << ' ' << std::setw(width) << "delta_std"
		  << ' ' << std::setw(width) << "delta_boost"
		  << '\n';
	int i_min = 0;
	for (int i = i_min; i <= +500; ++i)
	  {
	    auto x = 0.10L * i;
	    auto ens = std::__detail::__expint_En_series(n, x);
	    auto enc = std::__detail::__expint_En_cont_frac(n, x);
	    auto enb = beast::expint(n, x);
	    std::cout << ' ' << std::setw(width) << n
		      << ' ' << std::setw(width) << x
		      << ' ' << std::setw(width) << ens
		      << ' ' << std::setw(width) << enc
		      << ' ' << std::setw(width) << (ens - enb) / std::abs(enb)
		      << ' ' << std::setw(width) << (enc - enb) / std::abs(enb)
		      << '\n';
	  }
      }
  }

int
main()
{
  test_expint<double>();
}
