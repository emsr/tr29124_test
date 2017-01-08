/*
$HOME/bin_tr29124/bin/g++ -std=gnu++17 -g -Wall -Wextra -Wno-psabi -I. -o test_reperiodized_trig test_reperiodized_trig.cpp wrap_boost.cpp -lquadmath
./test_reperiodized_trig > test_reperiodized_trig.txt

g++ -std=c++14 -o test_reperiodized_trig test_reperiodized_trig.cpp wrap_boost.cpp -lgsl -lgslcblas
./test_reperiodized_trig > test_reperiodized_trig.txt
*/

#include <iostream>
#include <iomanip>
#include <limits>

#include <ext/cmath>

#include "wrap_boost.h"

template<typename _Tp>
  void
  run_sin_cos_pi(_Tp proto = _Tp{})
  {
    const _Tp _S_pi = __gnu_cxx::__const_pi(proto);

    std::cout.precision(__gnu_cxx::__digits10(proto));
    std::cout << std::showpoint << std::scientific;
    auto width = 8 + std::cout.precision();

    std::cout << std::endl;
    std::cout << std::setw(width) << "x"
	      << std::setw(width) << "sin_pi (GCC)"
	      << std::setw(width) << "sin_pi (Boost)"
	      << std::setw(width) << "delta sin_pi"
	      << std::setw(width) << "delta sin(pi x)"
	      << std::setw(width) << "cos_pi (GCC)"
	      << std::setw(width) << "cos_pi (Boost)"
	      << std::setw(width) << "delta cos_pi"
	      << std::setw(width) << "delta cos(pi x)"
	      << std::setw(width) << "tan_pi (GCC)"
	      << std::setw(width) << "delta tan(pi x)"
	      << '\n';
    std::cout << std::setw(width) << "==============="
	      << std::setw(width) << "==============="
	      << std::setw(width) << "==============="
	      << std::setw(width) << "==============="
	      << std::setw(width) << "==============="
	      << std::setw(width) << "==============="
	      << std::setw(width) << "==============="
	      << std::setw(width) << "==============="
	      << std::setw(width) << "==============="
	      << std::setw(width) << "==============="
	      << std::setw(width) << "==============="
	      << '\n';
    const auto del = _Tp{1} / _Tp{10};
    for (int i = -1600; i <= +1600; ++i)
      {
	auto x = del * i;
	auto sin_pi_g = __gnu_cxx::sin_pi(x);
	auto sin_pi_b = beast::sin_pi(x);
	auto cos_pi_g = __gnu_cxx::cos_pi(x);
	auto cos_pi_b = beast::cos_pi(x);
	auto tan_pi_g = __gnu_cxx::tan_pi(x);
	std::cout << std::setw(width) << x
		  << std::setw(width) << sin_pi_g
		  << std::setw(width) << sin_pi_b
		  << std::setw(width) << sin_pi_g - sin_pi_b
		  << std::setw(width) << sin_pi_g - std::sin(_S_pi * x)
		  << std::setw(width) << cos_pi_g
		  << std::setw(width) << cos_pi_b
		  << std::setw(width) << cos_pi_g - cos_pi_b
		  << std::setw(width) << cos_pi_g - std::cos(_S_pi * x)
		  << std::setw(width) << tan_pi_g
		  << std::setw(width) << tan_pi_g - std::tan(_S_pi * x)
		  << '\n';
      }
    std::cout << std::endl;

    std::cout << std::endl;
    std::cout << std::setw(width) << "x"
	      << std::setw(width) << "sin_pi (GCC)"
	      << std::setw(width) << "sin_pi (Boost)"
	      << std::setw(width) << "delta sin_pi"
	      << std::setw(width) << "delta sin(pi x)"
	      << std::setw(width) << "cos_pi (GCC)"
	      << std::setw(width) << "cos_pi (Boost)"
	      << std::setw(width) << "delta cos_pi"
	      << std::setw(width) << "delta cos(pi x)"
	      << std::setw(width) << "tan_pi (GCC)"
	      << std::setw(width) << "delta tan(pi x)"
	      << '\n';
    std::cout << std::setw(width) << "==============="
	      << std::setw(width) << "==============="
	      << std::setw(width) << "==============="
	      << std::setw(width) << "==============="
	      << std::setw(width) << "==============="
	      << std::setw(width) << "==============="
	      << std::setw(width) << "==============="
	      << std::setw(width) << "==============="
	      << std::setw(width) << "==============="
	      << std::setw(width) << "==============="
	      << std::setw(width) << "==============="
	      << '\n';
    const auto bigdel = _Tp{1} / _Tp{4};
    for (int i = 0; i <= +3200; ++i)
      {
	auto x = bigdel * i;
	auto sin_pi_g = __gnu_cxx::sin_pi(x);
	auto sin_pi_b = beast::sin_pi(x);
	auto cos_pi_g = __gnu_cxx::cos_pi(x);
	auto cos_pi_b = beast::cos_pi(x);
	auto tan_pi_g = __gnu_cxx::tan_pi(x);
	std::cout << std::setw(width) << x
		  << std::setw(width) << sin_pi_g
		  << std::setw(width) << sin_pi_b
		  << std::setw(width) << sin_pi_g - sin_pi_b
		  << std::setw(width) << sin_pi_g - std::sin(_S_pi * x)
		  << std::setw(width) << cos_pi_g
		  << std::setw(width) << cos_pi_b
		  << std::setw(width) << cos_pi_g - cos_pi_b
		  << std::setw(width) << cos_pi_g - std::cos(_S_pi * x)
		  << std::setw(width) << tan_pi_g
		  << std::setw(width) << tan_pi_g - std::tan(_S_pi * x)
		  << '\n';
      }
    std::cout << std::endl;
  }


int
main()
{
  std::cout << "\ndouble\n=====\n\n";
  run_sin_cos_pi<double>();
}
