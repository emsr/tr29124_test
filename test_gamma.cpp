/*
$HOME/bin_tr29124/bin/g++ -std=gnu++17 -I. -o test_gamma test_gamma.cpp wrap_boost.cpp -lquadmath
LD_LIBRARY_PATH=$HOME/bin_tr29124/lib64:$LD_LIBRARY_PATH ./test_gamma > test_gamma.txt

$HOME/bin/bin/g++ -std=gnu++17 -I. -o test_gamma test_gamma.cpp wrap_boost.cpp -lquadmath
./test_gamma > test_gamma.txt

g++ -std=gnu++17 -DNO_LOGBQ -I. -o test_gamma test_gamma.cpp wrap_boost.cpp -lquadmath
./test_gamma > test_gamma.txt
*/

#include <bits/specfun.h>
#include <bits/float128.h>
#include <ext/math_const.h>
#include <limits>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <vector>
#include <string>
#include <complex>

#include "wrap_boost.h"

// I'm not sure why I need this here and not other places...
template<>
  constexpr std::array<float, 7>
  std::__detail::_GammaSpouge<float>::_S_cheby;
template<>
  constexpr std::array<double, 18>
  std::__detail::_GammaSpouge<double>::_S_cheby;
template<>
  constexpr std::array<long double, 22>
  std::__detail::_GammaSpouge<long double>::_S_cheby;

template<>
  constexpr std::array<float, 7>
  std::__detail::_GammaLanczos<float>::_S_cheby;
template<>
  constexpr std::array<double, 10>
  std::__detail::_GammaLanczos<double>::_S_cheby;
template<>
  constexpr std::array<long double, 11>
  std::__detail::_GammaLanczos<long double>::_S_cheby;

template<typename _Tp, typename _Gamma>
  void
  test_gamma(_Gamma gamma)
  {
    using _Val = _Tp;
    using _Real = std::__detail::__num_traits_t<_Val>;

    std::cout.precision(std::numeric_limits<_Real>::digits10);
    std::cout << std::showpoint << std::scientific;
    auto width = 8 + std::cout.precision();

    std::cout << '\n'
	      << ' ' << std::setw(width) << "a"
	      << ' ' << std::setw(width) << "log_gamma"
	      << ' ' << std::setw(width) << "std::lgamma"
	      << ' ' << std::setw(width) << "__log_gamma"
	      << ' ' << std::setw(width) << "delta_std"
	      << ' ' << std::setw(width) << "delta_boost"
	      << '\n';
    int i_min = -200;
    for (int i = i_min; i <= +500; ++i)
      {
	auto a = _Tp{0.10L} * i;
	auto gam = gamma(a - _Tp{1});
	auto gam0 = std::lgamma(a);
	auto glgam = std::__detail::__log_gamma(a);
	auto blgam = beast::lgamma(a);
	std::cout << ' ' << std::setw(width) << a
		  << ' ' << std::setw(width) << gam
		  << ' ' << std::setw(width) << gam0
		  << ' ' << std::setw(width) << glgam
		  << ' ' << std::setw(width) << (gam - gam0) / std::abs(gam0)
		  << ' ' << std::setw(width) << (glgam - blgam) / std::abs(blgam)
		  << '\n';
      }
  }

int
main()
{

  std::cout << "\n\nLanczos Algorithm\n\n";
  test_gamma<double>(std::__detail::__log_gamma1p_lanczos<double>);

  std::cout << "\n\nSpouge Algorithm\n\n";
  test_gamma<double>(std::__detail::__log_gamma1p_spouge<double>);
}
