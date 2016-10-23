/*
$HOME/bin_tr29124/bin/g++ -std=gnu++17 -I. -o test_beta test_beta.cpp wrap_boost.cpp -lquadmath
./test_beta > test_beta.txt

$HOME/bin/bin/g++ -std=gnu++17 -I. -o test_beta test_beta.cpp wrap_boost.cpp -lquadmath
./test_beta > test_beta.txt

g++ -std=gnu++17 -DNO_LOGBQ -I. -o test_beta test_beta.cpp wrap_boost.cpp -lquadmath
./test_beta > test_beta.txt
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

template<typename _Tp>
  void
  test_beta()
  {
    using _Val = _Tp;
    using _Real = std::__detail::__num_traits_t<_Val>;

    std::cout.precision(std::numeric_limits<_Real>::digits10);
    std::cout << std::showpoint << std::scientific;
    auto width = 8 + std::cout.precision();

    std::cout << '\n'
	      << ' ' << std::setw(width) << "a"
	      << ' ' << std::setw(width) << "b"
	      << ' ' << std::setw(width) << "__beta"
	      << ' ' << std::setw(width) << "boost::beta"
	      << ' ' << std::setw(width) << "delta_boost"
	      << '\n';
    int i_min = 1;
    for (int i = i_min; i <= +500; ++i)
      {
	auto a = _Tp{0.1L} * i;
	int j_min = 1;
	for (int j = j_min; j <= +500; ++j)
	  {
	    auto b = _Tp{0.1L} * j;
	    auto gbet = std::__detail::__beta(a, b);
	    auto bbet = beast::beta(a, b);
	    std::cout << ' ' << std::setw(width) << a
		      << ' ' << std::setw(width) << b
		      << ' ' << std::setw(width) << gbet
		      << ' ' << std::setw(width) << bbet
		      << ' ' << std::setw(width) << (gbet - bbet) / std::abs(bbet)
		      << '\n';
	  }
      }
  }

template<typename _Tp>
  void
  plot_beta(std::string filename)
  {
    using _Val = _Tp;
    using _Real = std::__detail::__num_traits_t<_Val>;

    auto data = std::ofstream(filename);

    data.precision(std::numeric_limits<_Real>::digits10);
    data << std::showpoint << std::scientific;
    auto width = 8 + data.precision();

    int i_min = -500;
    for (int i = i_min; i <= +500; ++i)
      {
	auto a = _Tp{0.2L} * i;
	int j_min = -500;
	data << '\n';
	for (int j = j_min; j <= +500; ++j)
	  {
	    auto b = _Tp{0.2L} * j;
	    auto gbet = std::__detail::__beta(a, b);
	    data << ' ' << std::setw(width) << a
		 << ' ' << std::setw(width) << b
		 << ' ' << std::setw(width) << gbet
		 << '\n';
	  }
      }
  }

int
main()
{
  test_beta<float>();

  test_beta<double>();

  test_beta<long double>();

  // Beta seems to be either really tiny or really huge.
  // Maybe graph log_beta.
  //plot_beta<double>("plot/beta_double.txt");
}
