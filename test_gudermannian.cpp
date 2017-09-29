/*
$HOME/bin/bin/g++ -std=c++17 -g -I. -Wall -Wextra -o test_gudermannian test_gudermannian.cpp
./test_gudermannian > test_gudermannian.txt
*/

#include <iostream>
#include <iomanip>
#include <limits>
#include <cmath>

#include "ext/math_const.h"

template<typename _Tp>
  _Tp
  gd_series(_Tp x)
  {
    const auto xx = x * x;
    auto term = x;
    auto sum = term;
    for (int i = 1; i < 10; ++i)
      {
	term *= xx * __gnu_cxx::euler<_Tp>(2 * i)
	     / _Tp(2 * i) / _Tp(2 * i + 1);
	sum += term;
      }
    return sum;
  }

template<typename _Tp>
  _Tp
  gd_asin(_Tp x)
  {
    return std::asin(std::tanh(x));
  }

template<typename _Tp>
  _Tp
  gd(_Tp x)
  {
    const auto _S_pi = __const_2_pi(x);

    if (std::isnan(x))
      return x;
    else if (x == -std::numeric_limits<_Tp>::infinity())
      return -_S_pi / _Tp{2};
    else if (x == +std::numeric_limits<_Tp>::infinity())
      return +_S_pi / _Tp{2};
    else
      return gd_asin(x);
  }

template<typename _Tp>
  void
  test_gudermannian()
  {
    std::cout.precision(std::numeric_limits<_Tp>::digits10);
    auto w = std::cout.precision() + 8;

    const auto _S_pi = __gnu_cxx::__const_2_pi<_Tp>();
    const auto del = _S_pi / 400;
    for (int i = -200; i <= 200; ++i)
      {
	auto x = i * del;
	auto gd_trig = gd_asin(x);
	auto gd_ser = gd_series(x);
	auto gd_diff = gd_ser - gd_trig;
	std::cout << ' ' << std::setw(w) << x
		  << ' ' << std::setw(w) << gd_ser
		  << ' ' << std::setw(w) << gd_trig
		  << ' ' << std::setw(w) << gd_diff
		  << '\n';
      }
  }

int
main()
{
  test_gudermannian<float>();
}
