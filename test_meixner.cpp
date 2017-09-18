/*
$HOME/bin/bin/g++ -std=c++17 -g -Wall -Wextra -o test_meixner test_meixner.cpp -Lwrappers/debug -lwrap_burkhardt -lgfortran
LD_LIBRARY_PATH=$HOME/bin/lib64:wrappers/debug:$LD_LIBRARY_PATH ./test_meixner > test_meixner.txt
*/

#include <iostream>
#include <iomanip>
#include <limits>
#include <cmath>

#include "wrap_burkhardt.h"

/**
 * 
 */
template<typename _Tp>
  _Tp
  __meixner_recur(int n, _Tp beta, _Tp c, _Tp x)
  {
    auto Mnm2 = _Tp{1};
    if (n == 0)
      return Mnm2;

    auto Mnm1 = _Tp{1} + (_Tp{1} - _Tp{1} / c) * x / beta;
    if (n == 1)
      return Mnm1;

    const auto cc = _Tp{1} - c;
    auto nm1 = 1;
    auto cbnm1 = c * (beta + nm1);
    auto Mn = ((cbnm1 + nm1 - cc * x) * Mnm1 - nm1 * Mnm2) / cbnm1;
    for (int k = 3; k <= n; ++k)
      {
	nm1 = _Tp(k - 1);
	cbnm1 = c * (beta + nm1);
	Mnm2 = Mnm1;
	Mnm1 = Mn;
	Mn = ((cbnm1 + nm1 - cc * x) * Mnm1 - nm1 * Mnm2) / cbnm1;
      }

    return Mn;
  }

/**
 * 
 */
template<typename _Tp>
  _Tp
  __meixner(int n, _Tp beta, _Tp c, _Tp x)
  {
    if (std::isnan(beta))
      return beta;
    if (std::isnan(c))
      return c;
    else
      return __meixner_recur(n, beta, c, x);
  }

template<typename _Tp>
  void
  test_meixner(int n_max, _Tp beta, _Tp c)
  {
    std::cout.precision(std::numeric_limits<_Tp>::digits10);
    auto w = std::cout.precision() + 8;

    for (int n = 0; n <= n_max; ++n)
      {
	std::cout << '\n' << '\n' << " n = " << n << '\n';
	for (int i = 0; i <= 400; ++i)
	  {
	    auto x = i * _Tp{0.05L};
	    auto M = __meixner(n, beta, c, x);
	    auto M_test = burkhardt::meixner(n, beta, c, x);
	    std::cout << ' ' << std::setw(w) << x
		      << ' ' << std::setw(w) << M
		      << ' ' << std::setw(w) << M_test
		      << '\n';
	  }
      }
  }

int
main()
{
  test_meixner(10, 0.5f, 0.5f);
}
