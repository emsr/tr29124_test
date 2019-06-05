/*
$HOME/bin/bin/g++ -std=c++2a -g -Wall -Wextra -o test_meixner test_meixner.cpp -Lwrappers/debug -lwrap_burkhardt -lgfortran
LD_LIBRARY_PATH=$HOME/bin/lib64:wrappers/debug:$LD_LIBRARY_PATH ./test_meixner > test_meixner.txt
*/

#include <iostream>
#include <iomanip>
#include <limits>
#include <cmath>

#include <wrap_burkhardt.h>

/**
 * Compute the Meixner-Pollaczek polynomial by recursion:
 * @f[
 *    (c - 1)x M_n(x) = c(n + \beta) M_{n+1}(x)
 *                    - [n + (n + \beta)c] M_n(x)
 *                    + n M_{n-1}(x)
 * @f]
 * where @f$ M_n(x) = M_n(x; \beta, c) @f$.
 */
template<typename _Tp, typename _TpX>
  _Tp
  __meixner_recur(int n, _Tp beta, _Tp c, _TpX x)
  {
    auto Mnm1 = _Tp{1};
    if (n == 0)
      return Mnm1;

    auto Mn = _Tp{1} + (_Tp{1} - _Tp{1} / c) * _Tp(x) / beta;
    if (n == 1)
      return Mn;

    const auto cc = _Tp{1} - c;
    auto nn = _Tp{1};
    auto cbnn = c * (beta + nn);
    auto Mnp1 = ((cbnn + nn - cc * _Tp(x)) * Mn - nn * Mnm1) / cbnn;

    for (int k = 2; k < n; ++k)
      {
	nn = _Tp(k);
	cbnn = c * (beta + nn);
	Mnm1 = Mn;
	Mn = Mnp1;
	Mnp1 = ((cbnn + nn - cc * _Tp(x)) * Mn - nn * Mnm1) / cbnn;
      }

    return Mnp1;
  }

/**
 * Return the Meixner polynomial defined by
 * @f[
 *    M_n(x; \beta, c) = {}_2F_1(-n, -x; \beta; 1 - \frac{1}{c})
 * @f]
 */
template<typename _Tp, typename _TpX>
  _Tp
  __meixner(int n, _Tp beta, _Tp c, _TpX x)
  {
    if (std::isnan(beta))
      return beta;
    if (std::isnan(c))
      return c;
    if (std::isnan(x))
      return x;
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
  test_meixner(10, 10.f, 0.25f);
}
