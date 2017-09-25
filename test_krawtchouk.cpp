/*
$HOME/bin/bin/g++ -std=c++17 -g -Wall -Wextra -o test_krawtchouk test_krawtchouk.cpp -Lwrappers/debug -lwrap_burkhardt -lgfortran
LD_LIBRARY_PATH=$HOME/bin/lib64:wrappers/debug:$LD_LIBRARY_PATH ./test_krawtchouk > test_krawtchouk.txt
*/

#include <stdexcept>
#include <iostream>
#include <iomanip>
#include <limits>
#include <cmath>

#include "wrap_burkhardt.h"

/**
 * 
 */
template<typename _Tp, typename _TpX>
  _Tp
  __krawtchouk_recur(int n, _Tp p, int N, _TpX x)
  {
    auto Knm2 = _Tp{1};
    if (n == 0)
      return Knm2;

    auto Knm1 = _Tp{1} - _Tp(x) / p / _Tp(N);
    if (n == 1)
      return Knm1;

    const auto q = _Tp{1} - p;
    auto pnn = p * _Tp(N - 1);
    auto nm1q = _Tp{1} * q;
    auto Kn = ((pnn + nm1q - _Tp(x)) * Knm1 - nm1q * Knm2) / pnn;
    for (int k = 3; k <= n; ++k)
      {
	pnn = p * _Tp(N - k + 1);
	nm1q = _Tp(k - 1) * q;
	Knm2 = Knm1;
	Knm1 = Kn;
	Kn = ((pnn + nm1q - _Tp(x)) * Knm1 - nm1q * Knm2) / pnn;
      }

    return Kn;
  }

/**
 * 
 */
template<typename _Tp, typename _TpX>
  _Tp
  __krawtchouk(int n, _Tp p, int N, _TpX x)
  {
    if (std::isnan(p))
      return p;
    else if (std::isnan(x))
      return x;
    else if (n > N)
      std::__throw_domain_error(__N("__krawtchouk: "
				"Degree must be less than or equal to big N."));
    else
      return __krawtchouk_recur(n, p, N, x);
  }

/**
 * 
 */
template<typename _Tp>
  void
  test_krawtchouk(_Tp p, int N)
  {
    std::cout.precision(std::numeric_limits<_Tp>::digits10);
    auto w = std::cout.precision() + 8;

    for (int n = 0; n <= N; ++n)
      {
	std::cout << '\n' << '\n' << " n = " << n << '\n';
	for (int i = 0; i <= 400; ++i)
	  {
	    auto x = i * _Tp{0.05L};
	    auto K = __krawtchouk(n, p, N, x);
	    auto K_test = burkhardt::krawtchouk(n, p, N, x);
	    std::cout << ' ' << std::setw(w) << x
		      << ' ' << std::setw(w) << K
		      << ' ' << std::setw(w) << K_test
		      << '\n';
	  }
      }
  }

int
main()
{
  test_krawtchouk(0.5f, 10);
}
