/*
$HOME/bin/bin/g++ -std=c++2a -g -Wall -Wextra -o test_krawtchouk test_krawtchouk.cpp -Lwrappers/debug -lwrap_burkhardt -lgfortran
LD_LIBRARY_PATH=$HOME/bin/lib64:wrappers/debug:$LD_LIBRARY_PATH ./test_krawtchouk > test_krawtchouk.txt
*/

#include <stdexcept>
#include <iostream>
#include <iomanip>
#include <limits>
#include <cmath>

#include "wrap_burkhardt.h"

/**
 * Compute the Krawtchouk polynomial by recursion:
 * @f[
 *    -xK_n(x) = p(N - n) K_{n+1}(x)
 *             - [p(N - n) + n(1 - p)] K_n(x)
 *             + n(1 - p) K_{n-1}(x)
 * @f]
 * where @f$ K_n(x) = K_n(x; p, N) @f$.
 */
template<typename _Tp, typename _TpX>
  _Tp
  __krawtchouk_recur(int n, _Tp p, int N, _TpX x)
  {
    auto Knm1 = _Tp{1};
    if (n == 0)
      return Knm1;

    auto Kn = _Tp{1} - _Tp(x) / p / _Tp(N);
    if (n == 1)
      return Kn;

    const auto q = _Tp{1} - p;
    auto pnn = p * _Tp(N - 1);
    auto nq = q;
    auto Knp1 = ((pnn + nq - _Tp(x)) * Kn - nq * Knm1) / pnn;
    for (int k = 2; k < n; ++k)
      {
	pnn = p * _Tp(N - k);
	nq = _Tp(k) * q;
	Knm1 = Kn;
	Kn = Knp1;
	Knp1 = ((pnn + nq - _Tp(x)) * Kn - nq * Knm1) / pnn;
      }

    return Knp1;
  }

/**
 * Return the Krawtchouk polynomial defined by
 * @f[
 *    K_n(x; p, N) = {}_2F_1(-n, -x; -N; \frac{1}{p})
 * @f]
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
