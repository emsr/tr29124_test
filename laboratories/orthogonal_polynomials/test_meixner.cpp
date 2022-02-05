/**
 *
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
template<typename Tp, typename TpX>
  Tp
  meixner_recur(int n, Tp beta, Tp c, TpX x)
  {
    auto Mnm1 = Tp{1};
    if (n == 0)
      return Mnm1;

    auto Mn = Tp{1} + (Tp{1} - Tp{1} / c) * Tp(x) / beta;
    if (n == 1)
      return Mn;

    const auto cc = Tp{1} - c;
    auto nn = Tp{1};
    auto cbnn = c * (beta + nn);
    auto Mnp1 = ((cbnn + nn - cc * Tp(x)) * Mn - nn * Mnm1) / cbnn;

    for (int k = 2; k < n; ++k)
      {
	nn = Tp(k);
	cbnn = c * (beta + nn);
	Mnm1 = Mn;
	Mn = Mnp1;
	Mnp1 = ((cbnn + nn - cc * Tp(x)) * Mn - nn * Mnm1) / cbnn;
      }

    return Mnp1;
  }

/**
 * Return the Meixner polynomial defined by
 * @f[
 *    M_n(x; \beta, c) = {}_2F_1(-n, -x; \beta; 1 - \frac{1}{c})
 * @f]
 */
template<typename Tp, typename TpX>
  Tp
  meixner(int n, Tp beta, Tp c, TpX x)
  {
    if (std::isnan(beta))
      return beta;
    if (std::isnan(c))
      return c;
    if (std::isnan(x))
      return x;
    else
      return meixner_recur(n, beta, c, x);
  }

template<typename Tp>
  void
  test_meixner(int n_max, Tp beta, Tp c)
  {
    std::cout.precision(std::numeric_limits<Tp>::digits10);
    auto w = std::cout.precision() + 8;

    for (int n = 0; n <= n_max; ++n)
      {
	std::cout << '\n' << '\n' << " n = " << n << '\n';
	for (int i = 0; i <= 400; ++i)
	  {
	    auto x = i * Tp{0.05L};
	    auto M = meixner(n, beta, c, x);
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
