/**
 *
 */

#include <iostream>
#include <iomanip>
#include <limits>
#include <cmath>

#include <emsr/special_functions.h>
#include <wrap_burkhardt.h>

/**
 * Compute the Charlier polynomial by recursion:
 * @f[
 *    -x C_n(x) = a C_{n+1}(x) - (n + a) C_n(x) + n C_{n-1}(x)
 * @f]
 * where @f$ C_n(x) = C_n(x; a) @f$.
 */
template<typename Tp, typename TpX>
  Tp
  charlier_recur(int n, Tp a, TpX x)
  {
    auto Cnm1 = Tp{1};
    if (n == 0)
      return Cnm1;

    auto Cn = Tp{1} - Tp(x) / a;
    if (n == 1)
      return Cn;

    auto Cnp1 = ((Tp{1} + a - Tp(x)) * Cn - Cnm1) / a;
    for (int k = 2; k < n; ++k)
      {
	Cnm1 = Cn;
	Cn = Cnp1;
	Cnp1 = ((Tp(k) + a - Tp(x)) * Cn - Tp(k) * Cnm1) / a;
      }

    return Cnp1;
  }

/**
 * Return the Charlier polynomial defined by
 * @f[
 *    C_n(x; a) = {}_2F_0(-n, -x; ; -\frac{1}{a})
 * @f]
 */
template<typename Tp, typename TpX>
  Tp
  charlier(int n, Tp a, TpX x)
  {
    if (std::isnan(a))
      return a;
    else if (std::isnan(x))
      return x;
    else
      return charlier_recur(n, a, x);
  }

/**
 * 
 */
template<typename Tp>
  void
  test_charlier(int n_max, Tp a)
  {
    std::cout.precision(std::numeric_limits<Tp>::digits10);
    auto w = std::cout.precision() + 8;

    for (int n = 0; n <= n_max; ++n)
      {
	std::cout << '\n' << '\n' << " n = " << n << '\n';
	for (int i = 0; i <= 400; ++i)
	  {
	    auto x = i * Tp{0.05L};
	    auto C = charlier(n, a, x);
	    auto C_test = burkhardt::charlier(n, a, x);
	    std::cout << ' ' << std::setw(w) << x
		      << ' ' << std::setw(w) << C
		      << ' ' << std::setw(w) << C_test
		      << '\n';
	  }
      }
  }

int
main()
{
  test_charlier<float>(10, 2.0f);
}
