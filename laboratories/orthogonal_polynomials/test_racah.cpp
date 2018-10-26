/*
$HOME/bin/bin/g++ -std=c++17 -g -Wall -Wextra -o test_racah test_racah.cpp
./test_racah > test_racah.txt
*/

#include <stdexcept>
#include <iostream>
#include <iomanip>
#include <limits>
#include <cmath>

/**
 * Compute the Racah polynomial by recursion:
 * @f[
 *   \lambda(x)R_n(\lambda(x)) = A_n R_{n+1}(\lambda(x))
 *                               - (A_n + C_n) R_n(\lambda(x))
 *                               + C_n R_{n-1}(\lambda(x))
 * @f]
 * where
 * @f[
 *  A_n = \frac{(n + \alpha + 1)(n + \alpha + \beta + 1)
 *              (n + \beta + \delta + 1)(n + \gamma + 1)}
 *             {(2n + \alpha + \beta + 1)(2n + \alpha + \beta + 2)}
 * @f]
 * and
 * @f[
 *  C_n = \frac{n(n + \alpha + \beta - \gamma)(n + \alpha - \delta)(n + \beta)}
 *             {(2n + \alpha + \beta)(2n + \alpha + \beta + 1)}
 * @f]
 */
template<typename _Tp, typename _TpX>
  _Tp
  __racah_recur(int n, _Tp a, _Tp b, _Tp c, _Tp d, _TpX x)
  {
    auto Rnm1 = _Tp{1};
    if (n == 0)
      return Rnm1;

    auto lambda = _Tp(x) * (_Tp(x) + c + d + _Tp{1});

    auto Rn = _Tp{1} + (a + b + _Tp{2}) * lambda
		     / (a + _Tp{1}) / (b + d + _Tp{1}) / (c + _Tp{1});
    if (n == 1)
      return Rn;

    auto An = (a + _Tp{3}) * (a + b + _Tp{3}) * (b + d + _Tp{3}) * (c + _Tp{3})
	    / (a + b + _Tp{5}) / (a + b + _Tp{6});
    auto Cn = _Tp{2} * (a + b - c + _Tp{2}) * (a - d + _Tp{2}) * (b + _Tp{2})
		/ (a + b + _Tp{4}) / (a + b + _Tp{5});
    auto Rnp1 = ((An + Cn + lambda) * Rn - Cn * Rnm1) / An;

    for (int k = 2; k < n; ++k)
      {
	auto An = (a + _Tp(k + 1)) * (a + b + _Tp(k + 1))
		* (b + d + _Tp(k + 1)) * (c + _Tp(k + 1))
		/ (a + b + _Tp(2 * k + 1)) / (a + b + _Tp(2 * k + 2));

	auto Cn = _Tp(k) * (a + b - c + _Tp(k))
		* (a - d + _Tp(k)) * (b + _Tp(k))
		/ (a + b + _Tp(2 * k)) / (a + b + _Tp(2 * k + 1));

	Rnm1 = Rn;
	Rn = Rnp1;
        Rnp1 = ((An + Cn + lambda) * Rn - Cn * Rnm1) / An;
      }

    return Rnp1;
  }

/**
 * Return the Racah polynomial defined by
 * @f[
 *    R_n(\lambda(x), \alpha, \beta, \gamma, \delta)
 *      = {}_4F_3(-n, n + \alpha + \beta + 1, -x, x + \gamma + \delta + 1;
 *                \alpha + 1, \beta + \delta + 1, \gamma + 1; 1)
 * @f]
 * where @f$ \lambda(x) = x(x + \gamma + \delta + 1) @f$.
 */
template<typename _Tp, typename _TpX>
  _Tp
  __racah(int n, _Tp a, _Tp b, _Tp c, _Tp d, _TpX x)
  {
    if (std::isnan(a))
      return a;
    else if (std::isnan(b))
      return b;
    else if (std::isnan(c))
      return c;
    else if (std::isnan(d))
      return d;
    else if (std::isnan(x))
      return x;
    else
      return __racah_recur(n, a, b, c, d, x);
  }

/**
 * 
 */
template<typename _Tp>
  void
  test_racah(int n_max, _Tp a, _Tp b, _Tp c, _Tp d)
  {
    std::cout.precision(std::numeric_limits<_Tp>::digits10);
    auto w = std::cout.precision() + 8;

    auto lambda = [c, d](_Tp x)
		  { return x * (x + c + d + _Tp{1}); };

    for (int n = 0; n <= n_max; ++n)
      {
	std::cout << '\n' << '\n' << " n = " << n << '\n';
	for (int i = 0; i <= 200; ++i)
	  {
	    auto x = i * _Tp{0.05L};
	    auto R = __racah(n, a, b, c, d, x);
	    std::cout << ' ' << std::setw(w) << lambda(x)
		      << ' ' << std::setw(w) << R
		      << '\n';
	  }
      }
  }

int
main()
{
  test_racah(10, -11.0f, 18.0f, 4.0f, 3.0f);
}
