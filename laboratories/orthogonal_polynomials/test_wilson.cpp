/**
 *
 */

#include <stdexcept>
#include <iostream>
#include <iomanip>
#include <limits>
#include <cmath>

template<typename _Tp>
  struct
  __wilson_t
  {
    _Tp __value;
    _Tp __factor;
  };

/**
 * Compute the Wilson polynomial by recursion:
 * @f[
 *    -(a^2 + x^2)\bar{W}_n(x^2) = A_n \bar{W}_{n+1}(x^2)
 *                               - (A_n + C_n) \bar{W}_n(x^2)
 *                               + C_n \bar{W}_{n-1}(x^2)
 * @f]
 * where  and
 * @f[
 *    A_n = \frac{(n + a + b + c + d - 1)(n + a + b)(n + a + c)(n + a + d)}
 *               {(2n + a + b + c + d - 1)(2n + a + b + c + d)}
 * @f]
 * and
 * @f[
 *    C_n = \frac{n(n + b + c - 1)(n + c + d - 1)(n + d + b - 1)}
 *               {(2n + a + b + c + d - 2)(2n + a + b + c + d - 1)}
 * @f]
 */
template<typename _Tp, typename _TpX>
  __wilson_t<_Tp>
  __wilson_recur(int n, _Tp a, _Tp b, _Tp c, _Tp d, _TpX x)
  {
    auto Wnm1 = _Tp{1};
    if (n == 0)
      return {Wnm1, _Tp{1}};

    const auto abcd = a + b + c + d;
    const auto aa = a * a;
    const auto ab = a + b;
    const auto ac = a + c;
    const auto ad = a + d;
    const auto bc = b + c;
    const auto cd = c + d;
    const auto db = d + b;
    const auto xx = _Tp(x * x);

    auto fact = _Tp{1};

    auto fn = ab * ac * ad;
    fact *= fn;
    auto Wn = _Tp{1} - abcd * (aa + xx) / fn;
    if (n == 1)
      return {Wn, fact};

    fn = (ab + _Tp{1}) * (ac + _Tp{1}) * (ad + _Tp{1});
    fact *= fn;
    auto An = abcd * fn / (abcd + _Tp{1}) / (abcd + _Tp{2});
    auto Cn = bc * cd * db / abcd / (abcd + _Tp{1});

    auto Wnp1 = ((An + Cn - (aa + xx)) * Wn - Cn * Wnm1) / An;

    for (int k = 2; k < n; ++k)
      {
	const auto abcd2n = abcd + _Tp(2 * k);

	fn = (ab + _Tp(k)) * (ac + _Tp(k)) * (ad + _Tp(k));
	fact *= fn;

	auto An = (abcd + _Tp(k - 1)) * fn / (abcd2n - _Tp(1)) / abcd2n;
	auto Cn = _Tp(k)
		* (bc + _Tp(k - 1)) * (cd + _Tp(k - 1)) * (db + _Tp(k - 1))
		/ (abcd2n - _Tp(2)) / (abcd2n - _Tp(1));

	Wnm1 = Wn;
	Wn = Wnp1;
        Wnp1 = ((An + Cn - (aa + xx)) * Wn - Cn * Wnm1) / An;
      }

    return {Wnp1, fact};
  }

/**
 * Return the Wilson polynomial defined by
 * @f[
 *    \frac{W_n(x^2, a, b, c, d)}{(a + b)_n(a + c)_n(a + d)_n}
 *      = \bar{W}_{n}(x^2, a, b, c, d)
 *      = {}_4F_3(-n, n + a + b + c + d, a + ix, a - ix;
 *                                   a + b, a + c, a + d; 1)
 * @f]
 */
template<typename _Tp, typename _TpX>
  __wilson_t<_Tp>
  __wilson(int n, _Tp a, _Tp b, _Tp c, _Tp d, _TpX x)
  {
    if (std::isnan(a))
      return {a, _Tp{}};
    else if (std::isnan(b))
      return {b, _Tp{}};
    else if (std::isnan(c))
      return {c, _Tp{}};
    else if (std::isnan(d))
      return {d, _Tp{}};
    else if (std::isnan(x))
      return {x, _Tp{}};
    else
      return __wilson_recur(n, a, b, c, d, x);
  }

/**
 * 
 */
template<typename _Tp>
  void
  test_wilson(int n_max, _Tp a, _Tp b, _Tp c, _Tp d)
  {
    std::cout.precision(std::numeric_limits<_Tp>::digits10);
    auto w = std::cout.precision() + 8;

    for (int n = 0; n <= n_max; ++n)
      {
	std::cout << '\n' << '\n' << " n = " << n << '\n';
	for (int i = 0; i <= 400; ++i)
	  {
	    auto x = i * _Tp{0.05L};
	    auto W = __wilson(n, a, b, c, d, std::sqrt(x));
	    std::cout << ' ' << std::setw(w) << x
		      << ' ' << std::setw(w) << W.__value
		      << '\n';
	  }
      }
  }

int
main()
{
  test_wilson(10, 1.0f, 1.0f, 1.0f, 1.0f);
}
