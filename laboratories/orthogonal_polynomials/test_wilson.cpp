/**
 *
 */

#include <stdexcept>
#include <iostream>
#include <iomanip>
#include <limits>
#include <cmath>

template<typename Tp>
  struct
  wilson_t
  {
    Tp value;
    Tp factor;
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
template<typename Tp, typename TpX>
  wilson_t<Tp>
  wilson_recur(int n, Tp a, Tp b, Tp c, Tp d, TpX x)
  {
    auto Wnm1 = Tp{1};
    if (n == 0)
      return {Wnm1, Tp{1}};

    const auto abcd = a + b + c + d;
    const auto aa = a * a;
    const auto ab = a + b;
    const auto ac = a + c;
    const auto ad = a + d;
    const auto bc = b + c;
    const auto cd = c + d;
    const auto db = d + b;
    const auto xx = Tp(x * x);

    auto fact = Tp{1};

    auto fn = ab * ac * ad;
    fact *= fn;
    auto Wn = Tp{1} - abcd * (aa + xx) / fn;
    if (n == 1)
      return {Wn, fact};

    fn = (ab + Tp{1}) * (ac + Tp{1}) * (ad + Tp{1});
    fact *= fn;
    auto An = abcd * fn / (abcd + Tp{1}) / (abcd + Tp{2});
    auto Cn = bc * cd * db / abcd / (abcd + Tp{1});

    auto Wnp1 = ((An + Cn - (aa + xx)) * Wn - Cn * Wnm1) / An;

    for (int k = 2; k < n; ++k)
      {
	const auto abcd2n = abcd + Tp(2 * k);

	fn = (ab + Tp(k)) * (ac + Tp(k)) * (ad + Tp(k));
	fact *= fn;

	auto An = (abcd + Tp(k - 1)) * fn / (abcd2n - Tp(1)) / abcd2n;
	auto Cn = Tp(k)
		* (bc + Tp(k - 1)) * (cd + Tp(k - 1)) * (db + Tp(k - 1))
		/ (abcd2n - Tp(2)) / (abcd2n - Tp(1));

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
template<typename Tp, typename TpX>
  wilson_t<Tp>
  wilson(int n, Tp a, Tp b, Tp c, Tp d, TpX x)
  {
    if (std::isnan(a))
      return {a, Tp{}};
    else if (std::isnan(b))
      return {b, Tp{}};
    else if (std::isnan(c))
      return {c, Tp{}};
    else if (std::isnan(d))
      return {d, Tp{}};
    else if (std::isnan(x))
      return {x, Tp{}};
    else
      return wilson_recur(n, a, b, c, d, x);
  }

/**
 * 
 */
template<typename Tp>
  void
  test_wilson(int n_max, Tp a, Tp b, Tp c, Tp d)
  {
    std::cout.precision(std::numeric_limits<Tp>::digits10);
    auto w = std::cout.precision() + 8;

    for (int n = 0; n <= n_max; ++n)
      {
	std::cout << '\n' << '\n' << " n = " << n << '\n';
	for (int i = 0; i <= 400; ++i)
	  {
	    auto x = i * Tp{0.05L};
	    auto W = wilson(n, a, b, c, d, std::sqrt(x));
	    std::cout << ' ' << std::setw(w) << x
		      << ' ' << std::setw(w) << W.value
		      << '\n';
	  }
      }
  }

int
main()
{
  test_wilson(10, 1.0f, 1.0f, 1.0f, 1.0f);
}
