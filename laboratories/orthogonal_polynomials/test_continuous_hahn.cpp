/**
 *
 */

#include <stdexcept>
#include <iostream>
#include <iomanip>
#include <limits>
#include <cmath>
#include <complex>

template<typename Tp>
  struct
  continuous_hahn_t
  {
    std::complex<Tp> value;
    std::complex<Tp> factor;
  };

/**
 * Compute the continuous Hahn polynomial by recursion:
 * @f[
 *    (a + ix)p_n(x) = A_n p_{n+1}(x) - (A_n + C_n) p_n(x) + C_n p_{n-1}(x)
 * @f]
 * where @f$ p_n(x) = p_n(x; a, b, c, d) @f$ and
 * @f[
 *    A_n = \frac{(n + a + b + c + d - 1)(n + a + c)(n + a + d)}
 *               {(2n + a + b + c + d - 1)(2n + a + b + c + d)}
 * @f]
 * and
 * @f[
 *    C_n = \frac{n(n + b + c - 1)(n + b + d - 1)}
 *               {(2n + a + b + c + d - 2)(2n + a + b + c + d - 1)}
 * @f]
 */
template<typename Tp, typename TpX>
  continuous_hahn_t<Tp>
  continuous_hahn_recur(int n, Tp a, Tp b, Tp c, Tp d, TpX x)
  {
    auto pnm1 = std::complex<Tp>{1};
    if (n == 0)
      return {pnm1, Tp{1}};

    constexpr std::complex<Tp> i{0, 1};
    const auto abcd = a + b + c + d;
    const auto ac = a + c;
    const auto ad = a + d;
    const auto bc = b + c;
    const auto db = d + b;
    const auto ix = i * Tp(x);

    auto fact = std::complex<Tp>{1};

    auto fn = ac * ad;
    fact *= i * fn;
    auto pn = Tp{1} - abcd * (a + ix) / fn;
    if (n == 1)
      return {pn, fact};

    fn = (ac + Tp{1}) * (ad + Tp{1});
    fact *= i * fn / Tp{2};
    auto An = -abcd * fn / (abcd + Tp{1}) / (abcd + Tp{2});
    auto Cn = bc * db / abcd / (abcd + Tp{1});

    auto pnp1 = ((An + Cn + (a + ix)) * pn - Cn * pnm1) / An;

    for (int k = 2; k < n; ++k)
      {
	const auto abcd2n = abcd + Tp(2 * k);

	fn = (ac + Tp(k)) * (ad + Tp(k));
	fact *= i * fn / Tp(1 + k);

	auto An = -(abcd + Tp(k - 1)) * fn / (abcd2n - Tp(1)) / abcd2n;
	auto Cn = Tp(k)
		* (bc + Tp(k - 1)) * (db + Tp(k - 1))
		/ (abcd2n - Tp(2)) / (abcd2n - Tp(1));

	pnm1 = pn;
	pn = pnp1;
        pnp1 = ((An + Cn + (a + ix)) * pn - Cn * pnm1) / An;
      }

    return {pnp1, fact};
  }

/**
 * Compute the continuous Hahn polynomial defined by
 * @f[
 *    p_n(x; a, b, c, d) = i^n\frac{(a+c)_n(a+d)_n}{n!}
 *       {}_3F_2(-n, n + a + b + c + d - 1, a + ix;
 *               a + c, a + d; 1)
 * @f]
 */
template<typename Tp, typename TpX>
  continuous_hahn_t<Tp>
  continuous_hahn(int n, Tp a, Tp b, Tp c, Tp d, TpX x)
  {
    if (std::isnan(a))
      return {std::complex<Tp>(a), std::complex<Tp>{}};
    else if (std::isnan(b))
      return {std::complex<Tp>(b), std::complex<Tp>{}};
    else if (std::isnan(c))
      return {std::complex<Tp>(c), std::complex<Tp>{}};
    else if (std::isnan(d))
      return {std::complex<Tp>(d), std::complex<Tp>{}};
    else if (std::isnan(x))
      return {std::complex<Tp>(x), std::complex<Tp>{}};
    else
      return continuous_hahn_recur(n, a, b, c, d, x);
  }

template<typename Tp>
  void
  test_continuous_hahn(int n_max, Tp a, Tp b, Tp c, Tp d)
  {
    std::cout.precision(std::numeric_limits<Tp>::digits10);
    auto w = std::cout.precision() + 8;

    for (int n = 0; n <= n_max; ++n)
      {
	std::cout << '\n' << '\n' << " n = " << n << '\n';
	for (int i = -200; i <= 200; ++i)
	  {
	    auto x = i * Tp{0.01L};
	    auto p = continuous_hahn(n, a, b, c, d, x);
	    std::cout << ' ' << std::setw(w) << x
		      << ' ' << std::setw(w) << std::real(p.factor * p.value)
		      << ' ' << std::setw(w) << std::imag(p.factor * p.value)
		      << '\n';
	  }
      }
  }

int
main()
{
  test_continuous_hahn(10, 1.0f, 2.0f, 2.0f, 1.0f);
}
