/*
$HOME/bin/bin/g++ -std=c++17 -g -Wall -Wextra -o test_continuous_hahn test_continuous_hahn.cpp
./test_continuous_hahn > test_continuous_hahn.txt
*/

#include <stdexcept>
#include <iostream>
#include <iomanip>
#include <limits>
#include <cmath>
#include <complex>

template<typename _Tp>
  struct
  __continuous_hahn_t
  {
    std::complex<_Tp> __value;
    std::complex<_Tp> __factor;
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
template<typename _Tp, typename _TpX>
  __continuous_hahn_t<_Tp>
  __continuous_hahn_recur(int n, _Tp a, _Tp b, _Tp c, _Tp d, _TpX x)
  {
    auto pnm1 = std::complex<_Tp>{1};
    if (n == 0)
      return {pnm1, _Tp{1}};

    constexpr std::complex<_Tp> i{0, 1};
    const auto abcd = a + b + c + d;
    const auto ac = a + c;
    const auto ad = a + d;
    const auto bc = b + c;
    const auto db = d + b;
    const auto ix = i * _Tp(x);

    auto fact = std::complex<_Tp>{1};

    auto fn = ac * ad;
    fact *= i * fn;
    auto pn = _Tp{1} - abcd * (a + ix) / fn;
    if (n == 1)
      return {pn, fact};

    fn = (ac + _Tp{1}) * (ad + _Tp{1});
    fact *= i * fn / _Tp{2};
    auto An = -abcd * fn / (abcd + _Tp{1}) / (abcd + _Tp{2});
    auto Cn = bc * db / abcd / (abcd + _Tp{1});

    auto pnp1 = ((An + Cn + (a + ix)) * pn - Cn * pnm1) / An;

    for (int k = 2; k < n; ++k)
      {
	const auto abcd2n = abcd + _Tp(2 * k);

	fn = (ac + _Tp(k)) * (ad + _Tp(k));
	fact *= i * fn / _Tp(1 + k);

	auto An = -(abcd + _Tp(k - 1)) * fn / (abcd2n - _Tp(1)) / abcd2n;
	auto Cn = _Tp(k)
		* (bc + _Tp(k - 1)) * (db + _Tp(k - 1))
		/ (abcd2n - _Tp(2)) / (abcd2n - _Tp(1));

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
template<typename _Tp, typename _TpX>
  __continuous_hahn_t<_Tp>
  __continuous_hahn(int n, _Tp a, _Tp b, _Tp c, _Tp d, _TpX x)
  {
    if (std::isnan(a))
      return {std::complex<_Tp>(a), std::complex<_Tp>{}};
    else if (std::isnan(b))
      return {std::complex<_Tp>(b), std::complex<_Tp>{}};
    else if (std::isnan(c))
      return {std::complex<_Tp>(c), std::complex<_Tp>{}};
    else if (std::isnan(d))
      return {std::complex<_Tp>(d), std::complex<_Tp>{}};
    else if (std::isnan(x))
      return {std::complex<_Tp>(x), std::complex<_Tp>{}};
    else
      return __continuous_hahn_recur(n, a, b, c, d, x);
  }

template<typename _Tp>
  void
  test_continuous_hahn(int n_max, _Tp a, _Tp b, _Tp c, _Tp d)
  {
    std::cout.precision(std::numeric_limits<_Tp>::digits10);
    auto w = std::cout.precision() + 8;

    for (int n = 0; n <= n_max; ++n)
      {
	std::cout << '\n' << '\n' << " n = " << n << '\n';
	for (int i = -200; i <= 200; ++i)
	  {
	    auto x = i * _Tp{0.01L};
	    auto p = __continuous_hahn(n, a, b, c, d, x);
	    std::cout << ' ' << std::setw(w) << x
		      << ' ' << std::setw(w) << std::real(p.__factor * p.__value)
		      << ' ' << std::setw(w) << std::imag(p.__factor * p.__value)
		      << '\n';
	  }
      }
  }

int
main()
{
  test_continuous_hahn(10, 1.0f, 2.0f, 2.0f, 1.0f);
}
