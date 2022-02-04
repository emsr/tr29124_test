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
  continuous_dual_hahn_t
  {
    Tp value;
    Tp factor;
  };

/**
 * Compute the continuous dual Hahn polynomial by recursion:
 * @f[
 *    -(a^2 + x^2) S_n(x^2) = A_n S_{n+1}(x^2)
 *                            - (A_n + C_n) S_n(x^2)
 *                            + C_n S_{n-1}(x^2)
 * @f]
 * where @f$ S_n(x^2) = S_n(x^2; a, b, c) @f$ and
 * @f[
 *    A_n = (n + a + b)(n + a + c)
 * @f]
 * and
 * @f[
 *    C_n = n(n + b + c - 1)
 * @f]
 */
template<typename Tp, typename _TpX>
  continuous_dual_hahn_t<Tp>
  continuous_dual_hahn_recur(int n, Tp a, Tp b, Tp c, _TpX x)
  {
    auto Snm1 = Tp{1};
    if (n == 0)
      return {Snm1, Tp{1}};

    const auto aa = a * a;
    const auto ab = a + b;
    const auto ac = a + c;
    const auto bc = b + c;
    const auto xx = Tp(x * x);

    auto fact = Tp{1};

    auto fn = ab * ac;
    fact *= fn;
    auto Sn = Tp{1} - (aa + xx) / fn;
    if (n == 1)
      return {Sn, fact};

    fn = (ab + Tp{1}) * (ac + Tp{1});
    fact *= fn;
    auto An = fn;
    auto Cn = bc;

    auto Snp1 = ((An + Cn - (aa + xx)) * Sn - Cn * Snm1) / An;

    for (int k = 2; k < n; ++k)
      {
	fn = (ab + Tp(k)) * (ac + Tp(k));
	fact *= fn;
	auto An = fn;
	auto Cn = Tp(k) * (bc + Tp(k - 1));

	Snm1 = Sn;
	Sn = Snp1;
        Snp1 = ((An + Cn - (aa + xx)) * Sn - Cn * Snm1) / An;
      }

    return {Snp1, fact};
  }

/**
 * Compute the continuous dual Hahn polynomial defined by
 * @f[
 *    S_n() = \frac{1}{(a + b)_n(a + c)_n}
 *            {}_3F_2(-n, a + ix, a - ix, a + b, a + c; 1)
 * @f]
 */
template<typename Tp, typename _TpX>
  continuous_dual_hahn_t<Tp>
  continuous_dual_hahn(int n, Tp a, Tp b, Tp c, _TpX x)
  {
    if (std::isnan(a))
      return {a, Tp{}};
    else if (std::isnan(b))
      return {b, Tp{}};
    else if (std::isnan(c))
      return {c, Tp{}};
    else if (std::isnan(x))
      return {x, Tp{}};
    else
      return continuous_dual_hahn_recur(n, a, b, c, x);
  }

/**
 * 
 */
template<typename Tp>
  void
  test_continuous_dual_hahn(int n_max, Tp a, Tp b, Tp c)
  {
    std::cout.precision(std::numeric_limits<Tp>::digits10);
    auto w = std::cout.precision() + 8;

    for (int n = 0; n <= n_max; ++n)
      {
	std::cout << '\n' << '\n' << " n = " << n << '\n';
	for (int i = 0; i <= 220; ++i)
	  {
	    auto x = i * Tp{0.05L};
	    auto S = continuous_dual_hahn(n, a, b, c, x);
	    std::cout << ' ' << std::setw(w) << x
		      << ' ' << std::setw(w) << S.value
		      << '\n';
	  }
      }
  }

int
main()
{
  test_continuous_dual_hahn(10, 2.0f, 2.0f, 2.0f);
}
