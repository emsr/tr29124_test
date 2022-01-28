/**
 *
 */

#include <limits>
#include <iostream>
#include <iomanip>
#include <vector>
#include <cmath>

#include <emsr/rational.h>
#include <emsr/fp_type_util.h>
#include <emsr/numeric_limits.h>

namespace emsr
{
namespace detail
{

  /**
   * Computes the sequence of Bernoulli numbers @f$B_{2m}@f$ (m > 0)
   * using the Akiyama-Tanigawa algorithm.
   * This might be unstable.
   */
  template<typename _Rat>
    std::vector<_Rat>
    bernoulli_a_t(std::size_t len)
    {
      auto n = 2 * len + 1;
      std::vector<_Rat> t;
      std::vector<_Rat> a;

      t.emplace_back(1LL);

      for (std::size_t m = 1; m < n; ++m)
	{
	  t.push_back(_Rat(1, m + 1));
	  for (std::size_t j = m; j > 0; --j)
            t[j - 1] = _Rat(j) * (t[j - 1] - t[j]);

	  // Get all Bernoulli numbers by deleting the 'if' clause.
	  if ((m & 1) == 0)
	    a.push_back(t[0]);
	}

      return a;
    }

  /**
   * Scales the even Bernoulli numbers @f$ B_{2m} @f$ with weights
   * @f$ (-1)^m/((2m - 1)2m) @f$, m > 0.
   */
  template<typename _Rat>
    std::vector<_Rat>
    weights(std::vector<_Rat> b)
    {
      int sgn = 1;
      for (std::size_t m = 0; m < b.size(); ++m)
	{
	  b[m] *= _Rat(sgn, (2 * m + 1) * (2 * m + 2));
	  sgn = -sgn;
	}

      return b;
    }

  // Computes Rutishauser's Quotient-Difference (QD) algorithm
  template<typename _Rat>
    std::vector<_Rat>
    quotient_difference(std::vector<_Rat> s)
    {
      auto len = s.size();
      auto zero = _Rat(0);
      std::vector<std::vector<_Rat>> m;
      std::vector<_Rat> r;

      for (std::size_t n = 0; n < len; ++n)
	{
	  m.push_back(std::vector<_Rat>());
	  m.back().push_back(zero);
	}
      for (std::size_t n = 0; n < len - 1; ++n)
	m[n].push_back(s[n + 1] / s[n]);

      r.push_back(s[0]);
      r.push_back(m[0][1]);

      for (std::size_t k = 2; k < len; ++k)
	{
	  for (std::size_t n = 0; n < len - k ; ++n)
            {
              auto a = m[n + 1][k - 2];
              auto b = m[n + 1][k - 1];
              auto c = m[n][k - 1];
              m[n].push_back((k & 1) == 0
				 ? a + b - c
				 : a * b / c);
            }
	  r.push_back(m[0][k]);
	}

      return r;
    }

  // Computes the Stieltjes continued fraction for the
  // Gamma function using Rutishauser's QD-algorithm.
  template<typename _Rat>
    std::vector<_Rat>
    stieltjes_cont_frac_seq(int len)
    { return quotient_difference(weights(bernoulli_a_t<_Rat>(len))); }

  /**
   *
   */
  template<typename Tp>
    Tp
    binet_asymp(Tp z)
    {
      using Val = Tp;
      using Real = emsr::num_traits_t<Val>;
      constexpr auto s_eps = std::numeric_limits<Real>::epsilon();

      // Weighted Bernoulli numbers: (-1)^k B_2k / ((2*k+1)*(2*k+2))
      constexpr std::size_t s_n = 12;
      constexpr Real
      s_b[s_n]
      {
	Real{1LL}             / Real{12LL},
	Real{1LL}             / Real{360LL},
	Real{1LL}             / Real{1260LL},
	Real{1LL}             / Real{1680LL},
	Real{1LL}             / Real{1188LL},
	Real{691LL}           / Real{360360LL},
	Real{1LL}             / Real{156LL},
	Real{3617LL}          / Real{122400LL},
	Real{43867LL}         / Real{244188LL},
	Real{174611LL}        / Real{125400LL},
	Real{77683LL}         / Real{5796LL},
	Real{236364091LL}     / Real{1506960LL}
      };

      auto z2 = Real{1} / (z * z);
      auto J = Val{};
      auto zk = Val{1};
      for (auto b : s_b)
	{
	  auto term = b * zk;
	  J += term;
	  if (std::abs(term) < s_eps * std::abs(J))
	    break;
	  zk *= z2;
	}
      J /= z;

      return J;
    }

  /**
   *
   */
  template<typename Tp>
    Tp
    binet_cont_frac(Tp z)
    {
      using Val = Tp;
      using Real = emsr::num_traits_t<Val>;

      // Stieltjes partial numerators.
      constexpr std::size_t s_n = 6;
      constexpr Real
      s_a[s_n]
      {
	Real{1LL}            / Real{12LL},
	Real{1LL}            / Real{30LL},
	Real{53LL}           / Real{210LL},
	Real{195LL}          / Real{371LL},
	Real{22999LL}        / Real{22737LL},
	Real{29944523LL}     / Real{19733142LL}
      };

      // Backward recurrence.
      auto w = Val{}; // The tail function.
      auto J = w;
      for (std::ptrdiff_t k = s_n - 1; k >= 0; --k)
        J = s_a[k] / (z + J);

      return J;
    }

  /**
   *
   */
  template<typename Tp>
    Tp
    binet(Tp z)
    {
      using Val = Tp;
      using Real = emsr::num_traits_t<Val>;

      constexpr auto s_switchover = Real{10}; /// @todo Find Binet function switch.

      if (std::isnan(z))
	return emsr::quiet_NaN<Tp>();
      else if (z < s_switchover)
	return binet_cont_frac(z);
      else
	return binet_asymp(z);
    }

} // namespace detail
} // namespace emsr


namespace emsr
{

  /**
   *
   */
  template<typename Tp>
    Tp
    lgamma_scaled(Tp z)
    { return emsr::detail::binet(z); }

  /**
   *
   */
  template<typename Tp>
    Tp
    tgamma_scaled(Tp z)
    { return std::exp(emsr::detail::binet(z)); }

} // namespace emsr


template<typename Tp>
  void
  test()
  {
    using _Rat = emsr::Rational<Tp>;

    using Real = long double;

    std::cout.precision(std::numeric_limits<Real>::digits10);
    auto width = std::cout.precision() + 6;

    std::cout << "\nBernoulli numbers\n";
    auto bern = emsr::detail::bernoulli_a_t<_Rat>(20);
    for (auto& b : bern)
      std::cout << b << '\n';

    std::cout << "\nWeighted Bernoulli numbers\n";
    auto wts = emsr::detail::weights(bern);
    for (auto& w : wts)
      std::cout << w << '\n';

    std::cout << "\nStieltjes partial numerators\n";
    int len = 7;
    auto cf = emsr::detail::stieltjes_cont_frac_seq<_Rat>(len);
    for (int k = 0; k < len; ++k)
      std::cout << k + 1 << ": " << cf[k]
		<< " = " << emsr::Rational_cast<Real>(cf[k]) << '\n';

    std::cout << "\nBinet asymptotic\n";
    for (int k = 1; k <= 5000; ++k)
      {
	auto x = 0.1L * k;
	auto j_as = emsr::detail::binet_asymp(x);
	auto j_cf = emsr::detail::binet_cont_frac(x);
	std::cout << ' ' << std::setw(width) << x
		  << ' ' << std::setw(width) << j_as
		  << ' ' << std::setw(width) << j_cf
		  << ' ' << std::setw(width) << (j_as - j_cf) / j_cf << '\n';
      }
  }

int
main()
{
  test<long long>();
}
