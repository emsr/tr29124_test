/**
 *
 */

#include <limits>
#include <iostream>
#include <iomanip>
#include <vector>
#include <cmath>

#include <emsr/rational.h>

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
  template<typename _Tp>
    _Tp
    binet_asymp(_Tp z)
    {
      using _Val = _Tp;
      using _Real = emsr::num_traits_t<_Val>;
      constexpr auto s_eps = std::numeric_limits<_Real>::epsilon();

      // Weighted Bernoulli numbers: (-1)^k B_2k / ((2*k+1)*(2*k+2))
      constexpr std::size_t s_n = 12;
      constexpr _Real
      s_b[s_n]
      {
	_Real{1LL}             / _Real{12LL},
	_Real{1LL}             / _Real{360LL},
	_Real{1LL}             / _Real{1260LL},
	_Real{1LL}             / _Real{1680LL},
	_Real{1LL}             / _Real{1188LL},
	_Real{691LL}           / _Real{360360LL},
	_Real{1LL}             / _Real{156LL},
	_Real{3617LL}          / _Real{122400LL},
	_Real{43867LL}         / _Real{244188LL},
	_Real{174611LL}        / _Real{125400LL},
	_Real{77683LL}         / _Real{5796LL},
	_Real{236364091LL}     / _Real{1506960LL}
      };

      auto z2 = _Real{1} / (z * z);
      auto J = _Val{};
      auto zk = _Val{1};
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
  template<typename _Tp>
    _Tp
    binet_cont_frac(_Tp z)
    {
      using _Val = _Tp;
      using _Real = emsr::num_traits_t<_Val>;

      // Stieltjes partial numerators.
      constexpr std::size_t s_n = 6;
      constexpr _Real
      s_a[s_n]
      {
	_Real{1LL}            / _Real{12LL},
	_Real{1LL}            / _Real{30LL},
	_Real{53LL}           / _Real{210LL},
	_Real{195LL}          / _Real{371LL},
	_Real{22999LL}        / _Real{22737LL},
	_Real{29944523LL}     / _Real{19733142LL}
      };

      // Backward recurrence.
      auto w = _Val{}; // The tail function.
      auto J = w;
      for (std::ptrdiff_t k = s_n - 1; k >= 0; --k)
        J = s_a[k] / (z + J);

      return J;
    }

  /**
   *
   */
  template<typename _Tp>
    _Tp
    binet(_Tp z)
    {
      using _Val = _Tp;
      using _Real = emsr::num_traits_t<_Val>;

      constexpr auto s_switchover = _Real{10}; /// @todo Find Binet function switch.

      if (std::isnan(z))
	return emsr::quiet_NaN<_Tp>();
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
  template<typename _Tp>
    _Tp
    lgamma_scaled(_Tp z)
    { return emsr::detail::binet(z); }

  /**
   *
   */
  template<typename _Tp>
    _Tp
    tgamma_scaled(_Tp z)
    { return std::exp(emsr::detail::binet(z)); }

} // namespace emsr


template<typename _Tp>
  void
  test()
  {
    using _Rat = emsr::Rational<_Tp>;

    using _Real = long double;

    std::cout.precision(std::numeric_limits<_Real>::digits10);
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
		<< " = " << emsr::Rational_cast<_Real>(cf[k]) << '\n';

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
