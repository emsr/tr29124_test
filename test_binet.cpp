/*
 $HOME/bin_tr29124/bin/g++ -g -D__STDCPP_WANT_MATH_SPEC_FUNCS__ -I. -o test_binet test_binet.cpp

 ./test_binet > test_binet.txt

 $HOME/bin/bin/g++ -g -D__STDCPP_WANT_MATH_SPEC_FUNCS__ -I. -o test_binet test_binet.cpp

 */

#include <limits>
#include <iostream>
#include <iomanip>
#include <vector>
#include "rational.h"
#include <cmath>

namespace std
{
namespace __detail
{

  /**
   * Computes the sequence of Bernoulli numbers @f$B_{2m}@f$ (m > 0)
   * using the Akiyama-Tanigawa algorithm.
   * This might be unstable.
   */
  template<typename _Rat>
    std::vector<_Rat>
    __bernoulli_a_t(std::size_t __len)
    {
      auto __n = 2 * __len + 1;
      std::vector<_Rat> __t;
      std::vector<_Rat> __a;

      __t.emplace_back(1LL);

      for (std::size_t __m = 1; __m < __n; ++__m)
	{
	  __t.push_back(_Rat(1, __m + 1));
	  for (std::size_t __j = __m; __j > 0; --__j)
            __t[__j - 1] = _Rat(__j) * (__t[__j - 1] - __t[__j]);

	  // Get all Bernoulli numbers by deleting the 'if' clause.
	  if ((__m & 1) == 0)
	    __a.push_back(__t[0]);
	}

      return __a;
    }

  /**
   * Scales the even Bernoulli numbers @f$ B_{2m} @f$ with weights
   * @f$ (-1)^m/((2m - 1)2m) @f$, m > 0.
   */
  template<typename _Rat>
    std::vector<_Rat>
    __weights(std::vector<_Rat> __b)
    {
      int __sgn = 1;
      for (std::size_t __m = 0; __m < __b.size(); ++__m)
	{
	  __b[__m] *= _Rat(__sgn, (2 * __m + 1) * (2 * __m + 2));
	  __sgn = -__sgn;
	}

      return __b;
    }

  // Computes Rutishauser's Quotient-Difference (QD) algorithm
  template<typename _Rat>
    std::vector<_Rat>
    __quotient_difference(std::vector<_Rat> __s)
    {
      auto __len = __s.size();
      auto __zero = _Rat(0);
      std::vector<std::vector<_Rat>> __m;
      std::vector<_Rat> __r;

      for (std::size_t __n = 0; __n < __len; ++__n)
	{
	  __m.push_back(std::vector<_Rat>());
	  __m.back().push_back(__zero);
	}
      for (std::size_t __n = 0; __n < __len - 1; ++__n)
	__m[__n].push_back(__s[__n + 1] / __s[__n]);

      __r.push_back(__s[0]);
      __r.push_back(__m[0][1]);

      for (std::size_t __k = 2; __k < __len; ++__k)
	{
	  for (std::size_t __n = 0; __n < __len - __k ; ++__n)
            {
              auto __a = __m[__n + 1][__k - 2];
              auto __b = __m[__n + 1][__k - 1];
              auto __c = __m[__n][__k - 1];
              __m[__n].push_back((__k & 1) == 0
				 ? __a + __b - __c
				 : __a * __b / __c);
            }
	  __r.push_back(__m[0][__k]);
	}

      return __r;
    }

  // Computes the Stieltjes continued fraction for the
  // Gamma function using Rutishauser's QD-algorithm.
  template<typename _Rat>
    std::vector<_Rat>
    __stieltjes_cont_frac_seq(int __len)
    { return __quotient_difference(__weights(__bernoulli_a_t<_Rat>(__len))); }

  /**
   *
   */
  template<typename _Tp>
    _Tp
    __binet_asymp(_Tp __z)
    {
      using _Val = _Tp;
      using _Real = std::__detail::__num_traits_t<_Val>;
      constexpr auto _S_eps = std::numeric_limits<_Real>::epsilon();

      // Weighted Bernoulli numbers: (-1)^k B_2k / ((2*k+1)*(2*k+2))
      constexpr std::size_t _S_n = 17;
      constexpr _Real
      _S_b[_S_n]
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
	_Real{236364091LL}     / _Real{1506960LL},
	_Real{657931LL}        / _Real{300LL},
	_Real{3392780147LL}    / _Real{93960LL},
	_Real{1723168255201LL} / _Real{2492028LL},
	_Real{7709321041217LL} / _Real{505920LL},
	_Real{151628697551LL}  / _Real{396LL}
      };

      auto __z2 = _Real{1} / (__z * __z);
      auto __J = _Val{};
      auto __zk = _Val{1};
      for (auto __b : _S_b)
	{
	  auto __term = __b * __zk;
	  __J += __term;
	  if (std::abs(__term) < _S_eps * std::abs(__J))
	    break;
	  __zk *= __z2;
	}
      __J /= __z;

      return __J;
    }

  /**
   *
   */
  template<typename _Tp>
    _Tp
    __binet_cont_frac(_Tp __z)
    {
      using _Val = _Tp;
      using _Real = std::__detail::__num_traits_t<_Val>;

      // Stieltjes partial numerators.
      constexpr std::size_t _S_n = 7;
      constexpr _Real
      _S_a[_S_n]
      {
	_Real{1LL}            / _Real{12LL},
	_Real{1LL}            / _Real{30LL},
	_Real{53LL}           / _Real{210LL},
	_Real{195LL}          / _Real{371LL},
	_Real{22999LL}        / _Real{22737LL},
	_Real{29944523LL}     / _Real{19733142LL},
	_Real{109535241009LL} / _Real{48264275462LL}
      };

      // Backward recurrence.
      auto __w = _Val{}; // The tail function.
      auto __J = __w;
      for (std::ptrdiff_t __k = _S_n - 1; __k >= 0; --__k)
        __J = _S_a[__k] / (__z + __J);

      return __J;
    }

  /**
   *
   */
  template<typename _Tp>
    _Tp
    __binet(_Tp __z)
    {
      using _Val = _Tp;
      using _Real = std::__detail::__num_traits_t<_Val>;

      constexpr auto _S_switchover = _Real{10}; /// @todo Find Binet function switch.

      if (std::__detail::__isnan(__z))
	return __gnu_cxx::__quiet_NaN<_Tp>();
      else if (__z < _S_switchover)
	return __binet_cont_frac(__z);
      else
	return __binet_asymp(__z);
    }

} // namespace __detail
} // namespace std


namespace __gnu_cxx
{

  /**
   *
   */
  template<typename _Tp>
    _Tp
    lgamma_scaled(_Tp __z)
    { return std::__detail::__binet(__z); }

  /**
   *
   */
  template<typename _Tp>
    _Tp
    tgamma_scaled(_Tp __z)
    { return std::exp(std::__detail::__binet(__z)); }

} // namespace__gnu_cxx


template<typename _Tp>
  void
  test()
  {
    using _Rat = __gnu_cxx::_Rational<_Tp>;

    using _Real = long double;

    std::cout.precision(std::numeric_limits<_Real>::digits10);
    auto width = std::cout.precision() + 6;

    std::cout << "\nBernoulli numbers\n";
    auto bern = std::__detail::__bernoulli_a_t<_Rat>(20);
    for (auto& b : bern)
      std::cout << b << '\n';

    std::cout << "\nWeighted Bernoulli numbers\n";
    auto wts = std::__detail::__weights(bern);
    for (auto& w : wts)
      std::cout << w << '\n';

    std::cout << "\nStieltjes partial numerators\n";
    int len = 7;
    auto cf = std::__detail::__stieltjes_cont_frac_seq<_Rat>(len);
    for (int k = 0; k < len; ++k)
      std::cout << k + 1 << ": " << cf[k]
		<< " = " << __gnu_cxx::_Rational_cast<_Real>(cf[k]) << '\n';

    std::cout << "\nBinet asymptotic\n";
    for (int k = 1; k <= 5000; ++k)
      {
	auto x = 0.1L * k;
	auto j_as = std::__detail::__binet_asymp(x);
	auto j_cf = std::__detail::__binet_cont_frac(x);
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
