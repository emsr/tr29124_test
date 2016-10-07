/*
 $HOME/bin_tr29124/bin/g++ -g -D__STDCPP_WANT_MATH_SPEC_FUNCS__ -I. -o test_binet_float test_binet_float.cpp -lquadmath

 ./test_binet_float > test_binet_float.txt

 $HOME/bin/bin/g++ -g -std=gnu++17 -I. -o test_binet_float test_binet_float.cpp -lquadmath

 */

#include <limits>
#include <iostream>
#include <iomanip>
#include <vector>
#include <bits/float128.h>
#include <bits/complex128.h>
#include <cmath>
#include <ext/polynomial.h>

namespace std
{
namespace __detail
{

  /**
   * Computes the sequence of Bernoulli numbers @f$B_{2m}@f$ (m > 0)
   * using the Akiyama-Tanigawa algorithm.
   * This might be unstable.
   */
  template<typename _Real>
    std::vector<_Real>
    __bernoulli_a_t(std::size_t __len)
    {
      auto __n = 2 * __len + 1;
      std::vector<_Real> __t;
      std::vector<_Real> __a;

      __t.emplace_back(1LL);

      for (std::size_t __m = 1; __m < __n; ++__m)
	{
	  __t.push_back(_Real{1} / _Real(__m + 1));
	  for (int __j = __m; __j > 0; --__j)
            __t[__j - 1] = _Real(__j) * (__t[__j - 1] - __t[__j]);

	  // Get all Bernoulli numbers by deleting the 'if' clause.
	  if ((__m & 1) == 0)
	    __a.push_back(__t[0]);
	}

      return __a;
    }

  /**
   * Computes the sequence of even Bernoulli numbers @f$ B_{2m} @f$ (m > 0)
   * using the old Riemann zeta function series.
   */
  template<typename _Real>
    std::vector<_Real>
    __bernoulli_vec(std::size_t __len)
    {
      std::vector<_Real> __a;
      for (std::size_t __m = 1; __m <= __len; ++__m)
	__a.push_back(__bernoulli<_Real>(2 * __m));
      return __a;
    }

  /**
   * Scales the even Bernoulli numbers @f$ B_{2m} @f$ with weights
   * @f$ (-1)^m/((2m - 1)2m) @f$, m > 0.
   */
  template<typename _Real>
    std::vector<_Real>
    __weights(std::vector<_Real> __b)
    {
      int __sgn = 1;
      for (std::size_t __m = 0; __m < __b.size(); ++__m)
	{
	  __b[__m] *= _Real(__sgn) / _Real(2 * __m + 1) / _Real(2 * __m + 2);
	  __sgn = -__sgn;
	}

      return __b;
    }

  /**
   * Computes the partial numerators for the Binet function
   * using Rutishauser's Quotient-Difference (QD) algorithm.
   */
  template<typename _Real>
    std::vector<_Real>
    __quotient_difference(std::vector<_Real> __s)
    {
      auto __len = __s.size();
      auto __zero = _Real{0};
      std::vector<std::vector<_Real>> __m;
      std::vector<_Real> __r;

      for (std::size_t __n = 0; __n < __len; ++__n)
	{
	  __m.push_back(std::vector<_Real>());
	  __m.back().push_back(__zero); // e[k+1,0] = 0, k >= 0
	}
      for (std::size_t __n = 0; __n < __len - 1; ++__n)
	__m[__n].push_back(__s[__n + 1] / __s[__n]); // q[k,1] = c[k+1]/c[k], k >= 0

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

  /**
   * For a series specified by coefficients @f$ c_k @f$:
   * @f[
   *  \Lambda(x) = \sum_{k=0}^{N} c_k x^k
   * @f]
   * form the coefficients @f$ d_k @f$ for the inverse
   * @f[
   *  \frac{1}{\Lambda(x)} = \sum_{k=0}^{N} d_k x^k
   * @f]
   */
  template<typename _Real>
    std::vector<_Real>
    __inverse_series(std::vector<_Real> __c)
    {
      if (__c.size() == 0)
	return std::vector<_Real>{};
      else if (__c[0] == _Real{0})
	std::__throw_domain_error("__inverse_series: "
				  "first (constant) coefficient is zero.");
      else if (__c.size() == 1)
	return std::vector<_Real>{{_Real{1} / __c[0]}};
      else
	{
	  auto __n = __c.size();
	  auto __m = __n + 1;
	  std::vector<_Real> __lambda(__n);
	  __lambda[0] = _Real{0};
	  for (unsigned __i = 1; __i < __n; ++__i)
	    __lambda[__i] = __c[__i] / __c[0];

	  std::vector<_Real> __lambdak(__m);
	  std::vector<_Real> __d(__m);
	  __d[0] = __lambdak[0] = _Real{1}; // k == 0
	  for (unsigned __k = 1; __k < __m; ++__k)
	    {
	      std::vector<_Real> __work(__m);
	      for (unsigned __i = 1; __i < __n; ++__i)
		for (unsigned __j = __k - 1; __j < __m; ++__j)
		  if (__i + __j < __m)
		    __work[__i + __j] += __lambda[__i] * __lambdak[__j];
	      std::swap(__work, __lambdak);
	      auto __sign = (__k % 2 == 0 ? +1 : -1);
	      for (unsigned __j = __k; __j < __m; ++__j)
		__d[__j] += __sign * __lambdak[__j];
	    }
	  for (unsigned __j = 0; __j < __m; ++__j)
	    __d[__j] /= __c[0];
	  return __d;
	}
    }

  /**
   * Computes the partial numerators for the Binet function
   * using Rutishauser's Quotient-Difference (QD) algorithm.
   * Use the more numerically stable progressive version.
   */
  template<typename _Real>
    std::vector<_Real>
    __quotient_difference_prog(std::vector<_Real> __s)
    {
      auto __n = __s.size();
      auto __zero = _Real{0};
      std::vector<_Real> __r;

      std::vector<std::vector<_Real>> __q;
      __q.push_back(std::vector<_Real>{});
      __q.back().push_back(-__s[1] / __s[0]);
      for (unsigned __k = 1; __k < __n; ++__k)
	__q.back().push_back(__zero);

      std::vector<std::vector<_Real>> __e;
      __e.push_back(std::vector<_Real>{});
      __e.back().push_back(_Real{0});
      for (unsigned __k = 1; __k < __n - 2; ++__k)
	__e.back().push_back(__s[__k + 1] / __s[__k]);

      __r.push_back(__e[0][0]);
      __r.push_back(__q[0][1]);

      for (unsigned __k = 1; __k < __n - 2; ++__k)
	{
	  __q.push_back(std::vector<_Real>{});
	  for (unsigned __l = 1; __l < __n - 2; ++__l)
	    __q.back().push_back(__q[__k - 1][__l] + __e[__k - 1][__l] - __e[__k - 1][__l - 1]);

	  __e.push_back(std::vector<_Real>{});
	  for (unsigned __l = 1; __l < __n - 2; ++__l)
	    __e.back().push_back(__q[__k - 1][__l + 1] * __e[__k - 1][__l] / __q.back().back());

	  __r.push_back(__q[0][__k]);
	  __r.push_back(__e[0][__k]);
	}

      return __r;
    }

  // Computes the Stieltjes continued fraction for the
  // Gamma function using Rutishauser's QD-algorithm.
  template<typename _Real>
    std::vector<_Real>
    __stieltjes_cont_frac_seq(int __len)
    { return __quotient_difference(__weights(__bernoulli_vec<_Real>(__len))); }

  // Computes the Stieltjes continued fraction for the
  // Gamma function using the progressive QD-algorithm.
  template<typename _Real>
    std::vector<_Real>
    __stieltjes_cont_frac_seq_prog(int __len)
    { return __quotient_difference_prog(__inverse_series(__weights(__bernoulli_vec<_Real>(__len)))); }

  /**
   * Compute the Binet function using the asymptotic series
   */
  template<typename _Tp>
    _Tp
    __binet_asymp(_Tp __z)
    {
      using _Val = _Tp;
      using _Real = std::__detail::__num_traits_t<_Val>;
      constexpr auto _S_eps = std::numeric_limits<_Real>::epsilon();

      // Weighted Bernoulli numbers: (-1)^k B_{2k} / ((2k + 1)(2k + 2)), k > 0.
      constexpr std::size_t _S_n = 50;
      constexpr _Real
      _S_b[_S_n]
      {
	0.0833333333333333333L,
	0.00277777777777777778L,
	0.000793650793650793651L,
	0.000595238095238095238L,
	0.000841750841750841751L,
	0.00191752691752691753L,
	0.00641025641025641026L,
	0.0295506535947712418L,
	0.179644372368830573L,
	1.39243221690590112L,
	13.402864044168392L,
	156.848284626002017L,
	2193.10333333333333L,
	36108.7712537249894L,
	691472.268851313067L,
	15238221.5394074162L,
	382900751.391414141L,
	10882266035.7843911L,
	347320283765.002252L,
	12369602142269.2744L,
	488788064793079.335L,
	21320333960919373.9L,
	1.02177529652570008e+18L,
	5.35754721733002036e+19L,
	3.06157826370488341e+21L,
	1.8999917426399204e+23L,
	1.27633740338288341e+25L,
	9.2528471761204163e+26L,
	7.21882259518561029e+28L,
	6.04518340599585696e+30L,
	5.42067047157009454e+32L,
	5.19295781531408194e+34L,
	5.30365885511970059e+36L,
	5.76332534816496401e+38L,
	6.65115571484845393e+40L,
	8.13737835813668053e+42L,
	1.05369669533571418e+45L,
	1.44181805999622062e+47L,
	2.08173565220895654e+49L,
	3.16702266348866618e+51L,
	5.07000646121113734e+53L,
	8.52997282030055187e+55L,
	1.50641728093405986e+58L,
	2.78934947038316368e+60L,
	5.40935043528604149e+62L,
	1.09753378215085198e+65L,
	2.32748762026184791e+67L,
	5.15392916206532138e+69L,
	1.19062102308902264e+72L,
	2.86689389602966737e+74L,
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
      constexpr std::size_t _S_n = 50;
      constexpr _Real
      _S_a[_S_n]
      {
	0.0833333333333333333L,
	0.0333333333333333333L,
	0.252380952380952381L,
	0.525606469002695418L,
	1.01152306812684171L,
	1.5174736491532874L,
	2.26948897420495996L,
	3.00991738325939818L,
	4.02688719234390119L,
	5.00276808075403014L,
	6.28391137081578202L,
	7.49591912238403412L,
	9.04066023436772668L,
	10.4893036545094816L,
	12.2971936103862081L,
	13.9828769539924246L,
	16.0535514167049548L,
	17.9766073998701922L,
	20.3097620274420059L,
	22.4704716399319049L,
	25.0658465489497183L,
	27.4644518250188972L,
	30.3218212316992515L,
	32.9585339299094923L,
	36.0776989314482468L,
	38.9527066819767291L,
	42.3334900442987634L,
	45.4469608486042283L,
	49.0892031317448499L,
 	52.4412887468308445L,
	56.3448453518945826L,
	59.9356839007347235L,
	64.1004227533901L,
	67.9301408283552804L,
	72.3559403887712393L,
	76.4246551820146945L,
	81.1114015206437931L,
	85.4192263186344977L,
	90.3668004955012448L,
	94.9138753967955118L,
	100.122079613541441L,
	104.908742524024746L,
	110.376919022080575L,
	115.404525862345152L,
	121.129836041189223L,
	126.404299775283082L,
	132.374557373786024L,
	137.920563636156107L,
	144.086972602240632L,
	149.997077690266067L,
      };

      // Backward recurrence.
      auto __w = _Val{}; // The tail function.
      auto __J = __w;
      for (std::ptrdiff_t __k = _S_n - 1; __k >= 0; --__k)
        __J = _S_a[__k] / (__z + __J);

      return __J;
    }

  /**
   * Compute the Binet function @f$ J(z) @f$ defined by
   * @f[
   *    J(z) = log\left(\Gamma(z)\right) + z
   *         - \left(z-\frac{1}{2}\right) log(z) - log(2\pi)
   * @f]
   * or
   * @f[
   *    \Gamma(z) = \sqrt{2\pi}z^{z-\frac{1}{2}}e^{-z}e^{J(z)}
   * @f]
   * where @f$ \Gamma(z) @f$ is the gamma function.
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
   * Return the Binet function @f$ J(z) @f$ or the scaled log gamma function
   * @f$ \Gamma^*(z) @f$ for @c float argument @f$ z @f$.
   *
   * @see lgamma_scaled for details.
  float
  lgamma_scaledf(float __z)
  { return std::__detail::__binet<float>(__z); }
   */

  /**
   * Return the Binet function @f$ J(z) @f$ or the scaled log gamma function
   * @f$ \Gamma^*(z) @f$ for <tt>long double</tt> argument @f$ z @f$.
   *
   * @see lgamma_scaled for details.
   */
  long double
  lgamma_scaledl(long double __z)
  { return std::__detail::__binet<long double>(__z); }

  /**
   * Return the Binet function @f$ J(z) @f$ or the scaled log gamma function
   * @f$ log(\Gamma^*(z)) @f$ defined by
   * @f[
   *    J(z) = log(\Gamma^*(z)) = log\left(\Gamma(z)\right) + z
   *         - \left(z-\frac{1}{2}\right) log(z) - log(2\pi)
   * @f]
   * or
   * @f[
   *    \Gamma(z) = \sqrt{2\pi}z^{z-\frac{1}{2}}e^{-z}e^{J(z)}
   * @f]
   * where @f$ \Gamma(z) @f$ is the gamma function.
   */
  template<typename _Tp>
    _Tp
    lgamma_scaled(_Tp __z)
    { return std::__detail::__binet(__z); }


  /**
   * Return the Binet function or the scaled log gamma
   * @f$ \Gamma^*(z) @f$ for @c float argument @f$ z @f$.
   *
   * @see tgamma_scaled for details.
  float
  tgamma_scaledf(float __z)
  { return std::exp(std::__detail::__binet<float>(__z)); }
   */

  /**
   * Return the Binet function or the scaled log gamma
   * @f$ \Gamma^*(z) @f$ for <tt>long double</tt> argument @f$ z @f$.
   *
   * @see tgamma_scaled for details.
   */
  long double
  tgamma_scaledl(long double __z)
  { return std::exp(std::__detail::__binet<long double>(__z)); }

  /**
   * Return the scaled log gamma @f$ \Gamma^*(z) @f$ defined by
   * @f[
   *    \Gamma^*(z) = \Gamma(z)/\sqrt{2\pi}z^{z-\frac{1}{2}}e^{-z}
   * @f]
   * or
   * @f[
   *    \Gamma(z) = \sqrt{2\pi}z^{z-\frac{1}{2}}e^{-z}\Gamma^*(z)
   * @f]
   * where @f$ \Gamma(z) @f$ is the gamma function.
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
    using _Real = _Tp;

    std::cout.precision(std::numeric_limits<_Real>::digits10);
    auto width = std::cout.precision() + 6;

    std::cout << "\nBernoulli numbers\n";
    auto bern = std::__detail::__bernoulli_vec<_Real>(50);
    for (auto& b : bern)
      std::cout << ' ' << std::setw(width) << b << '\n';

    std::cout << "\nWeighted Bernoulli numbers\n";
    auto wts = std::__detail::__weights(bern);
    for (auto& w : wts)
      std::cout << ' ' << std::setw(width) << w << '\n';

    std::cout << "\nStieltjes partial numerators\n";
    int len = 100;
    auto cf = std::__detail::__stieltjes_cont_frac_seq<_Real>(len);
    for (int k = 0; k < cf.size(); ++k)
      std::cout << ' ' << std::setw(2) << k + 1 << ": "
		<< ' ' << std::setw(width) << cf[k]
		<< ' ' << std::setw(width) << cf[k] / ((k+1) * (k+1) / _Tp{16}) << '\n';

    std::cout << "\nStieltjes partial numerators (progressive)\n";
    auto cfp = std::__detail::__stieltjes_cont_frac_seq_prog<_Real>(len);
    for (int k = 0; k < cfp.size(); ++k)
      std::cout << ' ' << std::setw(2) << k + 1 << ": "
		<< ' ' << std::setw(width) << cfp[k]
		<< ' ' << std::setw(width) << cfp[k] / ((k+1) * (k+1) / _Tp{16}) << '\n';

    std::cout << "\nBinet asymptotic\n";
    for (int k = 1; k <= 5000; ++k)
      {
	auto x = _Tp{0.1Q} * k;
	auto j_as = std::__detail::__binet_asymp(x);
	auto j_cf = std::__detail::__binet_cont_frac(x);
	std::cout << ' ' << std::setw(4) << x
		  << ' ' << std::setw(width) << j_as
		  << ' ' << std::setw(width) << j_cf
		  << ' ' << std::setw(width) << (j_as - j_cf) / j_cf
		  << '\n';
      }
  }

// Test polynomial inversion..
template<typename _Tp>
  void
  test_exp()
  {
    std::cout.precision(std::numeric_limits<_Tp>::digits10);
    auto width = std::cout.precision() + 6;

    std::vector<_Tp> coef;
    _Tp fact = 1;
    coef.push_back(_Tp{1} / fact);
    for (int i = 1; i <= 20; ++i)
      coef.push_back(_Tp{1} / (fact *= _Tp(i)));

    std::cout << "\n exp(x) oefficients:\n";
    for (auto cf : coef)
      std::cout << std::setw(width) << cf << '\n';

    std::vector<_Tp> inv = std::__detail::__inverse_series(coef);
    std::cout << "\n Inverse (hopefully exp(-x)) coefficients:\n";
    for (auto cf : inv)
      std::cout << std::setw(width) << cf << '\n';

/*
    std::cout << "\nTest exp(x)\n";
    __gnu_cxx::_Polynomial<_Tp> expoly(std::begin(coef), std::end(coef));
    __gnu_cxx::_Polynomial<_Tp> rat_numer(std::begin(coef), std::begin(coef) + 10);
    __gnu_cxx::_Polynomial<_Tp> rat_denom(std::begin(inv), std::begin(inv) + 10);
    for (int k = 0; k <= 500; ++k)
      {
	auto x = _Tp{0.1Q} * k;
	auto pexp = expoly(x);
	auto rexp = std::sqrt(rat_numer(x) / rat_denom(x));
	std::cout << ' ' << std::setw(4) << x
		  << ' ' << std::setw(width) << pexp
		  << ' ' << std::setw(width) << pexp - std::exp(x)
		  << ' ' << std::setw(width) << rexp
		  << ' ' << std::setw(width) << rexp - std::exp(x)
		  << '\n';
      }
*/
  }

int
main()
{
  std::cout << "\nTest polynomial inversion\n";
  test_exp<long double>();

  //std::cout << "\nfloat\n=====\n\n";
  //test<float>();

  std::cout << "\ndouble\n======\n";
  test<double>();

  std::cout << "\nlong double\n===========\n";
  test<long double>();

  std::cout << "\n__float128\n===========\n";
  test<__float128>();
}
