/**
 *
 */

#include <cmath>
#include <ext/float128_io.h>
#include <limits>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <vector>
#include <string>
#include <complex>

#include <emsr/summation.h>

#include <wrap_boost.h>

/**
 * @f[
 *    B_{2n}^{(2x)}(x) = -2x\sum_{i=0}^{n-1}
 *        \frac{\binom{2n-1}{2i+2}B_{2i+2}}{2i+2}B_{2n-2i-2}^{(2x)}(x)
 * @f]
 */
template<typename _Tp>
  _Tp
  __bernoulli_2n_2x(unsigned int __n, _Tp __x)
  {
    std::vector<_Tp> _B;

    //_Tp __fact =  * __bernoulli(0);
    _B.push_back(_Tp{1});
    if (__n == _B.size() - 1)
      return _B.back();

    _B.push_back(-__x / _Tp{6});
    if (__n == _B.size() - 1)
      return _B.back();

    _B.push_back(__x * (_Tp{1} + __x * _Tp{5}) / _Tp{60});
    if (__n == _B.size() - 1)
      return _B.back();

    _B.push_back(-__x * (_Tp{4}
		+ __x * (_Tp{21}
		+ __x * _Tp{35})) / _Tp{504});
    if (__n == _B.size() - 1)
      return _B.back();

    _B.push_back(__x * (_Tp{18}
	       + __x * (_Tp{101}
	       + __x * (_Tp{210}
	       + __x * _Tp{175}))) / _Tp{2160});
    if (__n == _B.size() - 1)
      return _B.back();

    _B.push_back(-__x * (_Tp{48}
		+ __x * (_Tp{286}
		+ __x * (_Tp{671}
		+ __x * (_Tp{770}
		+ __x * _Tp{385})))) / _Tp{3168});
    if (__n == _B.size() - 1)
      return _B.back();

    _B.push_back(__x * (_Tp{33168}
	       + __x * (_Tp{207974}
	       + __x * (_Tp{531531}
	       + __x * (_Tp{715715}
	       + __x * (_Tp{525525}
	       + __x * _Tp{175175}))))) / _Tp{786240});
    if (__n == _B.size() - 1)
      return _B.back();

    for (unsigned __k = _B.size(); __k <= __n; ++__k)
      {
	_B.push_back(_Tp{0});
	for (unsigned __i = 0; __i < _B.size() - 1; ++__i)
	  {
	    _B.back() += std::__detail::__binomial<_Tp>(2 * __k - 1, 2 * __i + 1)
			 * std::__detail::__bernoulli_2n<_Tp>(__i + 1)
			 * _B[__k - __i - 1] / _Tp(2 * __i + 2);
	  }
	_B.back() *= -_Tp{2} * __x;
      }

    return _B.back();
  }

  /**
   * Buhring equation modes.
   */
  enum __buhring_mode
  {
    automatic,
    equation2p7,
    equation3p1
  };

  /**
   * Factorial ratio.
   */
  template<typename _Tp>
    _Tp
    __factorial_ratio(long long __m, long long __n)
    {
      if (__m == __n)
	return _Tp{1};
      else if (__m < 0 && __n >= 0)
	return std::numeric_limits<_Tp>::quiet_NaN();
      else if (__m >= 0 && __n < 0)
	return _Tp{0};
      else if ((__m == 0 || __m == 1)
	    && __n < std::__detail::_S_num_factorials<_Tp>)
	return _Tp{1} / std::__detail::__factorial<_Tp>(__n);
      else if ((__n == 0 || __n == 1)
	    && __m < std::__detail::_S_num_factorials<_Tp>)
	return std::__detail::__factorial<_Tp>(__m);
      else
	{
	  // Try a running product.
	  auto __pmin = __m;
	  auto __pmax = __n;
	  if (__pmax < __pmin)
	    std::swap(__pmin, __pmax);
	  long long __p = 1;
	  for (; __pmin < __pmax; ++__pmin)
	    __p *= __pmin;
	  if (__m < __n)
	    return _Tp{1} / _Tp{__p};
	  else
	    return _Tp{__p};
	}
    }


  /**
   * Return the four-gamma ratio.
   * @f[
   *    \frac{\Gamma(a+n)\Gamma(b+n)}{\Gamma(c+n)\Gamma(a+b-c+n)}
   *      = 1 + sum_{m=1}^{M}\frac{(c-a)_m(c-b)_m}{m!(1+c-a-b-n)_m}
   *          + O(n^{-M-1})
   *      = 1 + sum_{m=1}^{M}\frac{(a-c)_m(b-c)_m}{m!(1-c-n)_m}
   *          + O(n^{-M-1})
   * @f]
   */
  template<typename _Tn, typename _Tp>
    _Tp
    __gamma_ratio_buhring(_Tn __n, _Tp __a, _Tp __b, _Tp __c,
			  int _S_M = 20, __buhring_mode mode = automatic)
    {
      const auto _S_eps = 10 * __gnu_cxx::__epsilon(__a);
      if (mode == equation2p7
	|| (mode == automatic && std::real(1 - __c - __n) < _Tp{0}))
	{
	  //emsr::WenigerDeltaSum<emsr::BasicSum<_Tp>> __sum;
	  emsr::BasicSum<_Tp> __sum;
	  auto __term = _Tp{1};
	  __sum += _Tp{__term};
	  auto __ca = __c - __a;
	  auto __cb = __c - __b;
	  auto __cabn = _Tp{1} + __c - __a - __b - _Tp(__n);
	  for (int __m = 1; __m <= _S_M; ++__m)
	    {
	      auto __prev = __term;
	      if (__cabn == _Tp{0})
		break;
	      __term *= __ca / __m * __cb / __cabn;
	      __sum += __term;
	      if (mode == automatic
		&& (std::abs(__term) < _S_eps * __sum()
		 || std::abs(__term) > std::abs(__prev)))
		break;
	      ++__ca;
	      ++__cb;
	      ++__cabn;
	    }
	  return __sum();
	}
      else if (mode == equation3p1
	   || (mode == automatic
	    && std::real(1 + __c - __a - __b - _Tp(__n)) < _Tp{0}))
	{
	  //emsr::WenigerDeltaSum<emsr::BasicSum<_Tp>> __sum;
	  emsr::BasicSum<_Tp> __sum;
	  auto __term = _Tp{1};
	  __sum += _Tp{__term};
	  auto __ac = __a - __c;
	  auto __bc = __b - __c;
	  auto __cn = _Tp{1} - __c - _Tp(__n);
	  for (int __m = 1; __m <= _S_M; ++__m)
	    {
	      auto __prev = __term;
	      if (__cn == _Tp{0})
		break;
	      __term *= __ac / __m * __bc / __cn;
	      __sum += __term;
	      if (mode == automatic
		&& (std::abs(__term) < _S_eps * __sum()
		 || std::abs(__term) > std::abs(__prev)))
		break;
	      if (mode == automatic && std::abs(__term) > std::abs(__prev))
		break;
	      ++__ac;
	      ++__bc;
	      ++__cn;
	    }
	  return __sum();
	}
      else
	{
	  // Average 2.7 and 3.1 sums.
	  emsr::BasicSum<_Tp> __sum2p7;
	  auto __term2p7 = _Tp{1};
	  __sum2p7 += _Tp{__term2p7};
	  auto __ca = __c - __a;
	  auto __cb = __c - __b;
	  auto __cabn = _Tp{1} + __c - __a - __b - _Tp(__n);

	  emsr::BasicSum<_Tp> __sum3p1;
	  auto __term3p1 = _Tp{1};
	  __sum3p1 += _Tp{__term3p1};
	  auto __ac = __a - __c;
	  auto __bc = __b - __c;
	  auto __cn = _Tp{1} - __c - _Tp(__n);

	  bool __conv2p7 = false;
	  bool __conv3p1 = false;
	  for (int __m = 1; __m <= _S_M; ++__m)
	    {
	      if (__cabn == _Tp{0})
		__conv2p7 = true;
	      else if (!__conv2p7)
		{
		  auto __prev2p7 = __term2p7;
		  __term2p7 *= __ca / __m * __cb / __cabn;
		  __sum2p7 += __term2p7;
		  ++__ca;
		  ++__cb;
		  ++__cabn;
		  if (std::abs(__term2p7) < _S_eps * __sum2p7()
		    || std::abs(__term2p7) > std::abs(__prev2p7))
		    __conv2p7 = true;
		}

	      if (__cn == _Tp{0})
		__conv3p1 = true;
	      else if (!__conv3p1)
		{
		  auto __prev3p1 = __term3p1;
		  __term3p1 *= __ac / __m * __bc / __cn;
		  __sum3p1 += __term3p1;
		  ++__ac;
		  ++__bc;
		  ++__cn;
		  if (std::abs(__term3p1) < _S_eps * __sum3p1()
		    || std::abs(__term3p1) > std::abs(__prev3p1))
		    __conv3p1 = true;
		}

	      if (__conv2p7 && __conv3p1)
		break;
	    }
	  return (__sum2p7() + __sum3p1()) / _Tp{2}
		* __gnu_cxx::sin_pi(__c + _Tp(__n))
		* __gnu_cxx::sin_pi(__a + __b - __c + _Tp(__n))
		/ __gnu_cxx::sin_pi(__a + __n)
		* __gnu_cxx::sin_pi(__b + _Tp(__n));
	}
    }

  /**
   * Juggle gamma special cases - from hyperg.
   */
  template<typename _Tp>
    _Tp
    __hyperg_reflect(_Tp __a, _Tp __b, _Tp __c)
    {
      const auto __d = __c - __a - __b;
      const auto __intd = emsr::fp_is_integer(__d);
      auto __F1 = _Tp{0};
      if (__intd)
	{
	  if (__intd() <= 0)
	    {
	      __F1 = _Tp{0};
	    }
	  else
	    {
	      __F1 = _Tp{0};
	    }
	}
      else
	{
	  __F1 = _Tp{0};
	}
    }

/**
 * Return the gamma ratio by asymptotic series:
 * @f[
 *    \frac{\Gamma(z+a)}{\Gamma(z+b)}
 *  ~ z^{a-b}\sum_{n=0}^{\infty}\frac{(b-a)_n}{n!}
 *         {}_2F_0(-n, z+b;;1/z)
 * @f]
 * The hypergeometric function @f$ {}_2F_0(a,b;;z) @f$
 * with nonpositive numerator parameter n
 * is a polynomial with a simple TTRR:
 * @f[
 *  z {}_2F_0(-n-1,z+b;;1/z) = -(n+b) {}_2F_0(-n,z+b;;1/z)
 *                           + n {}_2F_0(-n+1,z+b;;1/z)
 * @f]
 * where @f$ {}_2F_0(-0,z+b;;1/z) = 1 @f$
 * and @f$ {}_2F_0(-1,z+b;;1/z) = -b/z @f$.
 */
template<typename _Tp>
  _Tp
  gamma_ratio_asymp_2f0(_Tp __a, _Tp __b, _Tp __z)
  {
    const auto _S_max_iter = 1000;
    auto __fact = _Tp{1};
    auto _Fnm1 = _Tp{1};
    auto __sum = __fact * _Fnm1;
    __fact *= (__b - __a);
    auto _Fn = -__b / __z;
    auto __term = __fact * _Fn;
    __sum += __term;
    auto __prev_term = std::abs(__term);
    for (auto __n = 2; __n < _S_max_iter; ++__n)
      {
	__fact *= (__b - __a + __n - 1) / __n;
	auto _Fnp1 = (-(__n + __b) * _Fn + __n * _Fnm1) / __z;
	__term = __fact * _Fnp1;
	if (std::abs(__term) > __prev_term)
	  break;
	__sum += __term;
	_Fnm1 = _Fn;
	_Fn = _Fnp1;
      }
    return __sum * std::pow(__z, __a - __b);
  }

/**
 * Return the gamma ratio by asymptotic series:
 * @f[
 *    \frac{\Gamma(z+a)}{\Gamma(z+b)}
 *  ~ z^{(a-b)}\sum_{n=0}^{\infty}\frac{(a-b)_n}{n!}
 *         {}_1F_1(-n; z+a; z)
 * @f]
 * The hypergeometric function @f$ {}_1F_1(a;b;z) @f$
 * with nonpositive integral numerator parameter n
 * is a polynomial with a simple TTRR:
 * @f[
 *  (z+a+n) {}_1F_1(-n-1;z+a;z) = -(2n+a) {}_1F_1(-n;z+a;z)
 *                              - n {}_1F_1(-n+1;z+a;z)
 * @f]
 * where @f$ {}_1F_1(-0;z+a;z) = 1 @f$ and @f$ {}_1F_1(-1;z+a;z) = a/(z+a) @f$.
 */
template<typename _Tp>
  _Tp
  gamma_ratio_asymp_1f1(_Tp __z, _Tp __a, _Tp __b)
  {
    return _Tp{0};
  }

/**
 * Return the gamma ratio by asymptotic series:
 * @f[
 *    \frac{\Gamma(z + a)}{\Gamma(z + b)}
 *    ~ z^{a-b}\sum_{n=0}^{\infty} \frac{(a-b+1-n)_n B_n^{a-b+1}(a)}{n!}z^{-n}
 * @f]
 * where @f$@f$ is the Norlund or generalized Bernoulli polynomial.
 */
template<typename _Tp>
  _Tp
  gamma_ratio_asymp_erdelyi_tricomi(_Tp __z, _Tp __a, _Tp __b)
  {
    return _Tp{0};
  }

/**
 * Return the ratio
 * @f[
 *   \frac{\Gamma(z+a)}{\Gamma(z)}
 * @f]
 * for large @f$ a @f$.
 */
template<typename _Tp>
  _Tp
  __gamma_ratio_asymp(_Tp __z, _Tp __a)
  {
    return _Tp{0};
  }

template<typename _Tp>
  void
  test_gamma_ratio(_Tp proto = _Tp{})
  {
    //using _Val = _Tp;
    //using _Real = emsr::num_traits_t<_Val>;

    std::vector<_Tp> parm{_Tp{0.25}, _Tp{0.5}, _Tp{1}, _Tp{2}, _Tp{5}};

    std::cout.precision(__gnu_cxx::__digits10(proto));
    std::cout << std::showpoint << std::scientific;
    auto width = 8 + std::cout.precision();

    std::cout << '\n'
	      << ' ' << std::setw(width) << "a"
	      << ' ' << std::setw(width) << "b"
	      << ' ' << std::setw(width) << "c"
	      << ' ' << std::setw(width) << "z"
	      << ' ' << std::setw(width) << "hyperg0"
	      << ' ' << std::setw(width) << "hyperg"
	      << '\n';
    int i_min = -200;
    for (auto a : parm)
      for (auto b : parm)
	for (auto c : parm)
	  for (int i = i_min; i <= +200; ++i)
	    {
	      auto n = _Tp{0.1Q} * i;
	      //auto gamrat0 = __gamma_ratio_log(a, b, c, z);
	      auto gamrat = __gamma_ratio_buhring(n, a, b, c);
	      std::cout << ' ' << std::setw(width) << n
			<< ' ' << std::setw(width) << a
			<< ' ' << std::setw(width) << b
			<< ' ' << std::setw(width) << c
			//<< ' ' << std::setw(width) << gamrat0
			<< ' ' << std::setw(width) << gamrat
			//<< ' ' << std::setw(width) << gamrat - gamrat0
			<< '\n';
	    }
  }

template<typename _Tp>
  void
  test_gamma_ratio_buhring(_Tp proto = _Tp{})
  {
    std::cout.precision(__gnu_cxx::__digits10(proto));
    std::cout << std::showpoint << std::scientific;
    auto width = 8 + std::cout.precision();

    _Tp a = -11.7;
    _Tp b = -11.2;
    _Tp c = -11.4;

    std::cout << '\n';

    std::cout << '\n';
    for (int M = 1; M <= 10; ++M)
      {
	auto b10_2p7 = __gamma_ratio_buhring(10, a, b, c, M, equation2p7);
	auto b10_3p1 = __gamma_ratio_buhring(10, a, b, c, M, equation3p1);
	std::cout << std::setw(2) << M
		  << std::setw(width) << b10_2p7
		  << std::setw(width) << b10_3p1
		  << '\n';
      }

    std::cout << '\n';
    for (int M = 1; M <= 10; ++M)
      {
	auto b20_2p7 = __gamma_ratio_buhring(20, a, b, c, M, equation2p7);
	auto b20_3p1 = __gamma_ratio_buhring(20, a, b, c, M, equation3p1);
	std::cout << std::setw(2) << M
		  << std::setw(width) << b20_2p7
		  << std::setw(width) << b20_3p1
		  << '\n';
      }

    a = 11.7;
    b = 11.2;
    c = 11.4;

    std::cout << '\n';

    std::cout << '\n';
    for (int M = 1; M <= 10; ++M)
      {
	auto b10_2p7 = __gamma_ratio_buhring(-15, a, b, c, M, equation2p7);
	auto b10_3p1 = __gamma_ratio_buhring(-15, a, b, c, M, equation3p1);
	std::cout << std::setw(2) << M
		  << std::setw(width) << b10_2p7
		  << std::setw(width) << b10_3p1
		  << '\n';
      }

    std::cout << '\n';
    for (int M = 1; M <= 10; ++M)
      {
	auto b20_2p7 = __gamma_ratio_buhring(-5, a, b, c, M, equation2p7);
	auto b20_3p1 = __gamma_ratio_buhring(-5, a, b, c, M, equation3p1);
	std::cout << std::setw(2) << M
		  << std::setw(width) << b20_2p7
		  << std::setw(width) << b20_3p1
		  << '\n';
      }
  }

template<typename _Tp>
  void
  test_bernoulli_2n_2x(_Tp proto = _Tp{})
  {
    std::cout.precision(__gnu_cxx::__digits10(proto));
    std::cout << std::showpoint << std::scientific;
    auto w = 8 + std::cout.precision();

    std::cout << '\n';
    std::cout << '\n';
    for (int i = 0; i <= 100; ++i)
      {
	auto x = i * 0.1;
	std::cout << ' ' << std::setw(w) << x;
	for (int n = 0; n <= 10; ++n)
	  std::cout << ' ' << std::setw(w) << __bernoulli_2n_2x(n, x);
	std::cout << '\n';
      }
  }


/**
 * @f[
 *   \frac{\Gamma(z + a)}{\Gamma(z + b)} = (z + a - \rho)^{a-b}
 *      \sum_{k=0}^{N-1}\frac{(b - a)_{2k}(z +a - \rho)^{-2k}}{(2k)!}
 *          B_{(2\rho)}^{2k}(\rho)
 * @f]
 */
template<typename _Tp>
  _Tp
  __gamma_ratio_asymp_field(_Tp __z, _Tp __a, _Tp __b)
  {
    const auto __bma = __b - __a;
    const auto __rho = (_Tp{1} + __a + __b) / _Tp{2};
    const auto __arg = __z + __a - __rho;
    const auto __arg2 = __arg * __arg;
    auto __fact = _Tp{1};
    auto __ratio = _Tp{1};
    for (int __k = 1; __k < 7; ++__k)
      {
	__fact *= ((__bma + 2 * __k - 2) / _Tp(2 * __k - 1))
		* ((__bma + 2 * __k - 1) / _Tp(2 * __k))
		/ __arg2;
	__ratio += __fact * __bernoulli_2n_2x(__k, __rho);
      }
    return __ratio / std::pow(__arg, __bma);
  }


/**
 * @f[
 *    \binom{\nu}{k} = \frac{(-1)^k (k-\nu/2)^{-(\nu+1)}}{\Gamma(-\nu)}
 *    \sum_{i=0}^{\infty}
 *      \frac{B_{2i}^{(2\rho)}(\rho)(\nu+1)_{2i}}
 *           {(2i)!(k-\nu/2)^{2i}}
 * @f]
 * where @f$ \rho = -\nu/2 @f$.
 */
template<typename _Tp>
  _Tp
  __binomial_asymp_field(_Tp __nu, unsigned int __k)
  {
    constexpr auto _S_max_iter = 1000U;
    const auto _S_eps = __gnu_cxx::__epsilon(__nu);
    const auto __rho = -__nu / _Tp{2};
    const auto __pocharg = __nu + _Tp{1};
    const auto __powarg = _Tp(__k) + __rho;
    const auto __powarg2 = __powarg * __powarg;
    auto __sum = __bernoulli_2n_2x(0, __rho);
    auto __fact = _Tp{1};
    for (int __i = 1; __i < _S_max_iter; ++__i)
      {
	__fact *= ((__pocharg + _Tp(2 * __i - 2)) / _Tp(2 * __i - 1))
		* ((__pocharg + _Tp(2 * __i - 1)) / _Tp(2 * __i))
		/ __powarg2;
	const auto __term = __bernoulli_2n_2x(__i, __rho) * __fact;
	__sum += __term;
	if (std::abs(__term) < _S_eps * std::abs(__sum))
	  break;
      }
    return _Tp((__k & 1) ? -1 : +1) / std::pow(__powarg, __pocharg)
	 * __sum * __gamma_recip(-__nu);
  }

int
main()
{
  test_gamma_ratio_buhring(1.0);

  test_gamma_ratio(1.0f);
  
  test_gamma_ratio(1.0);

  test_gamma_ratio(1.0l);

  //test_gamma_ratio(1.0q);

  test_bernoulli_2n_2x(1.0);
}
