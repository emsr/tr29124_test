/**
 *
 */

#include <cmath>
#include <limits>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <vector>
#include <string>
#include <complex>

#include <emsr/summation.h>
#include <emsr/float128_io.h>
#include <emsr/special_functions.h>

#include <wrap_boost.h>

/**
 * @f[
 *    B_{2n}^{(2x)}(x) = -2x\sum_{i=0}^{n-1}
 *        \frac{\binom{2n-1}{2i+2}B_{2i+2}}{2i+2}B_{2n-2i-2}^{(2x)}(x)
 * @f]
 */
template<typename _Tp>
  _Tp
  bernoulli_2n_2x(unsigned int n, _Tp x)
  {
    std::vector<_Tp> _B;

    //_Tp fact =  * bernoulli(0);
    _B.push_back(_Tp{1});
    if (n == _B.size() - 1)
      return _B.back();

    _B.push_back(-x / _Tp{6});
    if (n == _B.size() - 1)
      return _B.back();

    _B.push_back(x * (_Tp{1} + x * _Tp{5}) / _Tp{60});
    if (n == _B.size() - 1)
      return _B.back();

    _B.push_back(-x * (_Tp{4}
		+ x * (_Tp{21}
		+ x * _Tp{35})) / _Tp{504});
    if (n == _B.size() - 1)
      return _B.back();

    _B.push_back(x * (_Tp{18}
	       + x * (_Tp{101}
	       + x * (_Tp{210}
	       + x * _Tp{175}))) / _Tp{2160});
    if (n == _B.size() - 1)
      return _B.back();

    _B.push_back(-x * (_Tp{48}
		+ x * (_Tp{286}
		+ x * (_Tp{671}
		+ x * (_Tp{770}
		+ x * _Tp{385})))) / _Tp{3168});
    if (n == _B.size() - 1)
      return _B.back();

    _B.push_back(x * (_Tp{33168}
	       + x * (_Tp{207974}
	       + x * (_Tp{531531}
	       + x * (_Tp{715715}
	       + x * (_Tp{525525}
	       + x * _Tp{175175}))))) / _Tp{786240});
    if (n == _B.size() - 1)
      return _B.back();

    for (unsigned k = _B.size(); k <= n; ++k)
      {
	_B.push_back(_Tp{0});
	for (unsigned i = 0; i < _B.size() - 1; ++i)
	  {
	    _B.back() += emsr::detail::binomial<_Tp>(2 * k - 1, 2 * i + 1)
			 * emsr::detail::bernoulli_2n<_Tp>(i + 1)
			 * _B[k - i - 1] / _Tp(2 * i + 2);
	  }
	_B.back() *= -_Tp{2} * x;
      }

    return _B.back();
  }

  /**
   * Buhring equation modes.
   */
  enum buhring_mode
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
    factorial_ratio(long long m, long long n)
    {
      if (m == n)
	return _Tp{1};
      else if (m < 0 && n >= 0)
	return std::numeric_limits<_Tp>::quiet_NaN();
      else if (m >= 0 && n < 0)
	return _Tp{0};
      else if ((m == 0 || m == 1)
	    && n < emsr::detail::s_num_factorials<_Tp>)
	return _Tp{1} / emsr::detail::factorial<_Tp>(n);
      else if ((n == 0 || n == 1)
	    && m < emsr::detail::s_num_factorials<_Tp>)
	return emsr::detail::factorial<_Tp>(m);
      else
	{
	  // Try a running product.
	  auto pmin = m;
	  auto pmax = n;
	  if (pmax < pmin)
	    std::swap(pmin, pmax);
	  long long p = 1;
	  for (; pmin < pmax; ++pmin)
	    p *= pmin;
	  if (m < n)
	    return _Tp{1} / _Tp{p};
	  else
	    return _Tp{p};
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
    gamma_ratio_buhring(_Tn n, _Tp a, _Tp b, _Tp c,
			  int s_M = 20, buhring_mode mode = automatic)
    {
      const auto s_eps = 10 * emsr::epsilon(a);
      if (mode == equation2p7
	|| (mode == automatic && std::real(1 - c - n) < _Tp{0}))
	{
	  //emsr::WenigerDeltaSum<emsr::BasicSum<_Tp>> sum;
	  emsr::BasicSum<_Tp> sum;
	  auto term = _Tp{1};
	  sum += _Tp{term};
	  auto ca = c - a;
	  auto cb = c - b;
	  auto cabn = _Tp{1} + c - a - b - _Tp(n);
	  for (int m = 1; m <= s_M; ++m)
	    {
	      auto prev = term;
	      if (cabn == _Tp{0})
		break;
	      term *= ca / m * cb / cabn;
	      sum += term;
	      if (mode == automatic
		&& (std::abs(term) < s_eps * sum()
		 || std::abs(term) > std::abs(prev)))
		break;
	      ++ca;
	      ++cb;
	      ++cabn;
	    }
	  return sum();
	}
      else if (mode == equation3p1
	   || (mode == automatic
	    && std::real(1 + c - a - b - _Tp(n)) < _Tp{0}))
	{
	  //emsr::WenigerDeltaSum<emsr::BasicSum<_Tp>> sum;
	  emsr::BasicSum<_Tp> sum;
	  auto term = _Tp{1};
	  sum += _Tp{term};
	  auto ac = a - c;
	  auto bc = b - c;
	  auto cn = _Tp{1} - c - _Tp(n);
	  for (int m = 1; m <= s_M; ++m)
	    {
	      auto prev = term;
	      if (cn == _Tp{0})
		break;
	      term *= ac / m * bc / cn;
	      sum += term;
	      if (mode == automatic
		&& (std::abs(term) < s_eps * sum()
		 || std::abs(term) > std::abs(prev)))
		break;
	      if (mode == automatic && std::abs(term) > std::abs(prev))
		break;
	      ++ac;
	      ++bc;
	      ++cn;
	    }
	  return sum();
	}
      else
	{
	  // Average 2.7 and 3.1 sums.
	  emsr::BasicSum<_Tp> sum2p7;
	  auto term2p7 = _Tp{1};
	  sum2p7 += _Tp{term2p7};
	  auto ca = c - a;
	  auto cb = c - b;
	  auto cabn = _Tp{1} + c - a - b - _Tp(n);

	  emsr::BasicSum<_Tp> sum3p1;
	  auto term3p1 = _Tp{1};
	  sum3p1 += _Tp{term3p1};
	  auto ac = a - c;
	  auto bc = b - c;
	  auto cn = _Tp{1} - c - _Tp(n);

	  bool conv2p7 = false;
	  bool conv3p1 = false;
	  for (int m = 1; m <= s_M; ++m)
	    {
	      if (cabn == _Tp{0})
		conv2p7 = true;
	      else if (!conv2p7)
		{
		  auto prev2p7 = term2p7;
		  term2p7 *= ca / m * cb / cabn;
		  sum2p7 += term2p7;
		  ++ca;
		  ++cb;
		  ++cabn;
		  if (std::abs(term2p7) < s_eps * sum2p7()
		    || std::abs(term2p7) > std::abs(prev2p7))
		    conv2p7 = true;
		}

	      if (cn == _Tp{0})
		conv3p1 = true;
	      else if (!conv3p1)
		{
		  auto prev3p1 = term3p1;
		  term3p1 *= ac / m * bc / cn;
		  sum3p1 += term3p1;
		  ++ac;
		  ++bc;
		  ++cn;
		  if (std::abs(term3p1) < s_eps * sum3p1()
		    || std::abs(term3p1) > std::abs(prev3p1))
		    conv3p1 = true;
		}

	      if (conv2p7 && conv3p1)
		break;
	    }
	  return (sum2p7() + sum3p1()) / _Tp{2}
		* emsr::sin_pi(c + _Tp(n))
		* emsr::sin_pi(a + b - c + _Tp(n))
		/ emsr::sin_pi(a + n)
		* emsr::sin_pi(b + _Tp(n));
	}
    }

  /**
   * Juggle gamma special cases - from hyperg.
   */
  template<typename _Tp>
    _Tp
    hyperg_reflect(_Tp a, _Tp b, _Tp c)
    {
      const auto d = c - a - b;
      const auto intd = emsr::fp_is_integer(d);
      auto F1 = _Tp{0};
      if (intd)
	{
	  if (intd() <= 0)
	    {
	      F1 = _Tp{0};
	    }
	  else
	    {
	      F1 = _Tp{0};
	    }
	}
      else
	{
	  F1 = _Tp{0};
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
  gamma_ratio_asymp_2f0(_Tp a, _Tp b, _Tp z)
  {
    const auto s_max_iter = 1000;
    auto fact = _Tp{1};
    auto _Fnm1 = _Tp{1};
    auto sum = fact * _Fnm1;
    fact *= (b - a);
    auto _Fn = -b / z;
    auto term = fact * _Fn;
    sum += term;
    auto prev_term = std::abs(term);
    for (auto n = 2; n < s_max_iter; ++n)
      {
	fact *= (b - a + n - 1) / n;
	auto _Fnp1 = (-(n + b) * _Fn + n * _Fnm1) / z;
	term = fact * _Fnp1;
	if (std::abs(term) > prev_term)
	  break;
	sum += term;
	_Fnm1 = _Fn;
	_Fn = _Fnp1;
      }
    return sum * std::pow(z, a - b);
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
  gamma_ratio_asymp_1f1(_Tp z, _Tp a, _Tp b)
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
  gamma_ratio_asymp_erdelyi_tricomi(_Tp z, _Tp a, _Tp b)
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
  gamma_ratio_asymp(_Tp z, _Tp a)
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

    std::cout.precision(emsr::digits10(proto));
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
	      //auto gamrat0 = gamma_ratio_log(a, b, c, z);
	      auto gamrat = gamma_ratio_buhring(n, a, b, c);
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
    std::cout.precision(emsr::digits10(proto));
    std::cout << std::showpoint << std::scientific;
    auto width = 8 + std::cout.precision();

    _Tp a = -11.7;
    _Tp b = -11.2;
    _Tp c = -11.4;

    std::cout << '\n';

    std::cout << '\n';
    for (int M = 1; M <= 10; ++M)
      {
	auto b10_2p7 = gamma_ratio_buhring(10, a, b, c, M, equation2p7);
	auto b10_3p1 = gamma_ratio_buhring(10, a, b, c, M, equation3p1);
	std::cout << std::setw(2) << M
		  << std::setw(width) << b10_2p7
		  << std::setw(width) << b10_3p1
		  << '\n';
      }

    std::cout << '\n';
    for (int M = 1; M <= 10; ++M)
      {
	auto b20_2p7 = gamma_ratio_buhring(20, a, b, c, M, equation2p7);
	auto b20_3p1 = gamma_ratio_buhring(20, a, b, c, M, equation3p1);
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
	auto b10_2p7 = gamma_ratio_buhring(-15, a, b, c, M, equation2p7);
	auto b10_3p1 = gamma_ratio_buhring(-15, a, b, c, M, equation3p1);
	std::cout << std::setw(2) << M
		  << std::setw(width) << b10_2p7
		  << std::setw(width) << b10_3p1
		  << '\n';
      }

    std::cout << '\n';
    for (int M = 1; M <= 10; ++M)
      {
	auto b20_2p7 = gamma_ratio_buhring(-5, a, b, c, M, equation2p7);
	auto b20_3p1 = gamma_ratio_buhring(-5, a, b, c, M, equation3p1);
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
    std::cout.precision(emsr::digits10(proto));
    std::cout << std::showpoint << std::scientific;
    auto w = 8 + std::cout.precision();

    std::cout << '\n';
    std::cout << '\n';
    for (int i = 0; i <= 100; ++i)
      {
	auto x = i * 0.1;
	std::cout << ' ' << std::setw(w) << x;
	for (int n = 0; n <= 10; ++n)
	  std::cout << ' ' << std::setw(w) << bernoulli_2n_2x(n, x);
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
  gamma_ratio_asymp_field(_Tp z, _Tp a, _Tp b)
  {
    const auto bma = b - a;
    const auto rho = (_Tp{1} + a + b) / _Tp{2};
    const auto arg = z + a - rho;
    const auto arg2 = arg * arg;
    auto fact = _Tp{1};
    auto ratio = _Tp{1};
    for (int k = 1; k < 7; ++k)
      {
	fact *= ((bma + 2 * k - 2) / _Tp(2 * k - 1))
		* ((bma + 2 * k - 1) / _Tp(2 * k))
		/ arg2;
	ratio += fact * bernoulli_2n_2x(k, rho);
      }
    return ratio / std::pow(arg, bma);
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
  binomial_asymp_field(_Tp nu, unsigned int k)
  {
    constexpr auto s_max_iter = 1000U;
    const auto s_eps = emsr::epsilon(nu);
    const auto rho = -nu / _Tp{2};
    const auto pocharg = nu + _Tp{1};
    const auto powarg = _Tp(k) + rho;
    const auto powarg2 = powarg * powarg;
    auto sum = bernoulli_2n_2x(0, rho);
    auto fact = _Tp{1};
    for (int i = 1; i < s_max_iter; ++i)
      {
	fact *= ((pocharg + _Tp(2 * i - 2)) / _Tp(2 * i - 1))
		* ((pocharg + _Tp(2 * i - 1)) / _Tp(2 * i))
		/ powarg2;
	const auto term = bernoulli_2n_2x(i, rho) * fact;
	sum += term;
	if (std::abs(term) < s_eps * std::abs(sum))
	  break;
      }
    return _Tp((k & 1) ? -1 : +1) / std::pow(powarg, pocharg)
	 * sum * gamma_recip(-nu);
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
