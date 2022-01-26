/**
 *
 */

#include <limits>
#include <iostream>
#include <iomanip>
#include <cmath>

#include <emsr/continued_fractions.h>
#include <emsr/polynomial.h>
#include <emsr/horner.h>
#include <emsr/special_functions.h>

#include <wrap_boost.h>

/**
 * 
 */
template<typename _Tp>
  _Tp
  polygamma_series(unsigned int m, _Tp x)
  {
    constexpr int s_max_iter = 1000;
    const auto s_eps = emsr::epsilon(x);
    auto sum = _Tp{0};
    for (int k = 0; k < s_max_iter; ++k)
      {
	auto term = std::pow(x + _Tp(k), _Tp(m + 1));
	sum += term;
	if (std::abs(term) < s_eps * std::abs(sum))
	  break;
      }
    return (m & 1 ? +1 : -1) * emsr::factorial<_Tp>(m) * sum;
  }

/**
 * 
 */
template<typename _Tp>
  _Tp
  trigamma_cont_frac(_Tp x)
  {
    const auto s_2pi = emsr::tau_v<_Tp>;
    const auto s_12pi = _Tp{1} / s_2pi / _Tp{6};

    auto a
      = [s_12pi](std::size_t k, _Tp)
	{
	  if (k == 1)
	    return s_12pi;
	  else
	    {
	      auto kk = _Tp(k * k);
	      return kk * (kk - _Tp{1}) / _Tp{4} / (_Tp{4} * kk - _Tp{1});
	    }
	};
    using _AFun = decltype(a);

    auto b
      = [](std::size_t k, _Tp x)
	{ return k == 0 ? _Tp{0} : (k & 1 ? x * x : _Tp{1}); };
    using _BFun = decltype(b);

    auto w
      = [a, b](std::size_t k, _Tp x)
	{
	  auto arg = _Tp{4} * a(k + 1, x) / x / x + _Tp{1};
	  return b(k, x) * (std::sqrt(arg) - _Tp{1}) / _Tp{2};
	};
    using _WFun = decltype(w);

    emsr::LentzContinuedFraction<_Tp, _AFun, _BFun, _WFun>
      G1Frac(a, b, w);

    auto g1 = G1Frac(x);
    auto rx = _Tp{1} / x;

    return rx + rx * rx / _Tp{2} + s_2pi * rx * g1;
  }

/**
 * 
 */
template<typename _Tp>
  _Tp
  tetragamma_cont_frac(_Tp x)
  {
    const auto s_2pi = emsr::tau_v<_Tp>;
    const auto s_8pi2 = _Tp{1} / s_2pi / s_2pi / _Tp{2};

    auto a
      = [s_8pi2](std::size_t k, _Tp)
	{
	  if (k == 1)
	    return s_8pi2;
	  else
	    {
	      auto j = k & 1 ? (k - 1) / 2 : k / 2;
	      auto aa = _Tp(j) * _Tp(j + 1) / _Tp{2} / _Tp(2 * j + 1);
	      return (k & 1 ? aa * _Tp(j + 1) : aa * _Tp(j));
	    }
	};
    using _AFun = decltype(a);

    auto b
      = [](std::size_t k, _Tp x)
	{ return k == 0 ? _Tp{0} : (k & 1 ? x * x : _Tp{1}); };
    using _BFun = decltype(b);

    auto w
      = [a, b](std::size_t k, _Tp x)
	{
	  auto arg = _Tp{4} * a(k + 1, x) / x / x + _Tp{1};
	  return b(k, x) * (std::sqrt(arg) - _Tp{1}) / _Tp{2};
	};
    using _WFun = decltype(w);

    emsr::LentzContinuedFraction<_Tp, _AFun, _BFun, _WFun>
      G2Frac(a, b, w);

    auto g2 = G2Frac(x);

    auto fg = s_2pi / x;
    auto xm2 = _Tp{1} / (x * x);

    return -xm2 - xm2 / x - fg * fg * g2;
  }

/**
 * 
 */
template<typename _Tp>
  void
  test_trigamma(_Tp proto)
  {
    std::cout.precision(emsr::digits10(proto));
    std::cout << std::showpoint << std::scientific;
    auto width = 8 + std::cout.precision();

    std::cout << '\n';
    std::cout << '\n';
    for (int i = 200; i <= 1000; ++i)
      {
	auto x = i * (0.01);
	std::cout << ' ' << x
		  << ' ' << std::setw(width) << trigamma_cont_frac(x)
		  << '\n';
      }

    std::cout << '\n';
    for (int i = 100; i <= 950; ++i)
      {
	auto x = i * (0.1);
	std::cout << ' ' << x
		  << ' ' << std::setw(width) << trigamma_cont_frac(x)
		  << '\n';
      }
  }

/**
 * 
 */
template<typename _Tp>
  void
  test_tetragamma(_Tp proto)
  {
    std::cout.precision(emsr::digits10(proto));
    std::cout << std::showpoint << std::scientific;
    auto width = 8 + std::cout.precision();

    std::cout << '\n';
    std::cout << '\n';
    for (int i = 200; i <= 1000; ++i)
      {
	auto x = i * (0.01);
	std::cout << ' ' << x
		  << ' ' << std::setw(width) << tetragamma_cont_frac(x)
		  << '\n';
      }

    std::cout << '\n';
    for (int i = 100; i <= 950; ++i)
      {
	auto x = i * (0.1);
	std::cout << ' ' << x
		  << ' ' << std::setw(width) << tetragamma_cont_frac(x)
		  << '\n';
      }
  }

/**
 * Get the polygamma reflection polynomial.
 *
 * From boost (c = cos(pi*x), polynomials are even - i.e. in c*c):
 *  {-1}
 *  c * {2}
 *  { -2, -4 }
 *  c * { 16, 8 }
 *  { -16, -88, -16 }
 *  c * { 272, 416, 32 }
 *  { -272, -2880, -1824, -64 }
 *  c * { 7936, 24576, 7680, 128 }
 *  { -7936, -137216, -185856, -31616, -256 }
 *  c * { 353792, 1841152, 1304832, 128512, 512 }
 *  { -353792, -9061376, -21253376, -8728576, -518656, -1024}
 *  c * { 22368256, 175627264, 222398464, 56520704, 2084864, 2048 }
 *  // if have long long
 *  { -22368256LL, -795300864LL, -2868264960LL, -2174832640LL, -357888000LL, -8361984LL, -4096LL }
 *  c * { 1903757312LL, 21016670208LL, 41731645440LL, 20261765120LL, 2230947840LL, 33497088LL, 8192LL }
 *  { -1903757312LL, -89702612992LL, -460858269696LL, -559148810240LL, -182172651520LL, -13754155008LL, -134094848LL, -16384LL }
 *  c * { 209865342976LL, 3099269660672LL, 8885192097792LL, 7048869314560LL, 1594922762240LL, 84134068224LL, 536608768LL, 32768LL }
 *  { -209865342976LL, -12655654469632LL, -87815735738368LL, -155964390375424LL, -84842998005760LL, -13684856848384LL, -511780323328LL, -2146926592LL, -65536LL }
 *  c * { 29088885112832LL, 553753414467584LL, 2165206642589696LL, 2550316668551168LL, 985278548541440LL, 115620218667008LL, 3100738912256LL, 8588754944LL, 131072LL }
 *  { -29088885112832LL, -2184860175433728LL, -19686087844429824LL, -48165109676113920LL, -39471306959486976LL, -11124607890751488LL, -965271355195392LL, -18733264797696LL, -34357248000LL, -262144LL }
 *  c * { 4951498053124096LL, 118071834535526400LL, 603968063567560704LL, 990081991141490688LL, 584901762421358592LL, 122829335169859584LL, 7984436548730880LL, 112949304754176LL, 137433710592LL, 524288LL }
 */
template<typename _Tp>
  emsr::Polynomial<_Tp>
  polygamma_poly(unsigned int m)
  {
    emsr::Polynomial<_Tp> a{_Tp{0}, _Tp{1}};
    emsr::Polynomial<_Tp> b{_Tp{1}, _Tp{0}, _Tp{-1}};
    if (m == 0)
      return a;
    else
      {
	auto P = a;
	auto Pp = P.derivative();
	for (unsigned int k = 0; k < m; ++k)
	  {
	    P = -(_Tp(k + 1) * a * P + b * Pp);
	    Pp = P.derivative();
	  }
	return P;
      }
  }

/**
 * 
 */
template<typename _Tp>
  void
  test_polygamma_poly(_Tp proto)
  {
    std::cout.precision(emsr::digits10(proto));
    std::cout << std::showpoint << std::scientific;

    std::cout << '\n' << '\n';
    for (int m = 0; m <= 20; ++m)
      {
	auto P = polygamma_poly<_Tp>(m);
	std::cout << P << '\n';
      }
  }

  /**
   * Return
   * @f[
   *    \frac{\pi^{m+1}}{sin^{m+1}(\pi x)} \frac{d^m}{dx^m} cot(\pi x)
   * @f]
   */
  template<typename _Tp>
    _Tp
    polygamma_reflect(unsigned int m, _Tp x)
    {
      const auto s_pi = emsr::pi_v<_Tp>;
      const auto c = emsr::detail::cos_pi(x);
      const auto cc = c * c;
      const auto s = emsr::detail::sin_pi(x);
      const auto fact = std::pow(s_pi / s, _Tp(m + 1));
      if (m == 0)
	return c * fact
	     * emsr::horner(cc, -1LL);
      else if (m == 1)
	return fact
	     * emsr::horner(cc, 2LL);
      else if (m == 2)
	return c * fact
	     * emsr::horner(cc, -2LL, -4LL);
      else if (m == 3)
	return fact
	     * emsr::horner(cc, 16LL, 8LL);
      else if (m == 4)
	return c * fact
	     * emsr::horner(cc, -16LL, -88LL, -16LL);
      else if (m == 5)
	return fact
	     * emsr::horner(cc, 272LL, 416LL, 32LL);
      else if (m == 6)
	return c * fact
	     * emsr::horner(cc, -272LL, -2880LL, -1824LL, -64LL);
      else
	{
	  auto poly = polygamma_poly<long long>(m);
	  return fact
		* (m & 1 ? poly.eval_odd(c) : poly.eval_even(c));
	}
    }

  /**
   *
   */
  template<typename _Tp>
    _Tp
    polygamma_hurwitz(unsigned int m, _Tp x)
    {
      const auto hzeta = emsr::detail::hurwitz_zeta(_Tp(m + 1), x);
      const auto ln_nfact = emsr::detail::log_gamma(_Tp(m + 1));
      auto result = std::exp(ln_nfact) * hzeta;
      if (m % 2 == 0)
	result = -result;
      return result;
    }

  /**
   * @brief  Return the polygamma function @f$ \psi^{(m)}(x) @f$.
   *
   * The polygamma function is related to the Hurwitz zeta function:
   * @f[
   *   \psi^{(m)}(x) = (-1)^{m+1} m! \zeta(m+1,x)
   * @f]
   */
  template<typename _Tp>
    _Tp
    polygamma(unsigned int m, _Tp x)
    {
      if (x <= _Tp{0})
	{
	  if (const auto n = emsr::fp_is_integer(x); n)
	    return emsr::infinity(x);
	  else
	    return _Tp(m & 1 ? -1 : +1) * polygamma(m, _Tp{1} - x)
		 + polygamma_reflect(m, x);
	}
      else if (m == 0)
	return emsr::detail::digamma(x);
      else
	return polygamma_hurwitz(m, x);
    }

/**
 * 
 */
template<typename _Tp>
  void
  test_polygamma(_Tp proto)
  {
    std::cout.precision(emsr::digits10(proto));
    std::cout << std::showpoint << std::scientific;
    auto w = 8 + std::cout.precision();

    std::cout << '\n' << '\n';
    for (int i = -400; i <= 400; ++i)
      {
	auto x = i * 0.01;
	std::cout << ' ' << std::setw(w) << x;
	for (int m = 0; m <= 10; ++m)
	  {
	    auto P = polygamma(m, x);
	    std::cout << ' ' << std::setw(w) << P;
	  }
	std::cout << '\n';
      }
  }

/**
 * 
 */
template<typename _Tp>
  void
  test_polygamma_reflect_vs_hurwitz(_Tp proto)
  {
    std::cout.precision(emsr::digits10(proto));
    std::cout << std::showpoint << std::scientific;
    auto w = 8 + std::cout.precision();

    std::cout << '\n';
    for (int m = 0; m <= 10; ++m)
      {
	std::cout << '\n';
	std::cout << " m = " << m << '\n';
	for (int i = -40; i <= 0; ++i)
	  {
	    auto x = i * 0.1;
	    auto psi_reflect = polygamma(m, x);
	    auto psi_hurwitz = polygamma_hurwitz(m, x);
	    auto psi_boost = _Tp{0};
	    try
	      {
	        psi_boost = beast::polygamma(m, x);
	      }
	    catch (...)
	      {
	        psi_boost = emsr::quiet_NaN(proto);
	      }
	    std::cout << ' ' << std::setw(w) << x
		      << ' ' << std::setw(w) << psi_reflect
		      << ' ' << std::setw(w) << psi_hurwitz
		      << ' ' << std::setw(w) << psi_boost
		      << ' ' << std::setw(w) << (psi_reflect - psi_hurwitz) / psi_hurwitz
		      << '\n';
	  }
      }
  }


int
main()
{
  test_trigamma(1.0);

  test_tetragamma(1.0);

  test_polygamma_poly(1LL);

  test_polygamma(1.0);

  test_polygamma_reflect_vs_hurwitz(1.0);
}
