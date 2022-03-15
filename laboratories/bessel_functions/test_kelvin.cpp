/**
 *
 */

#include <limits>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <cmath>

#include <emsr/summation.h>
#include <emsr/math_constants.h>
#include <emsr/numeric_limits.h>
#include <emsr/special_functions.h>

namespace emsr
{
namespace detail
{

  /**
   * Data structure for Kelvin functions.
   */
  template<typename _TpNu, typename Tp>
    struct _KelvinState
    {
      _TpNu nu;
      Tp x;
      Tp ber;
      Tp bei;
      Tp ker;
      Tp kei;
    };


  /**
   * Compute the Kelvin functions by series summation.
   *
   * @f[
   *    ber(x) = \sum_{k=0}^\infty \frac{[-(x/4)^4]^k}{[(1/2)_k]^2[(1)_k]^2}
   *           = \sum_{k=0}^\infty \frac{[-(x/2)^4]^k}{[(2k+1)!]^2[(2k)!]^2}
   * @f]
   * @f[
   *    bei(x) = \frac{x^2}{4}
   *             \sum_{k=0}^\infty \frac{[-(x/4)^4]^k}{[(3/2)_k]^2[(1)_k]^2}
   *           = \frac{x^2}{4}
   *             \sum_{k=0}^\infty \frac{[-(x/2)^4]^k}{
   * @f]
   * Note: @f$ (1)_k = k! @f$
   */
  template<typename Tp>
    Tp
    kelvin_bex_series(Tp x, int sign)
    {
      using _WijnSum = emsr::VanWijngaardenSum<Tp>;
      using _WenigerDeltaWijnSum = emsr::WenigerDeltaSum<_WijnSum>;

      const auto s_eps = emsr::epsilon(x);
      constexpr auto s_maxiter = 100;
      const auto y = x / Tp{2};
      const auto y2 = y * y;
      const auto y4 = -y2 * y2;
      auto term = Tp{1};
      _WenigerDeltaWijnSum bex;
      bex += term;
      for (auto k = 1; k < s_maxiter; ++k)
	{
	  const auto fact = Tp{1} / (2 * k * (2 * k + sign));
	  term *= y4 * fact * fact;
	  bex += term;
	  if (std::abs(term) < s_eps * std::abs(bex()))
	    break;
	}
      return bex();
    }


  /**
   * Compute the series sums required for the irregular Kelvin functions.
   * @f[
   *    ker(x): \sum_{k=0}^\infty \frac{[-(x/2)^4]^k}{[(2k)!]^2}H_{2k}
   *    kei(x): \sum_{k=0}^\infty \frac{[-(x/2)^4]^k}{[(2k+1)!]^2}H_{2k+1}
   * @f]
   */
  template<typename Tp>
    Tp
    kelvin_kex_series(Tp x, int sign)
    {
      using _BasicSum = emsr::BasicSum<Tp>;
      using _WijnSum = emsr::VanWijngaardenSum<Tp>;
      using _WenigerDeltaWijnSum = emsr::WenigerDeltaSum<_WijnSum>;

      const auto s_eps = emsr::epsilon(x);
      constexpr auto s_maxiter = 100;
      const auto y = x / Tp{2};
      const auto y2 = y * y;
      const auto y4 = -y2 * y2;
      auto term = Tp{1};
      _BasicSum _H_n;
      _WenigerDeltaWijnSum kex;
      kex += term;
      for (auto k = 1; k < s_maxiter; ++k)
	{
	  const auto tk = 2 * k;
	  const auto tkps = tk + sign;
	  const auto fact = Tp{1} / (tk * tkps);
	  term *= y4 * fact * fact;

	  _H_n += Tp{1} / tk + Tp{1} / tkps;

	  kex += term * _H_n();

	  if (std::abs(term) < s_eps * std::abs(kex()))
	    break;
	}
      return kex();
    }


  /**
   * Return the Kelvin function @f$ ber(x) @f$ for real argument @c x
   * computed by series summation.
   */
  template<typename Tp>
    inline Tp
    kelvin_ber_series(Tp x)
    { return kelvin_bex_series(x, -1); }


  /**
   * Return the Kelvin function @f$ bei(x) @f$ for real argument @c x
   * computed by series summation.
   */
  template<typename Tp>
    inline Tp
    kelvin_bei_series(Tp x)
    { return x * x * kelvin_bex_series(x, +1) / Tp{4}; }


  /**
   * Return the irregular Kelvin function @f$ ker(x) @f$ for real argument @c x
   * computed by series summation.
   *
   * @f[
   *    ker(x) = -\left[log\left(\frac{x}{2}\right) + \gamma_E\right] ber(x)
   *           + \frac{\pi}{4} bei(x)
   *           + \sum_{k=0}^\infty \frac{(-x^4/16)}{[(2k)!]^2}H_{2k}
   * @f]
   */
  template<typename Tp>
    Tp
    kelvin_ker_series(Tp x)
    {
      const auto s_gamma_e = emsr::egamma_v<Tp>;
      const auto s_pi_4 = emsr::pi_v<Tp> / Tp{4};
      const auto ker = kelvin_kex_series(x, -1);
      const auto ber = kelvin_bex_series(x, -1);
      const auto bei = kelvin_bex_series(x, +1);
      const auto x2 = x * x / Tp{4};
      const auto ln = std::log(x / Tp{2}) + s_gamma_e;
      return -ln * ber + s_pi_4 * x2 * bei + ker;
    }


  /**
   * Return the irregular Kelvin function @f$ kei(x) @f$ for real argument @c x
   * computed by series summation.
   *
   * @f[
   *    kei(x) = -\left[log\left(\frac{x}{2}\right) + \gamma_E\right] bei(x)
   *           - \frac{\pi}{4} ber(x)
   *           + \frac{x^2}{4}
   *               \sum_{k=0}^\infty \frac{(-x^4/16)}{[(2k+1)!]^2}H_{2k+1}
   * @f]
   */
  template<typename Tp>
    Tp
    kelvin_kei_series(Tp x)
    {
      const auto s_gamma_e = emsr::egamma_v<Tp>;
      const auto s_pi_4 = emsr::pi_v<Tp> / Tp{4};
      const auto kei = kelvin_kex_series(x, +1);
      const auto ber = kelvin_bex_series(x, -1);
      const auto bei = kelvin_bex_series(x, +1);
      const auto x2 = x * x / Tp{4};
      const auto ln = std::log(x / Tp{2}) + s_gamma_e;
      return -ln * x2 * bei - s_pi_4 * ber + x2 * kei;
    }


  /**
   * Compute the Kelvin functions by series expansion.
   * Computing all four functions is probably not that useful because of the
   * different stability breaks.
   *
   * The Kelvin functions and their expansions are:
   * @f[
   *    ber(x) = \sum_{k=0}^\infty  \frac{(-x^4/16)}{[(2k)!]^2}
   * @f]
   * @f[
   *    bei(x) = \frac{x^2}{4}\sum_{k=0}^\infty \frac{(-x^4/16)}{[(2k+1)!]^2}
   * @f]
   * @f[
   *    ker(x) = -\left[log\left(\frac{x}{2}\right) + \gamma_E\right] ber(x)
   *           + \frac{\pi}{4} bei(x)
   *           + \sum_{k=0}^\infty \frac{(-x^4/16)}{[(2k)!]^2}H_{2k}
   * @f]
   * @f[
   *    kei(x) = -\left[log\left(\frac{x}{2}\right) + \gamma_E\right] bei(x)
   *           - \frac{\pi}{4} ber(x)
   *           + \frac{x^2}{4}
   *               \sum_{k=0}^\infty \frac{(-x^4/16)}{[(2k+1)!]^2}H_{2k+1}
   * @f]
   * where
   * @f[
   *    H_n = \sum_{k=1}^n \frac{1}{k}
   * @f]
   * is the harmonic number.
   */
  template<typename Tp>
    _KelvinState<int, Tp>
    kelvin_series(Tp x)
    {
      using _BasicSum = emsr::BasicSum<Tp>;
      using _WijnSum = emsr::VanWijngaardenSum<Tp>;
      using _WenigerDeltaWijnSum = emsr::WenigerDeltaSum<_WijnSum>;

      constexpr auto s_maxiter = 100;
      const auto s_gamma_e = emsr::egamma_v<Tp>;
      const auto s_pi_4 = emsr::pi_v<Tp> / Tp{4};
      const auto s_eps = emsr::epsilon(x);
      const auto y = x / Tp{2};
      const auto x2 = y * y;
      const auto y4 = x2 * x2;
      if (x == Tp{0})
	{
	  const auto s_inf = emsr::infinity<Tp>();
	  return {0, x, Tp{1}, Tp{0}, s_inf, -s_pi_4};
	}
      else
	{
	  auto term = Tp{1};
	  _BasicSum _H_n;
	  _WenigerDeltaWijnSum ber, bei, ker, kei;
	  ber += term;
	  ker += term;
	  bei += term;
	  kei += term;
	  _H_n += Tp{1};
	  for (auto k = 1; k < s_maxiter; ++k)
	    {
	      const auto tk = Tp(2 * k);
	      const auto tkp1 = Tp(2 * k + 1);

	      const auto factr = Tp{1} / tk;
	      term *= -y4 * factr * factr;
	      ber += term;

	      _H_n += Tp{1} / tk;
	      ker += term * _H_n();

	      const auto facti = Tp{1} / tkp1;
	      term *= facti * facti;
	      bei += term;

	      _H_n += Tp{1} / tkp1;
	      kei += term * _H_n();

	      if (std::abs(term) < s_eps * std::abs(ber())
	       && std::abs(term) < s_eps * std::abs(bei()))
		break;
	    }
	  const auto ln = std::log(x / Tp{2}) + s_gamma_e;
	  return _KelvinState<int, Tp>{0, x, ber(), x2 * bei(),
	      -ln * ber() + s_pi_4 * x2 * bei() + ker(),
	      -ln * x2 * bei() - s_pi_4 * ber() + x2 * kei()};
	}
    }


  /**
   * Compute the Kelvin functions by asymptotic series expansion.
   *
   *
   *
   *
   */
  template<typename Tp>
    _KelvinState<int, Tp>
    kelvin_asymp(Tp x)
    {
      using _Cmplx = std::complex<Tp>;
      using _BasicSum = emsr::BasicSum<_Cmplx>;
      using _WenigerDeltaSum = emsr::WenigerDeltaSum<_BasicSum>;

      const auto s_j = _Cmplx{0, 1};
      const auto s_1d2 = Tp{1} / Tp{2};
      const auto s_pi = emsr::pi_v<Tp>;
      const auto s_pi_4 = emsr::pi_v<Tp> / Tp{4};
      const auto s_pi_8 = s_1d2 * s_pi_4;
      const auto s_3pi_4 = Tp{3} * s_pi_4;
      const auto s_sqrt_2 = emsr::sqrt2_v<Tp>;
      const auto s_sqrt_pi = emsr::sqrtpi_v<Tp>;
      const auto s_eps = emsr::epsilon(x);
      constexpr auto s_maxiter = 1000;
      const auto y = Tp{1} / (Tp{8} * x);
      auto term = Tp{1};
      const auto xrt2 = x / s_sqrt_2;
      auto barg = xrt2 - s_pi_8;
      auto karg = -xrt2 - s_pi_8;
      _WenigerDeltaSum be, ke;
      be += std::polar(term, barg);
      ke += std::polar(term, karg);
      for (auto k = 1; k < s_maxiter; ++k)
	{
	  barg -= s_pi_4;
	  karg += s_3pi_4;
	  auto fact = Tp(2 * k - 1);
	  auto next = y * fact * fact / Tp(k);
	  if (std::abs(next) > Tp{1})
	    break;
	  term *= -next;
	  be += std::polar(term, barg);
	  ke += std::polar(term, karg);

	  if (std::abs(term) < s_eps * std::abs(be()))
	    break;
	}
      const auto exp = std::exp(xrt2);
      const auto rt = std::sqrt(Tp{2} * x);
      const auto kfact = s_sqrt_pi / rt / exp;
      const auto kex = kfact * ke();
      const auto bfact = exp / s_sqrt_pi / rt;
      const auto bex = bfact * be() + s_j * kex / s_pi;
      return _KelvinState<int, Tp>{0, x,
				    std::real(bex), std::imag(bex),
				    std::real(kex), std::imag(kex)};
    }

  /**
   * Compute the Kelvin functions by asymptotic series expansion.
   */
  template<typename Tp>
    inline Tp
    kelvin_ber_asymp(Tp x)
    { return kelvin_asymp(x).ber; }

  /**
   * Compute the Kelvin functions by asymptotic series expansion.
   */
  template<typename Tp>
    inline Tp
    kelvin_bei_asymp(Tp x)
    { return kelvin_asymp(x).bei; }

  /**
   * Compute the Kelvin functions by asymptotic series expansion.
   */
  template<typename Tp>
    inline Tp
    kelvin_ker_asymp(Tp x)
    { return kelvin_asymp(x).ker; }

  /**
   * Compute the Kelvin functions by asymptotic series expansion.
   */
  template<typename Tp>
    inline Tp
    kelvin_kei_asymp(Tp x)
    { return kelvin_asymp(x).kei; }


  /**
   * Compute the Kelvin functions of integer order n by series expansion.
   */
  template<typename Tp>
    _KelvinState<int, Tp>
    kelvin_series(int n, Tp x)
    {
      using _Cmplx = std::complex<Tp>;
      using _BasicSum = emsr::BasicSum<_Cmplx>;
      //using _WenigerDeltaSum = emsr::WenigerDeltaSum<_BasicSum>;

      const auto s_j = _Cmplx{0, 1};
      const auto s_pi_2 = emsr::pi_v<Tp> / Tp{2};
      const auto s_pi_4 = emsr::pi_v<Tp> / Tp{4};
      const auto s_3pi_4 = Tp{3} * s_pi_4;
      const auto s_eps = emsr::epsilon(x);
      constexpr auto s_maxiter = 1000;
      if (n < 0)
	{
	  const auto _Cnp = Tp(1 - 2 * (n & 1));
	  auto _Kv = kelvin_series(-n, x);
	  return _KelvinState<int, Tp>{n, x,
	  				_Cnp * _Kv.ber, _Cnp * _Kv.bei,
					_Cnp * _Kv.ker, _Cnp * _Kv.kei};
	}
      else
	{
	  const auto y = x / Tp{2};
	  const auto y2 = y * y;

	  _BasicSum be;
	  auto bterm = Tp{1};
	  auto barg = s_3pi_4 * Tp(n);
	  be += std::polar(bterm, barg);
	  for (auto k = 1; k < s_maxiter; ++k)
	    {
	      bterm *= y2 / Tp(k + n) / Tp(k);
	      barg += s_pi_2;
	      be += std::polar(bterm, barg);
	      if (std::abs(bterm) < s_eps * std::abs(be()))
		break;
	    }

	  _BasicSum ke1;
	  if (n > 0)
	    {
	      auto kterm1 = Tp{1};
	      auto karg1 = s_3pi_4 * Tp(n);
	      ke1 += std::polar(kterm1, -karg1);
	      for (auto k = 1; k < n - 1; ++k)
		{
		  kterm1 *= y2 / Tp(n - 1 - k) / Tp(k);
		  karg1 += s_pi_2;
		  ke1 += std::polar(kterm1, -karg1);
		}
	      if (n > 1)
		{
		  kterm1 *= y2 / Tp(n - 1);
		  karg1 += s_pi_2;
		  ke1 += std::polar(kterm1, -karg1);
		}
	    }

	  _BasicSum ke2;
	  auto hsum2 = digamma<Tp>(1) + digamma<Tp>(1 + n);
	  auto kterm2 = Tp{1};
	  auto karg2 = s_3pi_4 * Tp(n);
	  ke2 += std::polar(hsum2 * kterm2, karg2);
	  for (auto k = 1; k < s_maxiter; ++k)
	    {
	      hsum2 += Tp{1} / Tp(1 + k)
		       + Tp{1} / Tp(1 + n + k);
	      kterm2 *= y2 / Tp(k) / Tp(n + k);
	      karg2 += s_pi_2;
	      ke2 += std::polar(hsum2 * kterm2, karg2);
	      if (std::abs(hsum2 * kterm2) < s_eps * std::abs(ke2()))
		break;
	    }

	  const auto nfact = factorial<Tp>(n);
	  const auto pow = std::pow(y, Tp(n));
	  const auto bex = pow * be() / nfact;
	  const auto fact = n > 0 ? factorial<Tp>(n - 1) : Tp{1};
	  const auto kex = -std::log(y) * bex - s_j * s_pi_4 * bex
			   + ke1() / pow / Tp{2} / fact
			   + pow * ke2() / nfact / Tp{2};
	  return _KelvinState<int, Tp>{n, x,
				     std::real(bex), std::imag(bex),
				     std::real(kex), std::imag(kex)};
	}
    }


  /**
   * Compute the Kelvin functions of order @f$ \nu @f$ by series expansion.
   */
  template<typename Tp>
    _KelvinState<Tp, Tp>
    kelvin_series(Tp nu, Tp x)
    {
      using _Cmplx = std::complex<Tp>;
      using _BasicSum = emsr::BasicSum<_Cmplx>;
      using _WenigerDeltaSum = emsr::WenigerDeltaSum<_BasicSum>;

      const auto s_j = _Cmplx{0, 1};
      const auto s_pi = emsr::pi_v<Tp>;
      const auto s_pi_2 = emsr::pi_v<Tp> / Tp{2};
      const auto s_pi_4 = emsr::pi_v<Tp> / Tp{4};
      const auto s_3pi_4 = Tp{3} * s_pi_4;
      constexpr auto s_maxiter = 1000;

      if (auto nuint = emsr::fp_is_integer(nu); nuint)
	{
	  const auto ret = kelvin_series(nuint(), x);
	  return _KelvinState<Tp, Tp>{Tp(ret.nu), ret.x,
					ret.ber, ret.bei,
					ret.ker, ret.kei};
	}
      else
	{
	  const auto y = x / Tp{2};
	  const auto y2 = y * y;

	  _WenigerDeltaSum bep;
	  auto termp = Tp{1};
	  auto argp = s_3pi_4 * nu;
	  bep += std::polar(termp, argp);
	  for (auto k = 1; k < s_maxiter; ++k)
	    {
	      termp *= y2 / Tp(k + nu) / Tp(k);
	      argp += s_pi_2;
	      bep += std::polar(termp, argp);
	    }

	  const auto bex = std::pow(y, nu) * bep()
			   / std::tgamma(Tp{1} + nu);

	  _WenigerDeltaSum bem;
	  auto termm = Tp{1};
	  auto argm = -s_3pi_4 * nu;
	  bem += std::polar(termm, argm);
	  for (auto k = 1; k < s_maxiter; ++k)
	    {
	      termm *= y2 / Tp(k - nu) / Tp(k);
	      argm += s_pi_2;
	      bem += std::polar(termm, argm);
	    }

	  const auto bey = std::pow(y, -nu) * bem()
			   / std::tgamma(Tp{1} - nu);

	  const auto nupi = nu * s_pi;
	  const auto csc = Tp{1} / std::sin(nupi);
	  const auto cot = Tp{1} / std::tan(nupi);
	  const auto kex = s_pi_2
			   * (csc * bey - cot * bex + s_j * bex);

	  return _KelvinState<Tp, Tp>{nu, x,
					std::real(bex), std::imag(bex),
					std::real(kex), std::imag(kex)};
	}
    }


  /**
   * Compute the Kelvin functions by asymptotic series expansion.
   */
  template<typename Tp>
    _KelvinState<Tp, Tp>
    kelvin_asymp(Tp nu, Tp x)
    {
      using _Cmplx = std::complex<Tp>;
      using _BasicSum = emsr::BasicSum<_Cmplx>;
      using _WenigerDeltaSum = emsr::WenigerDeltaSum<_BasicSum>;

      constexpr auto s_j = _Cmplx{0, 1};
      constexpr auto s_1d2 = Tp{1} / Tp{2};
      const auto s_pi = emsr::pi_v<Tp>;
      const auto s_pi_2 = emsr::pi_v<Tp> / Tp{2};
      const auto s_pi_4 = emsr::pi_v<Tp> / Tp{4};
      const auto s_pi_8 = s_1d2 * s_pi_4;
      const auto s_sqrt_2 = emsr::sqrt2_v<Tp>;
      const auto s_sqrt_pi = emsr::sqrtpi_v<Tp>;
      const auto s_eps = emsr::epsilon(x);
      constexpr auto s_maxiter = 1000;
      const auto y = Tp{1} / (Tp{2} * x);
      auto bterm = Tp{1};
      auto kterm = Tp{1};
      const auto xrt2 = x / s_sqrt_2;
      auto barg = xrt2 + nu * s_pi_2 - s_pi_8;
      auto karg = xrt2 + nu * s_pi_2 + s_pi_8;
      _WenigerDeltaSum be, ke;
      be += std::polar(bterm, barg);
      ke += std::polar(kterm, karg);
      for (auto k = 1; k < s_maxiter; ++k)
	{
	  barg -= s_pi_4;
	  karg += s_pi_4;
	  auto fact = (Tp(k) - s_1d2 - nu)
		      * (Tp(k) - s_1d2 + nu) / Tp(k);
	  auto next = y * fact;
	  if (std::abs(next) > Tp{1})
	    break;

	  bterm *= next;
	  kterm *= -bterm;
	  be += std::polar(bterm, barg);
	  ke += std::polar(kterm, karg);

	  if (std::abs(bterm) < s_eps * std::abs(be()))
	    break;
	}
      const auto exp = std::exp(xrt2);
      const auto rt = s_sqrt_2 * std::sqrt(x);
      const auto kfact = s_sqrt_pi / rt / exp;
      const auto kex = kfact * ke();
      const auto bfact = exp / s_sqrt_pi / rt;
      const auto bex = bfact * be() + s_j * kex / s_pi;
      return _KelvinState<Tp, Tp>{nu, x,
				std::real(bex), std::imag(bex),
				std::real(kex), -std::imag(kex)};
    }

  /**
   * Compute the Kelvin function @f$ ber(x)@f$.
   */
  template<typename Tp>
    inline Tp
    kelvin_ber(Tp x)
    {
      constexpr auto s_switch = Tp{26};
      if (std::isnan(x))
	return emsr::quiet_NaN<Tp>();
      else if (std::abs(x) < s_switch)
	return kelvin_ber_series(x);
      else
	return kelvin_ber_asymp(x);
    }

  /**
   * Compute the Kelvin function @f$ bei(x)@f$.
   */
  template<typename Tp>
    inline Tp
    kelvin_bei(Tp x)
    {
      constexpr auto s_switch = Tp{26};
      if (std::isnan(x))
	return emsr::quiet_NaN<Tp>();
      else if (std::abs(x) < s_switch)
	return kelvin_bei_series(x);
      else
	return kelvin_bei_asymp(x);
    }

  /**
   * Compute the Kelvin function @f$ ker(x)@f$.
   */
  template<typename Tp>
    inline Tp
    kelvin_ker(Tp x)
    {
      constexpr auto s_switch = Tp{5};
      if (std::isnan(x))
	return emsr::quiet_NaN<Tp>();
      else if (std::abs(x) < s_switch)
	return kelvin_ker_series(x);
      else
	return kelvin_ker_asymp(x);
    }

  /**
   * Compute the Kelvin function @f$ kei(x)@f$.
   */
  template<typename Tp>
    inline Tp
    kelvin_kei(Tp x)
    {
      constexpr auto s_switch = Tp{5};
      if (std::isnan(x))
	return emsr::quiet_NaN<Tp>();
      else if (std::abs(x) < s_switch)
	return kelvin_kei_series(x);
      else
	return kelvin_kei_asymp(x);
    }

  /**
   * Compute the Kelvin function @f$ ber_\nu(x)@f$.
   */
  template<typename Tp>
    inline Tp
    kelvin_ber(Tp nu, Tp x)
    {
      constexpr auto s_switch = Tp{26};
      if (std::isnan(nu) || std::isnan(x))
	return emsr::quiet_NaN<Tp>();
      else if (std::abs(x) < s_switch)
	return kelvin_series(nu, x).ber;
      else
	return kelvin_asymp(nu, x).ber;
    }

  /**
   * Compute the Kelvin function @f$ bei_\nu(x)@f$.
   */
  template<typename Tp>
    inline Tp
    kelvin_bei(Tp nu, Tp x)
    {
      constexpr auto s_switch = Tp{26};
      if (std::isnan(nu) || std::isnan(x))
	return emsr::quiet_NaN<Tp>();
      else if (std::abs(x) < s_switch)
	return kelvin_series(nu, x).bei;
      else
	return kelvin_asymp(nu, x).bei;
    }

  /**
   * Compute the Kelvin function @f$ ker_\nu(x)@f$.
   */
  template<typename Tp>
    inline Tp
    kelvin_ker(Tp nu, Tp x)
    {
      constexpr auto s_switch = Tp{5};
      if (std::isnan(nu) || std::isnan(x))
	return emsr::quiet_NaN<Tp>();
      else if (std::abs(x) < s_switch)
	return kelvin_series(nu, x).ker;
      else
	return kelvin_asymp(nu, x).ker;
    }

  /**
   * Compute the Kelvin function @f$ kei_\nu(x)@f$.
   */
  template<typename Tp>
    inline Tp
    kelvin_kei(Tp nu, Tp x)
    {
      const auto s_switch = Tp{5};
      if (std::isnan(nu) || std::isnan(x))
	return emsr::quiet_NaN<Tp>();
      else if (std::abs(x) < s_switch)
	return kelvin_series(nu, x).kei;
      else
	return kelvin_asymp(nu, x).kei;
    }

} // namespace detail
} // namespace emsr

namespace emsr
{


  //  Kelvin functions

  /**
   * Return the Kelvin function @f$ ber(x) @f$ for @c float argument @c x.
   *
   * @see kelvin_ber for details.
   */
  inline float
  kelvin_berf(float x)
  { return emsr::detail::kelvin_ber<float>(x); }

  /**
   * Return the Kelvin function @f$ ber(x) @f$
   * for <tt>long double</tt> argument @c x.
   *
   * @see kelvin_ber for details.
   */
  inline long double
  kelvin_berl(long double x)
  { return emsr::detail::kelvin_ber<long double>(x); }

  /**
   * Return the Kelvin function @f$ ber(x) @f$ for @c real argument @c x.
   *
   * Kelvin integrals @f$ ber(x) @f$ and @f$ bei(x) @f$ are given by
   * \f[
   *   J_0(\frac{x-ix}{\sqrt{2}}) = ber(x) + i bei(x)
   * \f]
   * where @f$ J_0(x) @f$ is the Bessel function of the first kind.
   *
   * @tparam Tp The floating-point type of the argument @c x.
   * @param  x  The argument of the Kelvin function
   */
  template<typename Tp>
    inline typename emsr::fp_promote_t<Tp>
    kelvin_ber(Tp x)
    {
      typedef typename emsr::fp_promote_t<Tp> type;
      return emsr::detail::kelvin_ber<type>(x);
    }

  /**
   * Return the Kelvin function @f$ bei(x) @f$ for @c float argument @c x.
   *
   * @see kelvin_bei for details.
   */
  inline float
  kelvin_beif(float x)
  { return emsr::detail::kelvin_bei<float>(x); }

  /**
   * Return the Kelvin function @f$ bei(x) @f$
   * for <tt>long double</tt> argument @c x.
   *
   * @see kelvin_bei for details.
   */
  inline long double
  kelvin_beil(long double x)
  { return emsr::detail::kelvin_bei<long double>(x); }

  /**
   * Return the Kelvin function @f$ bei(x) @f$ for @c real argument @c x.
   *
   * Kelvin integrals @f$ bei(x) @f$ and @f$ bei(x) @f$ are given by
   * \f[
   *   J_0(\frac{x-ix}{\sqrt{2}}) = ber(x) + i bei(x)
   * \f]
   * where @f$ J_0(x) @f$ is the Bessel function of the first kind.
   *
   * @tparam Tp The floating-point type of the argument @c x.
   * @param  x  The argument of the Kelvin function
   */
  template<typename Tp>
    inline typename emsr::fp_promote_t<Tp>
    kelvin_bei(Tp x)
    {
      typedef typename emsr::fp_promote_t<Tp> type;
      return emsr::detail::kelvin_bei<type>(x);
    }

  /**
   * Return the Kelvin function @f$ ker(x) @f$ for @c float argument @c x.
   *
   * @see kelvin_ker for details.
   */
  inline float
  kelvin_kerf(float x)
  { return emsr::detail::kelvin_ker<float>(x); }

  /**
   * Return the Kelvin function @f$ ker(x) @f$
   * for <tt>long double</tt> argument @c x.
   *
   * @see kelvin_ker for details.
   */
  inline long double
  kelvin_kerl(long double x)
  { return emsr::detail::kelvin_ker<long double>(x); }

  /**
   * Return the Kelvin function @f$ ker(x) @f$ for @c real argument @c x.
   *
   * Kelvin integrals @f$ ker(x) @f$ and @f$ kei(x) @f$ are given by
   * \f[
   *   J_0(\frac{x-ix}{\sqrt{2}}) = ker(x) + i kei(x)
   * \f]
   * where @f$ J_0(x) @f$ is the Bessel function of the first kind.
   *
   * @tparam Tp The floating-point type of the argument @c x.
   * @param  x  The argument of the Kelvin function
   */
  template<typename Tp>
    inline typename emsr::fp_promote_t<Tp>
    kelvin_ker(Tp x)
    {
      typedef typename emsr::fp_promote_t<Tp> type;
      return emsr::detail::kelvin_ker<type>(x);
    }

  /**
   * Return the Kelvin function @f$ kei(x) @f$ for @c float argument @c x.
   *
   * @see kelvin_kei for details.
   */
  inline float
  kelvin_keif(float x)
  { return emsr::detail::kelvin_kei<float>(x); }

  /**
   * Return the Kelvin function @f$ kei(x) @f$
   * for <tt>long double</tt> argument @c x.
   *
   * @see kelvin_kei for details.
   */
  inline long double
  kelvin_keil(long double x)
  { return emsr::detail::kelvin_kei<long double>(x); }

  /**
   * Return the Kelvin function @f$ kei(x) @f$ for @c real argument @c x.
   *
   * Kelvin integrals @f$ kei(x) @f$ and @f$ kei(x) @f$ are given by
   * \f[
   *   J_0(\frac{x-ix}{\sqrt{2}}) = ker(x) + i kei(x)
   * \f]
   * where @f$ J_0(x) @f$ is the Bessel function of the first kind.
   *
   * @tparam Tp The floating-point type of the argument @c x.
   * @param  x  The argument of the Kelvin function
   */
  template<typename Tp>
    inline typename emsr::fp_promote_t<Tp>
    kelvin_kei(Tp x)
    {
      typedef typename emsr::fp_promote_t<Tp> type;
      return emsr::detail::kelvin_kei<type>(x);
    }

  /**
   * Return the Kelvin function @f$ ber_\nu(x) @f$ for float
   * order @f$ \nu @f$ and argument @c x.
   *
   * @see kelvin_ber for details.
   */
  inline float
  kelvin_berf(float nu, float x)
  { return emsr::detail::kelvin_ber<float>(nu, x); }

  /**
   * Return the Kelvin function @f$ ber_\nu(x) @f$ for <tt>long double</tt>
   * order @f$ \nu @f$ and argument @c x.
   *
   * @see kelvin_ber for details.
   */
  inline long double
  kelvin_berl(long double nu, long double x)
  { return emsr::detail::kelvin_ber<long double>(nu, x); }

  /**
   * Return the Kelvin function @f$ ber_\nu(x) @f$ of real
   * order @f$ \nu @f$ and argument @f$ x @f$.
   * The Kelvin functions @f$ ber_\nu(x) @f$ and @f$ bei_\nu(x) @f$
   * are defined by:
   * @f[
   *    J_\nu(\frac{x-ix}{\sqrt{2}}) = J_\nu(xe^{i3\pi/4})
   *     = ber_\nu(x) + i bei_\nu(x)
   * @f]
   * where J_\nu(x) is the Bessel function of the first kind.
   *
   * @tparam Tp The real type of the argument
   * @param nu The real order
   * @param x The real argument
   */
  template<typename Tp>
    inline emsr::fp_promote_t<Tp>
    kelvin_ber(Tp nu, Tp x)
    {
      using type = emsr::fp_promote_t<Tp>;
      return emsr::detail::kelvin_ber<type>(nu, x);
    }

  /**
   * Return the Kelvin function @f$ bei_\nu(x) @f$ for float
   * order @f$ \nu @f$ and argument @c x.
   *
   * @see kelvin_bei for details.
   */
  inline float
  kelvin_beif(float nu, float x)
  { return emsr::detail::kelvin_bei<float>(nu, x); }

  /**
   * Return the Kelvin function @f$ bei_\nu(x) @f$ for <tt>long double</tt>
   * order @f$ \nu @f$ and argument @c x.
   *
   * @see kelvin_bei for details.
   */
  inline long double
  kelvin_beil(long double nu, long double x)
  { return emsr::detail::kelvin_bei<long double>(nu, x); }

  /**
   * Return the Kelvin function @f$ bei_\nu(x) @f$ of real
   * order @f$ \nu @f$ and argument @f$ x @f$.
   * The Kelvin functions @f$ bei_\nu(x) @f$ and @f$ bei_\nu(x) @f$
   * are defined by:
   * @f[
   *    J_\nu(\frac{x-ix}{\sqrt{2}}) = J_\nu(xe^{i3\pi/4})
   *     = ber_\nu(x) + i bei_\nu(x)
   * @f]
   * where J_\nu(x) is the Bessel function of the first kind.
   *
   * @tparam Tp The real type of the argument
   * @param nu The real order
   * @param x The real argument
   */
  template<typename Tp>
    inline emsr::fp_promote_t<Tp>
    kelvin_bei(Tp nu, Tp x)
    {
      using type = emsr::fp_promote_t<Tp>;
      return emsr::detail::kelvin_bei<type>(nu, x);
    }

  /**
   * Return the Kelvin function @f$ ker_\nu(x) @f$ for float
   * order @f$ \nu @f$ and argument @c x.
   *
   * @see kelvin_ker for details.
   */
  inline float
  kelvin_kerf(float nu, float x)
  { return emsr::detail::kelvin_ker<float>(nu, x); }

  /**
   * Return the Kelvin function @f$ ker_\nu(x) @f$ for <tt>long double</tt>
   * order @f$ \nu @f$ and argument @c x.
   *
   * @see kelvin_ker for details.
   */
  inline long double
  kelvin_kerl(long double nu, long double x)
  { return emsr::detail::kelvin_ker<long double>(nu, x); }

  /**
   * Return the Kelvin function @f$ ker_\nu(x) @f$ of real
   * order @f$ \nu @f$ and argument @f$ x @f$.
   * The Kelvin functions @f$ ker_\nu(x) @f$ and @f$ kei_\nu(x) @f$
   * are defined by:
   * @f[
   *    K_\nu(\frac{x+ix}{\sqrt{2}}) = K_\nu(xe^{i\pi/4})
   *     = ker_\nu(x) + i kei_\nu(x)
   * @f]
   * where K_\nu(x) is the regular modified Bessel function.
   *
   * @tparam Tp The real type of the argument
   * @param nu The real order
   * @param x The real argument
   */
  template<typename Tp>
    inline emsr::fp_promote_t<Tp>
    kelvin_ker(Tp nu, Tp x)
    {
      using type = emsr::fp_promote_t<Tp>;
      return emsr::detail::kelvin_ker<type>(nu, x);
    }

  /**
   * Return the Kelvin function @f$ kei_\nu(x) @f$ for float
   * order @f$ \nu @f$ and argument @c x.
   *
   * @see kelvin_kei for details.
   */
  inline float
  kelvin_keif(float nu, float x)
  { return emsr::detail::kelvin_kei<float>(nu, x); }

  /**
   * Return the Kelvin function @f$ kei_\nu(x) @f$ for <tt>long double</tt>
   * order @f$ \nu @f$ and argument @c x.
   *
   * @see kelvin_kei for details.
   */
  inline long double
  kelvin_keil(long double nu, long double x)
  { return emsr::detail::kelvin_kei<long double>(nu, x); }

  /**
   * Return the Kelvin function @f$ kei_\nu(x) @f$ of real
   * order @f$ \nu @f$ and argument @f$ x @f$.
   * The Kelvin functions @f$ kei_\nu(x) @f$ and @f$ kei_\nu(x) @f$
   * are defined by:
   * @f[
   *    K_\nu(\frac{x+ix}{\sqrt{2}}) = K_\nu(xe^{i\pi/4})
   *     = ker_\nu(x) + i kei_\nu(x)
   * @f]
   * where K_\nu(x) is the regular modified Bessel function.
   *
   * @tparam Tp The real type of the argument
   * @param nu The real order
   * @param x The real argument
   */
  template<typename Tp>
    inline emsr::fp_promote_t<Tp>
    kelvin_kei(Tp nu, Tp x)
    {
      using type = emsr::fp_promote_t<Tp>;
      return emsr::detail::kelvin_kei<type>(nu, x);
    }

} // namespace emsr


/**
 * Run Kelvin with individual sums: kelvin_ber_series, etc.
 */
template<typename Tp>
  void
  run_kelvin1(Tp proto = Tp{})
  {
    std::cout.precision(emsr::digits10(proto));
    std::cout << std::showpoint << std::scientific;
    auto w = 8 + std::cout.precision();

    emsr::detail::kelvin_ber_series(Tp{});

    std::cout << "\n\nPrint Kelvin functions computed by series expansions\n";
    std::cout << std::setw(w) << "x"
	      << std::setw(w) << "ber"
	      << std::setw(w) << "bei"
	      << std::setw(w) << "ker"
	      << std::setw(w) << "kei"
	      << '\n';
    const auto del = Tp{1}/Tp{10};
    for (int i = 0; i <= 200; ++i)
      {
	auto x = del * i;
	auto ber = emsr::detail::kelvin_ber_series(x);
	auto bei = emsr::detail::kelvin_bei_series(x);
	auto ker = emsr::detail::kelvin_ker_series(x);
	auto kei = emsr::detail::kelvin_kei_series(x);
	std::cout << std::setw(w) << x
		  << std::setw(w) << ber
		  << std::setw(w) << bei
		  << std::setw(w) << ker
		  << std::setw(w) << kei
		  << '\n';
      }
    std::cout << '\n' << std::flush;
  }


/**
 * Run Kelvin with kelvin_series which combines sums.
 */
template<typename Tp>
  void
  run_kelvin2(Tp proto = Tp{})
  {
    std::cout.precision(emsr::digits10(proto));
    std::cout << std::showpoint << std::scientific;
    auto w = 8 + std::cout.precision();

    emsr::detail::kelvin_series(Tp{});

    std::cout << "\n\nPrint Kelvin functions computed by series expansions\n";
    std::cout << std::setw(w) << "x"
	      << std::setw(w) << "ber"
	      << std::setw(w) << "bei"
	      << std::setw(w) << "ker"
	      << std::setw(w) << "kei"
	      << '\n';
    const auto del = Tp{1}/Tp{10};
    for (int i = 0; i <= 200; ++i)
      {
	auto x = del * i;
	auto ke = emsr::detail::kelvin_series(x);
	std::cout << std::setw(w) << ke.x
		  << std::setw(w) << ke.ber
		  << std::setw(w) << ke.bei
		  << std::setw(w) << ke.ker
		  << std::setw(w) << ke.kei
		  << '\n';
      }
    std::cout << '\n' << std::flush;
  }


/**
 * Diff Kelvin kelvin_asymp with kelvin_series.
 */
template<typename Tp>
  void
  diff_kelvin2(Tp proto = Tp{})
  {
    std::cout.precision(emsr::digits10(proto));
    std::cout << std::showpoint << std::scientific;
    auto w = 8 + std::cout.precision();

    std::cout << "\n\nRelative Difference Kelvin functions computed by series and asymptotic expansions\n";
    std::cout << ' ' << std::setw(w) << "x"
	      << ' ' << std::setw(w) << "ber asymp"
	      << ' ' << std::setw(w) << "ber series"
	      << ' ' << std::setw(w) << "bei asymp"
	      << ' ' << std::setw(w) << "bei series"
	      << ' ' << std::setw(w) << "ker asymp"
	      << ' ' << std::setw(w) << "ker series"
	      << ' ' << std::setw(w) << "kei asymp"
	      << ' ' << std::setw(w) << "kei series"
	      << ' ' << std::setw(w) << "ber asy-ser"
	      << ' ' << std::setw(w) << "bei asy-ser"
	      << ' ' << std::setw(w) << "ker asy-ser"
	      << ' ' << std::setw(w) << "kei asy-ser"
	      << '\n';
    const auto del = Tp{1}/Tp{10};
    for (int i = 50; i <= 400; ++i)
      {
	auto x = del * i;
	auto kes = emsr::detail::kelvin_series(x);
	auto kea = emsr::detail::kelvin_asymp(x);
	std::cout << ' ' << std::setw(w) << kes.x
		  << ' ' << std::setw(w) << kea.ber
		  << ' ' << std::setw(w) << kes.ber
		  << ' ' << std::setw(w) << kea.bei
		  << ' ' << std::setw(w) << kes.bei
		  << ' ' << std::setw(w) << kea.ker
		  << ' ' << std::setw(w) << kes.ker
		  << ' ' << std::setw(w) << kea.kei
		  << ' ' << std::setw(w) << kes.kei
		  << ' ' << std::setw(w) << (kea.ber - kes.ber) / std::abs(kes.ber)
		  << ' ' << std::setw(w) << (kea.bei - kes.bei) / std::abs(kes.bei)
		  << ' ' << std::setw(w) << (kea.ker - kes.ker) / std::abs(kes.ker)
		  << ' ' << std::setw(w) << (kea.kei - kes.kei) / std::abs(kes.kei)
		  << '\n';
      }
    std::cout << '\n' << std::flush;
  }


/**
 * Run Kelvin with kelvin_series(nu, x) which combines sums for general order.
 */
template<typename Tp>
  void
  run_kelvin3(Tp nu = Tp{0})
  {
    std::cout.precision(emsr::digits10(nu));
    std::cout << std::showpoint << std::scientific;
    auto w = 8 + std::cout.precision();

    std::cout << "\n\nPrint Kelvin functions computed by series expansions, nu = " << nu << '\n';
    std::cout << std::setw(w) << "x"
	      << std::setw(w) << "ber asy-ser"
	      << std::setw(w) << "bei asy-ser"
	      << std::setw(w) << "ker asy-ser"
	      << std::setw(w) << "kei asy-ser"
	      << '\n';
    const auto del = Tp{1}/Tp{10};
    for (int i = 0; i <= 200; ++i)
      {
	auto x = del * i;
	auto ke = emsr::detail::kelvin_series(nu, x);
	std::cout << std::setw(w) << ke.x
		  << std::setw(w) << ke.ber
		  << std::setw(w) << ke.bei
		  << std::setw(w) << ke.ker
		  << std::setw(w) << ke.kei
		  << '\n';
      }
    std::cout << '\n' << std::flush;
  }


/**
 * Diff Kelvin kelvin_asymp(nu, x) with kelvin_series(nu, x) for general order.
 */
template<typename Tp>
  void
  diff_kelvin3(Tp nu = Tp{0})
  {
    std::cout.precision(emsr::digits10(nu));
    std::cout << std::showpoint << std::scientific;
    auto w = 8 + std::cout.precision();

    std::cout << "\n\nDiff Kelvin functions computed by series expansions, nu = " << nu << '\n';
    std::cout << std::setw(w) << "x"
	      << std::setw(w) << "ber"
	      << std::setw(w) << "bei"
	      << std::setw(w) << "ker"
	      << std::setw(w) << "kei"
	      << '\n';
    const auto del = Tp{1} / Tp{10};
    for (int i = 50; i <= 400; ++i)
      {
	auto x = del * i;
	auto kes = emsr::detail::kelvin_series(nu, x);
	auto kea = emsr::detail::kelvin_asymp(nu, x);
	std::cout << ' ' << std::setw(w) << kes.x
		  << ' ' << std::setw(w) << (kea.ber - kes.ber) / std::abs(kes.ber)
		  << ' ' << std::setw(w) << (kea.bei - kes.bei) / std::abs(kes.bei)
		  << ' ' << std::setw(w) << (kea.ker - kes.ker) / std::abs(kes.ker)
		  << ' ' << std::setw(w) << (kea.kei - kes.kei) / std::abs(kes.kei)
		  << '\n';
      }
    std::cout << '\n' << std::flush;
  }


/**
 * Run Kelvin with kelvin_series(n, x) which combines sums for integral order.
 */
template<typename Tp>
  void
  run_kelvin4(int n = 0, Tp proto = Tp{})
  {
    std::cout.precision(emsr::digits10(proto));
    std::cout << std::showpoint << std::scientific;
    auto w = 8 + std::cout.precision();

    std::cout << "\n\nPrint Kelvin functions computed by series expansions, n = " << n << '\n';
    std::cout << std::setw(w) << "x"
	      << std::setw(w) << "ber"
	      << std::setw(w) << "bei"
	      << std::setw(w) << "ker"
	      << std::setw(w) << "kei"
	      << '\n';
    const auto del = Tp{1} / Tp{10};
    for (int i = 0; i <= 200; ++i)
      {
	auto x = del * i;
	auto ke = emsr::detail::kelvin_series(n, x);
	std::cout << std::setw(w) << ke.x
		  << std::setw(w) << ke.ber
		  << std::setw(w) << ke.bei
		  << std::setw(w) << ke.ker
		  << std::setw(w) << ke.kei
		  << '\n';
      }
    std::cout << '\n' << std::flush;
  }


/**
 * Plot the scaled Kelvin functions.
 */
template<typename Tp>
  void
  plot_kelvin(std::string filename, Tp proto = Tp{})
  {
    const auto s_sqrt_2 = emsr::sqrt2_v<Tp>;
    const auto s_sqrt_pi = emsr::sqrtpi_v<Tp>;

    auto data = std::ofstream(filename);

    data.precision(emsr::digits10(proto));
    data << std::showpoint << std::scientific;
    auto w = 8 + data.precision();

    data << "\n\n";
    data << std::setw(w) << "x"
	 << std::setw(w) << "ber"
	 << std::setw(w) << "bei"
	 << std::setw(w) << "ker"
	 << std::setw(w) << "kei"
	 << std::setw(w) << "Wronskian"
	 << '\n';
    const auto del = Tp{1} / Tp{100};
    for (int i = 0; i <= +4000; ++i)
      {
	const auto x = del * i;
	const auto ber = emsr::detail::kelvin_ber(x);
	const auto bei = emsr::detail::kelvin_bei(x);
	const auto ker = emsr::detail::kelvin_ker(x);
	const auto kei = emsr::detail::kelvin_kei(x);
	const auto exf = std::exp(x / s_sqrt_2);
	const auto rt2x = s_sqrt_2 * std::sqrt(x);
	data << std::setw(w) << x
	     << std::setw(w) << s_sqrt_pi * rt2x * ber / exf
	     << std::setw(w) << s_sqrt_pi * rt2x * bei / exf
	     << std::setw(w) << rt2x * exf * ker / s_sqrt_pi
	     << std::setw(w) << rt2x * exf * kei / s_sqrt_pi
	     << '\n';
      }
    data << "\n\n";
  }


/**
 * Plot the scaled Kelvin functions for orders 0, 1/2, 1, 3/2, 2, 5.
 */
template<typename Tp>
  void
  plot_kelvin_order(std::string filename, Tp proto = Tp{})
  {
    const auto s_sqrt_2 = emsr::sqrt2_v<Tp>;
    const auto s_sqrt_pi = emsr::sqrtpi_v<Tp>;

    auto data = std::ofstream(filename);

    data.precision(emsr::digits10(proto));
    data << std::showpoint << std::scientific;
    auto w = 8 + data.precision();

    data << "\n\n";
    data << std::setw(w) << "x"
	 << std::setw(w) << "ber"
	 << std::setw(w) << "bei"
	 << std::setw(w) << "ker"
	 << std::setw(w) << "kei"
	 << std::setw(w) << "Wronskian"
	 << '\n';

    std::vector<Tp> nuv{Tp{0}, Tp{1}/Tp{2}, Tp{1},
			   Tp{3}/Tp{2}, Tp{2}, Tp{5}};
    const auto del = Tp{1} / Tp{100};
    for (auto nu : nuv)
      {
	for (int i = 0; i <= +4000; ++i)
	  {
	    const auto x = del * i;
	    const auto ber = emsr::detail::kelvin_ber(nu, x);
	    const auto bei = emsr::detail::kelvin_bei(nu, x);
	    const auto ker = emsr::detail::kelvin_ker(nu, x);
	    const auto kei = emsr::detail::kelvin_kei(nu, x);
	    const auto exf = std::exp(x / s_sqrt_2);
	    const auto rt2x = s_sqrt_2 * std::sqrt(x);
	    data << std::setw(w) << x
		 << std::setw(w) << s_sqrt_pi * rt2x * ber / exf
		 << std::setw(w) << s_sqrt_pi * rt2x * bei / exf
		 << std::setw(w) << rt2x * exf * ker / s_sqrt_pi
		 << std::setw(w) << rt2x * exf * kei / s_sqrt_pi
		 << '\n';
	  }
	data << "\n\n";
      }
  }


int
main(int n_app_args, char** arg)
{
  std::string plot_data_dir = ".";
  if (n_app_args > 1)
    plot_data_dir = arg[1];

  run_kelvin4<double>();
  run_kelvin3<double>();
  run_kelvin2<double>();
  run_kelvin1<double>();

  diff_kelvin3<double>();
  diff_kelvin2<double>();

  // FIXME
  //plot_kelvin<float>(plot_data_dir + '/' + "kelvin_float.txt");
  //plot_kelvin<double>(plot_data_dir + '/' + "kelvin_double.txt");
  //plot_kelvin<long double>(plot_data_dir + '/' + "kelvin_long_double.txt");

  //plot_kelvin_order<float>(plot_data_dir + '/' + "kelvin_order_float.txt");
  //plot_kelvin_order<double>(plot_data_dir + '/' + "kelvin_order_double.txt");
  //plot_kelvin_order<long double>(plot_data_dir + '/' + "kelvin_order_long_double.txt");
}
