/**
 *
 */

#include <cmath>
#include <limits>
#include <iostream>
#include <iomanip>
#include <complex>
#include <emsr/math_constants.h>

  /**
   * Compute the Kelvin function @f$ ber(x) @f$ by series summation.
   *
   * @f[
   *    ber(x) = \sum_{k=0}^\infty \frac{[-(x/2)^4]^k}{[(2k)!]^2}
   * @f]
   */
  template<typename _Tp>
    std::pair<_Tp, _Tp>
    __kelvin_ber_series(_Tp __x)
    {
      constexpr auto _S_eps = std::numeric_limits<_Tp>::epsilon();
      constexpr int _S_max_iter = 100;

      if (__x == _Tp{0})
	return std::make_pair(_Tp{1}, _Tp{0});

      const auto __y = __x / _Tp{2};
      const auto __y2 = __y * __y;
      const auto __y4 = __y2 * __y2;

      auto __term = _Tp{1};
      auto __ber = __term;
      auto __berd = _Tp{0};
      for (int __k = 1; __k < _S_max_iter; ++__k)
	{
	  const auto __fact = _Tp{1} / _Tp(2 * __k - 1) / _Tp(2 * __k);
	  __term *= -__y4 * __fact * __fact;
	  __ber += __term;
	  __berd += _Tp(2 * __k) * __term / __y;
	  if (std::abs(__term) < _S_eps * std::abs(__ber))
	    break;
	}

      return std::make_pair(__ber, __berd);
    }

  /**
   * Compute the Kelvin function @f$ bei(x) @f$ by series summation.
   *
   * @f[
   *    bei(x) = (x/2)^2\sum_{k=0}^\infty \frac{[-(x/2)^4]^k}{[(2k+1)!]^2}
   * @f]
   */
  template<typename _Tp>
    std::pair<_Tp, _Tp>
    __kelvin_bei_series(_Tp __x)
    {
      constexpr auto _S_eps = std::numeric_limits<_Tp>::epsilon();
      constexpr int _S_max_iter = 100;

      if (__x == _Tp{0})
	return std::make_pair(_Tp{0}, _Tp{0});

      const auto __y = __x / _Tp{2};
      const auto __y2 = __y * __y;
      const auto __y4 = __y2 * __y2;

      auto __term = __y2;
      auto __bei = __term;
      auto __beid = __y;
      for (int __k = 1; __k < _S_max_iter; ++__k)
	{
	  const auto __fact = _Tp{1} / _Tp(2 * __k + 1) / _Tp(2 * __k);
	  __term *= -__y4 * __fact * __fact;
	  __bei += __term;
	  __beid += _Tp(2 * __k + 1) * __term / __y;
	  if (std::abs(__term) < _S_eps * std::abs(__bei))
	    break;
	}

      return std::make_pair(__bei, __beid);
    }

  template<typename _Tp>
    std::pair<_Tp, _Tp>
    __kelvin_ker_series(_Tp __x)
    {
      constexpr auto _S_eps = std::numeric_limits<_Tp>::epsilon();
      constexpr int _S_max_iter = 100;
      constexpr auto _S_gamma_e = emsr::egamma_v<_Tp>;
      constexpr auto _S_pi_4 = emsr::pi_v<_Tp> / _Tp{4};
      constexpr auto _S_inf = std::numeric_limits<_Tp>::infinity();

      if (__x == _Tp{0})
	return std::make_pair(_S_inf, -_S_inf);

      const auto __y = __x / _Tp{2};
      const auto __y2 = __y * __y;
      const auto __y4 = __y2 * __y2;

      // phi(1) = -gamma_E <=> H_0 = 0.
      auto _H = _Tp{0};
      auto __term = _Tp{1};
      auto __ker = _Tp{0}; 
      auto __kerd = _Tp{0};
      for (int __k = 1; __k < _S_max_iter; ++__k)
	{
	  _H += _Tp{1} / _Tp(2 * __k - 1) + _Tp{1} / _Tp(2 * __k);
	  const auto __fact = _Tp{1} / _Tp(2 * __k - 1) / _Tp(2 * __k);
	  __term *= -__y4 * __fact * __fact;
	  __ker += __term * _H;
	  __kerd += _Tp(2 * __k) * __term * _H / __y;
	  if (std::abs(__term * _H) < _S_eps * std::abs(__ker))
	    break;
	}

      const auto __ber = __kelvin_ber_series(__x);
      const auto __bei = __kelvin_bei_series(__x);
      __ker = __ker
	    - (std::log(__y) + _S_gamma_e) * __ber.first
	    + _S_pi_4 * __bei.first;
      __kerd = __kerd
	     - __ber.first / __x - (std::log(__y) + _S_gamma_e) * __ber.second
	     + _S_pi_4 * __bei.second;

      return std::make_pair(__ker, __kerd);
    }

  template<typename _Tp>
    std::pair<_Tp, _Tp>
    __kelvin_kei_series(_Tp __x)
    {
      constexpr auto _S_eps = std::numeric_limits<_Tp>::epsilon();
      constexpr int _S_max_iter = 100;
      constexpr auto _S_gamma_e = emsr::egamma_v<_Tp>;
      constexpr auto _S_pi_4 = emsr::pi_v<_Tp> / _Tp{4};

      if (__x == _Tp{0})
	return std::make_pair(-_S_pi_4, _Tp{0});

      const auto __y = __x / _Tp{2};
      const auto __y2 = __y * __y;
      const auto __y4 = __y2 * __y2;

      auto _H = _Tp{1};
      auto __term = __y2;
      auto __kei = __term;
      auto __keid = __y;
      for (int __k = 1; __k < _S_max_iter; ++__k)
	{
	  _H += _Tp{1} / _Tp(2 * __k) + _Tp{1} / _Tp(2 * __k + 1);
	  const auto __fact = _Tp{1} / _Tp(2 * __k) / _Tp(2 * __k + 1);
	  __term *= -__y4 * __fact * __fact;
	  __kei += __term * _H;
	  __keid += _Tp(2 * __k + 1) * __term * _H / __y;
	  if (std::abs(__term * _H) < _S_eps * std::abs(__kei))
	    break;
	}

      const auto __ber = __kelvin_ber_series(__x);
      const auto __bei = __kelvin_bei_series(__x);
      __kei = __kei
	    - (std::log(__y) + _S_gamma_e) * __bei.first
	    - _S_pi_4 * __ber.first;
      __keid = __keid
	     - __bei.first / __x - (std::log(__y) + _S_gamma_e) * __bei.second
	     - _S_pi_4 * __ber.second;

      return std::make_pair(__kei, __keid);
    }

  template<typename _TpNu, typename _Tp>
    struct Kelvin
    {
      _TpNu __nu;
      _Tp __x;
      _Tp __ber_value, __ber_deriv;
      _Tp __bei_value, __bei_deriv;
      _Tp __ker_value, __ker_deriv;
      _Tp __kei_value, __kei_deriv;

      _Tp
      __ber_deriv2() const
      { return deriv2(__ber_value, __ber_deriv); }
      _Tp
      __bei_deriv2() const
      { return deriv2(__bei_value, __bei_deriv); }
      _Tp
      __ker_deriv2() const
      { return deriv2(__ker_value, __ker_deriv); }
      _Tp
      __kei_deriv2() const
      { return deriv2(__kei_value, __kei_deriv); }

      _Tp
      deriv2(_Tp __val, _Tp __der)
      {
	if (__x == _Tp{0})
	  return std::numeric_limits<_Tp>::infinity();
	else
	  return -(__x * __der + (__x - __nu) * (__x + __nu) * __val)
	       / __x / __x;
      }
    };

  /**
   * Compute the Kelvin functions by series summation.
   *
   * @f[
   *    ber(x) = \sum_{k=0}^\infty \frac{-(x/2)^{4k}}{[(2k)!]^2}
   * @f]
   * @f[
   *    bei(x) = \frac{x^2}{4} \sum_{k=0}^\infty \frac{-(x/2)^{4k}}{[(2k+1)!]^2}
   * @f]
   * @f[
   *    ker(x) = -\left[log\left(\frac{x}{2}\right) + \gamma_E\right] ber(x)
   *           + \frac{\pi}{4} bei(x)
   *           + \sum_{k=0}^\infty \frac{-(x/2)^{4k}}{[(2k)!]^2}H_{2k}
   * @f]
   * @f[
   *    kei(x) = -\left[log\left(\frac{x}{2}\right) + \gamma_E\right] bei(x)
   *           - \frac{\pi}{4} ber(x)
   *           + \frac{x^2}{4}
   *               \sum_{k=0}^\infty \frac{-(x/2)^{4k}}{[(2k+1)!]^2}H_{2k+1}
   * @f]
   * where
   * @f[
   *    H_n = \sum_{k=1}^n \frac{1}{k}
   * @f]
   * is the harmonic number.
   */
  template<typename _Tp>
    Kelvin<int, _Tp>
    __kelvin_series(_Tp __x)
    {
      //using _BasicSum = emsr::BasicSum<_Tp>;
      //using _WijnSum = emsr::VanWijngaardenSum<_Tp>;
      //using _WenigerDeltaWijnSum = emsr::WenigerDeltaSum<_WijnSum>;

      constexpr auto _S_eps = std::numeric_limits<_Tp>::epsilon();
      constexpr int _S_max_iter = 100;
      constexpr auto _S_gamma_e = emsr::egamma_v<_Tp>;
      constexpr auto _S_pi_4 = emsr::pi_v<_Tp> / _Tp{4};
      constexpr auto _S_inf = std::numeric_limits<_Tp>::infinity();

      if (__x == _Tp{0})
	return {0, __x,
		_Tp{1}, _Tp{0}, _Tp{0}, _Tp{0},
		_S_inf, -_S_inf, -_S_pi_4, _Tp{0}};

      const auto __y = __x / _Tp{2};
      const auto __y2 = __y * __y;

      //_BasicSum _H_n;
      //_WenigerDeltaWijnSum __ber, __bei, __ker, __kei;
      //_WenigerDeltaWijnSum __berp, __beip, __kerp, __keip;
      auto __term = _Tp{1};
      auto __ber = __term;
      auto __berp = _Tp{0};
      auto _H = _Tp{0}; // phi(1) = -gamma_E <=> H_0 = 0.
      auto __ker = __term * _H; 
      auto __kerp = _Tp{0};

      __term = __y2;
      auto __bei = __term;
      auto __beip = __y;
      _H = _Tp{1};
      auto __kei = __term * _H;
      auto __keip = __y;

      for (int __k = 1; __k < _S_max_iter; ++__k)
	{
	  bool __done = true;

	  const auto __tk = _Tp(2 * __k);
	  const auto __factr = _Tp{1} / __tk;
	  __term *= -__y2 * __factr * __factr;
	  __ber += __term;
	  __berp += __tk * __term / __y;
	  __done = __done && std::abs(__term) < _S_eps * std::abs(__ber);
	  _H += __factr;
	  __ker += __term * _H;
	  __kerp += __tk * __term * _H / __y;
	  __done = __done && std::abs(__term * _H) < _S_eps * std::abs(__ker);

	  const auto __tkp1 = __tk + _Tp{1};
	  const auto __facti = _Tp{1} / __tkp1;
	  __term *= __y2 * __facti * __facti;
	  __bei += __term;
	  __beip += __tkp1 * __term / __y;
	  __done = __done && std::abs(__term) < _S_eps * std::abs(__bei);
	  _H += __facti;
	  __kei += __term * _H;
	  __keip += __tkp1 * __term * _H / __y;
	  __done = __done && std::abs(__term * _H) < _S_eps * std::abs(__kei);

	  if (__done)
	    break;
	}

      __ker = __ker
	    - (std::log(__y) + _S_gamma_e) * __ber
	    + _S_pi_4 * __bei;
      __kerp = __kerp
	     - __ber / __x - (std::log(__y) + _S_gamma_e) * __berp
	     + _S_pi_4 * __beip;

      __kei = __kei
	    - (std::log(__y) + _S_gamma_e) * __bei
	    - _S_pi_4 * __ber;
      __keip = __keip
	     - __bei / __x - (std::log(__y) + _S_gamma_e) * __beip
	     - _S_pi_4 * __berp;

      return {0, __x,
	      __ber, __berp, __bei, __beip,
	      __ker, __kerp, __kei, __keip};
    }

  /**
   * Compute the Kelvin functions by asymptotic series expansion.
   *
   *
   *
   *
   */
  template<typename _Tp>
    Kelvin<int, _Tp>
    __kelvin_asymp(_Tp __x)
    {
    }

  /**
   * Compute the Kelvin functions of integer order n by series expansion.
   */
  template<typename _Tp>
    Kelvin<int, _Tp>
    __kelvin_series(int __n, _Tp __x)
    {
return {__n, __x,
	_Tp{0}, _Tp{0}, _Tp{0}, _Tp{0},
	_Tp{0}, _Tp{0}, _Tp{0}, _Tp{0}};
    }

template<typename _Tp>
  void
  help_kelvin()
  {
    std::cout.precision(std::numeric_limits<_Tp>::digits10);
    std::cout << std::showpoint << std::scientific;
    auto w = 8 + std::cout.precision();

    std::cout << "\n\nPrint Kelvin functions computed by series expansions\n";
    std::cout << std::setw(w) << "x"
	      << std::setw(w) << "ber(x)"
	      << std::setw(w) << "bei(x)"
	      << std::setw(w) << "ker(x)"
	      << std::setw(w) << "kei(x)"
	      << std::setw(w) << "ber'(x)"
	      << std::setw(w) << "bei'(x)"
	      << std::setw(w) << "ker'(x)"
	      << std::setw(w) << "kei'(x)"
	      << '\n';
    const auto del = _Tp{1}/_Tp{10};
    for (int i = 0; i <= 200; ++i)
      {
	const auto x = del * i;
	const auto ber = __kelvin_ber_series(x);
	const auto bei = __kelvin_bei_series(x);
	const auto ker = __kelvin_ker_series(x);
	const auto kei = __kelvin_kei_series(x);
	std::cout << std::setw(w) << x
		  << std::setw(w) << ber.first
		  << std::setw(w) << bei.first
		  << std::setw(w) << ker.first
		  << std::setw(w) << kei.first
		  << std::setw(w) << ber.second
		  << std::setw(w) << bei.second
		  << std::setw(w) << ker.second
		  << std::setw(w) << kei.second
		  << '\n';
      }
    std::cout << '\n' << std::flush;
  }

template<typename _Tp>
  void
  help_kelvin2()
  {
    std::cout.precision(std::numeric_limits<_Tp>::digits10);
    std::cout << std::showpoint << std::scientific;
    auto w = 8 + std::cout.precision();

    std::cout << "\n\nPrint Kelvin functions computed by series expansions\n";
    std::cout << std::setw(w) << "x"
	      << std::setw(w) << "ber(x)"
	      << std::setw(w) << "bei(x)"
	      << std::setw(w) << "ker(x)"
	      << std::setw(w) << "kei(x)"
	      << std::setw(w) << "ber'(x)"
	      << std::setw(w) << "bei'(x)"
	      << std::setw(w) << "ker'(x)"
	      << std::setw(w) << "kei'(x)"
	      << '\n';
    const auto del = _Tp{1}/_Tp{10};
    for (int i = 0; i <= 200; ++i)
      {
	const auto x = del * i;
	const auto klv = __kelvin_series(x);
	std::cout << std::setw(w) << x
		  << std::setw(w) << klv.__ber_value
		  << std::setw(w) << klv.__bei_value
		  << std::setw(w) << klv.__ker_value
		  << std::setw(w) << klv.__kei_value
		  << std::setw(w) << klv.__ber_deriv
		  << std::setw(w) << klv.__bei_deriv
		  << std::setw(w) << klv.__ker_deriv
		  << std::setw(w) << klv.__kei_deriv
		  << '\n';
      }
    std::cout << '\n' << std::flush;
  }

template<typename _Tp>
  void
  help_kelvin3()
  {
    std::cout.precision(std::numeric_limits<_Tp>::digits10);
    std::cout << std::showpoint << std::scientific;
    auto w = 8 + std::cout.precision();

    std::cout << "\n\nPrint Kelvin functions computed by series expansions\n";
    std::cout << std::setw(w) << "x"
	      << std::setw(w) << "ber(x)"
	      << std::setw(w) << "bei(x)"
	      //<< std::setw(w) << "ber'(x)"
	      //<< std::setw(w) << "bei'(x)"
	      << '\n';
    const auto del = _Tp{1}/_Tp{10};
    const auto ph = std::polar(_Tp{1}, 3 * emsr::pi_v<_Tp> / 4);
    const auto pi = emsr::pi_v<_Tp>;
    for (int i = 0; i <= 200; ++i)
      {
	const auto x = del * i;
	const auto klv = std::__detail::__cyl_bessel_ij_series(0, ph * x, -1, 200);
	//const auto klv = std::__detail::__cyl_bessel_ij_series(0, x, -1, 200);

	std::cout << std::setw(w) << x
		  << std::setw(w) << std::real(klv)
		  << std::setw(w) << std::imag(klv)
		  //<< std::setw(w) << std::real(klv.__J_value)
		  //<< std::setw(w) << std::imag(klv.__J_value)
		  //<< std::setw(w) << std::real(klv.__J_deriv)
		  //<< std::setw(w) << std::imag(klv.__J_deriv)
		  << '\n';
      }
    std::cout << std::endl;
    for (int i = 100; i <= 400; ++i)
      {
	const auto x = del * i;
	const auto klv = std::__detail::__cyl_bessel_jn_asymp(0, ph * x);
	//const auto klv = std::__detail::__cyl_bessel_jn_asymp(0, x);

	std::cout << std::setw(w) << x
		  << std::setw(w) << std::real(klv.__J_value)
		  << std::setw(w) << std::imag(klv.__J_value)
		  << std::setw(w) << std::real(klv.__J_deriv)
		  << std::setw(w) << std::imag(klv.__J_deriv)
		  << std::setw(2 * w + 4) << pi * ph * x * klv.__Wronskian() / _Tp{2}
		  << '\n';
      }
    std::cout << std::endl;
  }

int
main()
{
  help_kelvin<double>();
  help_kelvin2<double>();
  help_kelvin3<double>();
}
