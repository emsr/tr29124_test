/**
 *
 */

#include <cmath>
#include <limits>
#include <iostream>
#include <iomanip>
#include <complex>

#include <emsr/math_constants.h>
#include <emsr/sf_bessel.h>

  /**
   * Compute the Kelvin function @f$ ber(x) @f$ by series summation.
   *
   * @f[
   *    ber(x) = \sum_{k=0}^\infty \frac{[-(x/2)^4]^k}{[(2k)!]^2}
   * @f]
   */
  template<typename _Tp>
    std::pair<_Tp, _Tp>
    kelvin_ber_series(_Tp x)
    {
      constexpr auto s_eps = std::numeric_limits<_Tp>::epsilon();
      constexpr int s_max_iter = 100;

      if (x == _Tp{0})
	return std::make_pair(_Tp{1}, _Tp{0});

      const auto y = x / _Tp{2};
      const auto y2 = y * y;
      const auto y4 = y2 * y2;

      auto term = _Tp{1};
      auto ber = term;
      auto berd = _Tp{0};
      for (int k = 1; k < s_max_iter; ++k)
	{
	  const auto fact = _Tp{1} / _Tp(2 * k - 1) / _Tp(2 * k);
	  term *= -y4 * fact * fact;
	  ber += term;
	  berd += _Tp(2 * k) * term / y;
	  if (std::abs(term) < s_eps * std::abs(ber))
	    break;
	}

      return std::make_pair(ber, berd);
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
    kelvin_bei_series(_Tp x)
    {
      constexpr auto s_eps = std::numeric_limits<_Tp>::epsilon();
      constexpr int s_max_iter = 100;

      if (x == _Tp{0})
	return std::make_pair(_Tp{0}, _Tp{0});

      const auto y = x / _Tp{2};
      const auto y2 = y * y;
      const auto y4 = y2 * y2;

      auto term = y2;
      auto bei = term;
      auto beid = y;
      for (int k = 1; k < s_max_iter; ++k)
	{
	  const auto fact = _Tp{1} / _Tp(2 * k + 1) / _Tp(2 * k);
	  term *= -y4 * fact * fact;
	  bei += term;
	  beid += _Tp(2 * k + 1) * term / y;
	  if (std::abs(term) < s_eps * std::abs(bei))
	    break;
	}

      return std::make_pair(bei, beid);
    }

  template<typename _Tp>
    std::pair<_Tp, _Tp>
    kelvin_ker_series(_Tp x)
    {
      constexpr auto s_eps = std::numeric_limits<_Tp>::epsilon();
      constexpr int s_max_iter = 100;
      constexpr auto s_gamma_e = emsr::egamma_v<_Tp>;
      constexpr auto s_pi_4 = emsr::pi_v<_Tp> / _Tp{4};
      constexpr auto s_inf = std::numeric_limits<_Tp>::infinity();

      if (x == _Tp{0})
	return std::make_pair(s_inf, -s_inf);

      const auto y = x / _Tp{2};
      const auto y2 = y * y;
      const auto y4 = y2 * y2;

      // phi(1) = -gamma_E <=> H_0 = 0.
      auto _H = _Tp{0};
      auto term = _Tp{1};
      auto ker = _Tp{0}; 
      auto kerd = _Tp{0};
      for (int k = 1; k < s_max_iter; ++k)
	{
	  _H += _Tp{1} / _Tp(2 * k - 1) + _Tp{1} / _Tp(2 * k);
	  const auto fact = _Tp{1} / _Tp(2 * k - 1) / _Tp(2 * k);
	  term *= -y4 * fact * fact;
	  ker += term * _H;
	  kerd += _Tp(2 * k) * term * _H / y;
	  if (std::abs(term * _H) < s_eps * std::abs(ker))
	    break;
	}

      const auto ber = kelvin_ber_series(x);
      const auto bei = kelvin_bei_series(x);
      ker = ker
	    - (std::log(y) + s_gamma_e) * ber.first
	    + s_pi_4 * bei.first;
      kerd = kerd
	     - ber.first / x - (std::log(y) + s_gamma_e) * ber.second
	     + s_pi_4 * bei.second;

      return std::make_pair(ker, kerd);
    }

  template<typename _Tp>
    std::pair<_Tp, _Tp>
    kelvin_kei_series(_Tp x)
    {
      constexpr auto s_eps = std::numeric_limits<_Tp>::epsilon();
      constexpr int s_max_iter = 100;
      constexpr auto s_gamma_e = emsr::egamma_v<_Tp>;
      constexpr auto s_pi_4 = emsr::pi_v<_Tp> / _Tp{4};

      if (x == _Tp{0})
	return std::make_pair(-s_pi_4, _Tp{0});

      const auto y = x / _Tp{2};
      const auto y2 = y * y;
      const auto y4 = y2 * y2;

      auto _H = _Tp{1};
      auto term = y2;
      auto kei = term;
      auto keid = y;
      for (int k = 1; k < s_max_iter; ++k)
	{
	  _H += _Tp{1} / _Tp(2 * k) + _Tp{1} / _Tp(2 * k + 1);
	  const auto fact = _Tp{1} / _Tp(2 * k) / _Tp(2 * k + 1);
	  term *= -y4 * fact * fact;
	  kei += term * _H;
	  keid += _Tp(2 * k + 1) * term * _H / y;
	  if (std::abs(term * _H) < s_eps * std::abs(kei))
	    break;
	}

      const auto ber = kelvin_ber_series(x);
      const auto bei = kelvin_bei_series(x);
      kei = kei
	    - (std::log(y) + s_gamma_e) * bei.first
	    - s_pi_4 * ber.first;
      keid = keid
	     - bei.first / x - (std::log(y) + s_gamma_e) * bei.second
	     - s_pi_4 * ber.second;

      return std::make_pair(kei, keid);
    }

  template<typename _TpNu, typename _Tp>
    struct Kelvin
    {
      _TpNu nu;
      _Tp x;
      _Tp ber_value, ber_deriv;
      _Tp bei_value, bei_deriv;
      _Tp ker_value, ker_deriv;
      _Tp kei_value, kei_deriv;

      _Tp
      ber_deriv2() const
      { return deriv2(ber_value, ber_deriv); }
      _Tp
      bei_deriv2() const
      { return deriv2(bei_value, bei_deriv); }
      _Tp
      ker_deriv2() const
      { return deriv2(ker_value, ker_deriv); }
      _Tp
      kei_deriv2() const
      { return deriv2(kei_value, kei_deriv); }

      _Tp
      deriv2(_Tp val, _Tp der)
      {
	if (x == _Tp{0})
	  return std::numeric_limits<_Tp>::infinity();
	else
	  return -(x * der + (x - nu) * (x + nu) * val)
	       / x / x;
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
    kelvin_series(_Tp x)
    {
      //using _BasicSum = emsr::BasicSum<_Tp>;
      //using _WijnSum = emsr::VanWijngaardenSum<_Tp>;
      //using _WenigerDeltaWijnSum = emsr::WenigerDeltaSum<_WijnSum>;

      constexpr auto s_eps = std::numeric_limits<_Tp>::epsilon();
      constexpr int s_max_iter = 100;
      constexpr auto s_gamma_e = emsr::egamma_v<_Tp>;
      constexpr auto s_pi_4 = emsr::pi_v<_Tp> / _Tp{4};
      constexpr auto s_inf = std::numeric_limits<_Tp>::infinity();

      if (x == _Tp{0})
	return {0, x,
		_Tp{1}, _Tp{0}, _Tp{0}, _Tp{0},
		s_inf, -s_inf, -s_pi_4, _Tp{0}};

      const auto y = x / _Tp{2};
      const auto y2 = y * y;

      //_BasicSum _H_n;
      //_WenigerDeltaWijnSum ber, bei, ker, kei;
      //_WenigerDeltaWijnSum berp, beip, kerp, keip;
      auto term = _Tp{1};
      auto ber = term;
      auto berp = _Tp{0};
      auto _H = _Tp{0}; // phi(1) = -gamma_E <=> H_0 = 0.
      auto ker = term * _H; 
      auto kerp = _Tp{0};

      term = y2;
      auto bei = term;
      auto beip = y;
      _H = _Tp{1};
      auto kei = term * _H;
      auto keip = y;

      for (int k = 1; k < s_max_iter; ++k)
	{
	  bool done = true;

	  const auto tk = _Tp(2 * k);
	  const auto factr = _Tp{1} / tk;
	  term *= -y2 * factr * factr;
	  ber += term;
	  berp += tk * term / y;
	  done = done && std::abs(term) < s_eps * std::abs(ber);
	  _H += factr;
	  ker += term * _H;
	  kerp += tk * term * _H / y;
	  done = done && std::abs(term * _H) < s_eps * std::abs(ker);

	  const auto tkp1 = tk + _Tp{1};
	  const auto facti = _Tp{1} / tkp1;
	  term *= y2 * facti * facti;
	  bei += term;
	  beip += tkp1 * term / y;
	  done = done && std::abs(term) < s_eps * std::abs(bei);
	  _H += facti;
	  kei += term * _H;
	  keip += tkp1 * term * _H / y;
	  done = done && std::abs(term * _H) < s_eps * std::abs(kei);

	  if (done)
	    break;
	}

      ker = ker
	    - (std::log(y) + s_gamma_e) * ber
	    + s_pi_4 * bei;
      kerp = kerp
	     - ber / x - (std::log(y) + s_gamma_e) * berp
	     + s_pi_4 * beip;

      kei = kei
	    - (std::log(y) + s_gamma_e) * bei
	    - s_pi_4 * ber;
      keip = keip
	     - bei / x - (std::log(y) + s_gamma_e) * beip
	     - s_pi_4 * berp;

      return {0, x,
	      ber, berp, bei, beip,
	      ker, kerp, kei, keip};
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
    kelvin_asymp(_Tp x)
    {
    }

  /**
   * Compute the Kelvin functions of integer order n by series expansion.
   */
  template<typename _Tp>
    Kelvin<int, _Tp>
    kelvin_series(int n, _Tp x)
    {
return {n, x,
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
	const auto ber = kelvin_ber_series(x);
	const auto bei = kelvin_bei_series(x);
	const auto ker = kelvin_ker_series(x);
	const auto kei = kelvin_kei_series(x);
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
	const auto klv = kelvin_series(x);
	std::cout << std::setw(w) << x
		  << std::setw(w) << klv.ber_value
		  << std::setw(w) << klv.bei_value
		  << std::setw(w) << klv.ker_value
		  << std::setw(w) << klv.kei_value
		  << std::setw(w) << klv.ber_deriv
		  << std::setw(w) << klv.bei_deriv
		  << std::setw(w) << klv.ker_deriv
		  << std::setw(w) << klv.kei_deriv
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
	const auto klv = emsr::detail::cyl_bessel_ij_series(0, ph * x, -1, 200);
	//const auto klv = emsr::detail::cyl_bessel_ij_series(0, x, -1, 200);

	std::cout << std::setw(w) << x
		  << std::setw(w) << std::real(klv)
		  << std::setw(w) << std::imag(klv)
		  //<< std::setw(w) << std::real(klv.J_value)
		  //<< std::setw(w) << std::imag(klv.J_value)
		  //<< std::setw(w) << std::real(klv.J_deriv)
		  //<< std::setw(w) << std::imag(klv.J_deriv)
		  << '\n';
      }
    std::cout << std::endl;
    for (int i = 100; i <= 400; ++i)
      {
	const auto x = del * i;
	const auto klv = emsr::detail::cyl_bessel_jn_asymp(0, ph * x);
	//const auto klv = emsr::detail::cyl_bessel_jn_asymp(0, x);

	std::cout << std::setw(w) << x
		  << std::setw(w) << std::real(klv.J_value)
		  << std::setw(w) << std::imag(klv.J_value)
		  << std::setw(w) << std::real(klv.J_deriv)
		  << std::setw(w) << std::imag(klv.J_deriv)
		  << std::setw(2 * w + 4) << pi * ph * x * klv.Wronskian() / _Tp{2}
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
