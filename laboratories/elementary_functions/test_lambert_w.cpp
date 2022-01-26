/**
 *
 */

#include <cmath>
#include <iostream>
#include <iomanip>

  /**
   * This is the third-order Halley's root finding algorithm for Lambert W.
   * The radius of convergence is 1/e but it staggers on pretty well up to above 3.
   */
  template<typename _Tp>
    _Tp
    lambert_w_series(_Tp z)
    {
      const auto _S_eps = emsr::epsilon(z);
      const auto _S_max_iter = 1000u;

      auto _W = z * (_Tp{1} - z);
      auto term = -z * z;
      for (auto k = 3u; k < _S_max_iter; ++k)
	{
	  term *= -z * std::pow(_Tp(k) / _Tp(k - 1), k - 2);
	  _W += term;
	  if (std::abs(term) < _S_eps * std::abs(_W))
	    break;
	}
      return _W;
    }

  /**
   * This is the second-order Newton root finding algorithm for Lambert W.
   */
  template<typename _Tp>
    _Tp
    lambert_w_newton(_Tp z, _Tp _W = _Tp{1})
    {
      const auto _S_eps = emsr::epsilon(z);
      const auto _S_max_iter = 1000u;

      auto wk = _W;
      for (auto k = 0u; k < _S_max_iter; ++k)
	{
          const auto expwk = std::exp(wk);
          const auto wexpwk = wk * expwk;
	  const auto wkp1 = wk - (wexpwk - z)
				   / (_Tp{1} + wk) / expwk;
	  const auto del = std::abs(wkp1 - wk);
	  wk = wkp1;
	  if (del < _S_eps)
	    break;
	}
      return wk;
    }

  /**
   * This is the third-order Halley root finding algorithm for Lambert W.
   */
  template<typename _Tp>
    _Tp
    lambert_w_halley(_Tp z, _Tp _W = _Tp{1})
    {
      const auto _S_eps = emsr::epsilon(z);
      const auto _S_max_iter = 1000u;

      auto wk = _W;
      for (auto k = 0u; k < _S_max_iter; ++k)
	{
          const auto expwk = std::exp(wk);
	  const auto fact = wk * expwk - z;
          const auto wkp1 = wk - fact
		      / ((wk + 1) * expwk - (wk + 2) * fact / (2 * wk + 2));
	  const auto del = std::abs(wkp1 - wk);
	  wk = wkp1;
	  if (del < _S_eps)
	    break;
	}
      return wk;
    }


  /**
   * This is the fifth-order Schroder's update for Lambert W.
   */
  template<typename _Tp>
    _Tp
    lambert_w_schroder(_Tp z, _Tp _W)
    {
      const auto y = z * std::exp(-_W);
      const auto f0 = _W - y;
      const auto f1 = _Tp{1} + y;
      const auto f2 = y;
      const auto f11 = f1 * f1;
      const auto f0y = f0 * y;
      const auto f00y = f0 * f0y;
      return _W - 4 *f0 * (6 * f1 * (f11 + f0y) + f00y)
		/ (f11 * (24 * f11 + 36 * f0y)
		 + 6 * f00y * (14 * y + 8 + f0));
    }


/**
 * This is the asymptotic log series for @f$ W_0(z) @f$
 * as @f$ z \rightarrow \infty @f$.
 * @f[
 *    W_0(z) = \xi - ln(\xi) + \frac{ln(\xi)}{\xi}
 *           - \frac{ln(\xi)}{\xi^2} + \frac{(ln(\xi))^2}{\xi^2}
 * @f]
 * where @f$ \xi = ln(z) @f$
 */
template<typename _Tp>
  _Tp
  lambert_w_0_log_series(_Tp z)
  {
    const auto xi = std::log(z);
    const auto lnxi = std::log(xi);
    return xi - lnxi * (_Tp{1} - (_Tp{1} / xi)
		* (_Tp{1} - (_Tp{1} / xi) * (_Tp{1} - lnxi / _Tp{2})));
  }


/**
 * This is the asymptotic log series for @f$ W_1(z) @f$
 * as @f$ z \rightarrow 0- @f$.
 * @f[
 *    W_1(z) = -\eta - ln(\eta) - \frac{ln(\eta)}{\eta}
 *           - \frac{ln(\eta)}{\eta^2} - \frac{(ln(\eta))^2}{\eta^2}
 * @f]
 * where @f$ \eta = ln(-1/z) @f$
 */
template<typename _Tp>
  _Tp
  lambert_w_1_log_series(_Tp z)
  {
    const auto eta = std::log(-_Tp{1} / z);
    const auto lneta = std::log(eta);
    return -eta - lneta * (_Tp{1} + (_Tp{1} / eta)
		* (_Tp{1} + (_Tp{1} / eta) * (_Tp{1} + lneta / _Tp{2})));
  }


template<typename _Tp>
  void
  test_lambert_w(_Tp proto = _Tp{})
  {
    std::cout.precision(emsr::digits10(proto));
    auto w = std::cout.precision() + 8;
    std::cout << std::showpoint << std::scientific;

    const int N0 = 100;
    const auto _S_e = emsr::e_v<_Tp>;
    const auto _S_1_e = _Tp{1} / _S_e;
    const auto del0 = (_S_e + _S_1_e) / N0;

    const int Nm1 = 50;
    const auto delm1 = _S_1_e / Nm1;

    std::cout << '\n';
    for (int i = 0; i <= N0; ++i)
      {
	auto z = -_S_1_e  + del0 * i;
        auto W_newton = lambert_w_newton(z);
        auto W_halley = lambert_w_halley(z);
        auto W_series = lambert_w_series(z);
	std::cout << ' ' << std::setw(w) << z
		  << ' ' << std::setw(w) << W_newton
		  << ' ' << std::setw(w) << W_halley
		  << ' ' << std::setw(w) << W_series
		  << '\n';
      }

    std::cout << '\n';
    std::cout << '\n';
    for (int i = 0; i <= Nm1; ++i)
      {
	auto z = -_S_1_e  + delm1 * i;
        auto W_newton = lambert_w_newton(z, _Tp{-2});
        auto W_halley = lambert_w_halley(z, _Tp{-2});
	std::cout << ' ' << std::setw(w) << z
		  << ' ' << std::setw(w) << W_newton
		  << ' ' << std::setw(w) << W_halley
		  << '\n';
      }

    std::cout << '\n';
    std::cout << '\n';
    auto term = _Tp{-1};
    for (auto k = 3u; k < 50; ++k)
      {
	term *= -std::pow(_Tp(k) / _Tp(k - 1), k - 2);
	std::cout << ' ' << term << '\n';
      }
  }

int
main()
{
  test_lambert_w(1.0);
}
