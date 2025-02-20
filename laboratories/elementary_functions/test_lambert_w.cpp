/**
 *
 */

#include <cmath>
#include <iostream>
#include <iomanip>

#include <emsr/numeric_limits.h>
#include <emsr/math_constants.h>

  /**
   * This is the third-order Halley's root finding algorithm for Lambert W.
   * The radius of convergence is 1/e but it staggers on pretty well up to above 3.
   */
  template<typename Tp>
    Tp
    lambert_w_series(Tp z)
    {
      const auto s_eps = emsr::epsilon(z);
      const auto s_max_iter = 1000u;

      auto _W = z * (Tp{1} - z);
      auto term = -z * z;
      for (auto k = 3u; k < s_max_iter; ++k)
	{
	  term *= -z * std::pow(Tp(k) / Tp(k - 1), k - 2);
	  _W += term;
	  if (std::abs(term) < s_eps * std::abs(_W))
	    break;
	}
      return _W;
    }

  /**
   * This is the second-order Newton root finding algorithm for Lambert W.
   */
  template<typename Tp>
    Tp
    lambert_w_newton(Tp z, Tp _W = Tp{1})
    {
      const auto s_eps = emsr::epsilon(z);
      const auto s_max_iter = 1000u;

      auto wk = _W;
      for (auto k = 0u; k < s_max_iter; ++k)
	{
          const auto expwk = std::exp(wk);
          const auto wexpwk = wk * expwk;
	  const auto wkp1 = wk - (wexpwk - z)
				   / (Tp{1} + wk) / expwk;
	  const auto del = std::abs(wkp1 - wk);
	  wk = wkp1;
	  if (del < s_eps)
	    break;
	}
      return wk;
    }

  /**
   * This is the third-order Halley root finding algorithm for Lambert W.
   */
  template<typename Tp>
    Tp
    lambert_w_halley(Tp z, Tp _W = Tp{1})
    {
      const auto s_eps = emsr::epsilon(z);
      const auto s_max_iter = 1000u;

      auto wk = _W;
      for (auto k = 0u; k < s_max_iter; ++k)
	{
          const auto expwk = std::exp(wk);
	  const auto fact = wk * expwk - z;
          const auto wkp1 = wk - fact
		      / ((wk + 1) * expwk - (wk + 2) * fact / (2 * wk + 2));
	  const auto del = std::abs(wkp1 - wk);
	  wk = wkp1;
	  if (del < s_eps)
	    break;
	}
      return wk;
    }


  /**
   * This is the fifth-order Schroder's update for Lambert W.
   */
  template<typename Tp>
    Tp
    lambert_w_schroder(Tp z, Tp _W)
    {
      const auto y = z * std::exp(-_W);
      const auto f0 = _W - y;
      const auto f1 = Tp{1} + y;
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
template<typename Tp>
  Tp
  lambert_w_0_log_series(Tp z)
  {
    const auto xi = std::log(z);
    const auto lnxi = std::log(xi);
    return xi - lnxi * (Tp{1} - (Tp{1} / xi)
		* (Tp{1} - (Tp{1} / xi) * (Tp{1} - lnxi / Tp{2})));
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
template<typename Tp>
  Tp
  lambert_w_1_log_series(Tp z)
  {
    const auto eta = std::log(-Tp{1} / z);
    const auto lneta = std::log(eta);
    return -eta - lneta * (Tp{1} + (Tp{1} / eta)
		* (Tp{1} + (Tp{1} / eta) * (Tp{1} + lneta / Tp{2})));
  }


template<typename Tp>
  void
  test_lambert_w(Tp proto = Tp{})
  {
    std::cout.precision(emsr::digits10(proto));
    auto w = std::cout.precision() + 8;
    std::cout << std::showpoint << std::scientific;

    const int N0 = 100;
    const auto s_e = emsr::e_v<Tp>;
    const auto s_1_e = Tp{1} / s_e;
    const auto del0 = (s_e + s_1_e) / N0;

    const int Nm1 = 50;
    const auto delm1 = s_1_e / Nm1;

    std::cout << '\n';
    for (int i = 0; i <= N0; ++i)
      {
	auto z = -s_1_e  + del0 * i;
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
	auto z = -s_1_e  + delm1 * i;
        auto W_newton = lambert_w_newton(z, Tp{-2});
        auto W_halley = lambert_w_halley(z, Tp{-2});
	std::cout << ' ' << std::setw(w) << z
		  << ' ' << std::setw(w) << W_newton
		  << ' ' << std::setw(w) << W_halley
		  << '\n';
      }

    std::cout << '\n';
    std::cout << '\n';
    auto term = Tp{-1};
    for (auto k = 3u; k < 50; ++k)
      {
	term *= -std::pow(Tp(k) / Tp(k - 1), k - 2);
	std::cout << ' ' << term << '\n';
      }
  }

int
main()
{
  test_lambert_w(1.0);
}
