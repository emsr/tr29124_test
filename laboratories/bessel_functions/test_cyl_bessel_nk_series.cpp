/**
 *
 */

#include <stdexcept>
#include <cmath>
#include <iostream>
#include <iomanip>

#include <emsr/math_constants.h>
#include <emsr/numeric_limits.h>
#include <emsr/special_functions.h>

  template<typename Tp>
    struct bessel_nk_series_t
    {
      Tp _Z_mu;
      Tp _Z_mup1;
    };

  template<typename Tp>
    emsr::gamma_temme_t<Tp>
    gamma_temme(Tp mu)
    {
      using gammat_t = emsr::gamma_temme_t<Tp>;
      const auto s_eps = emsr::epsilon(mu);
      const auto s_gamma_E = emsr::egamma_v<Tp>;

      if (std::abs(mu) < s_eps)
	return gammat_t{mu, Tp{1}, Tp{1}, -s_gamma_E, Tp{1}};
      else
	{
	  Tp gamp, gamm;
	  if (std::real(mu) <= Tp{0})
	    {
	      gamp = emsr::detail::gamma_reciprocal_series(Tp{1} + mu);
	      gamm = -emsr::detail::gamma_reciprocal_series(-mu) / mu;
	    }
	  else
	    {
	      gamp = emsr::detail::gamma_reciprocal_series(mu) / mu;
	      gamm = emsr::detail::gamma_reciprocal_series(Tp{1} - mu);
	    }
	  const auto gam1 = (gamm - gamp) / (Tp{2} * mu);
	  const auto gam2 = (gamm + gamp) / Tp{2};
	  return gammat_t{mu, gamp, gamm, gam1, gam2};
	}
    }

  template<typename Tp>
    bessel_nk_series_t<Tp>
    old_n(Tp nu, Tp x, int max_iter = 10000)
    {
      const auto s_eps = emsr::epsilon<Tp>();
      const auto s_pi = emsr::pi_v<Tp>;
      const int n = std::nearbyint(nu);
      const auto mu = nu - Tp(n);
      const auto mu2 = mu * mu;
      const auto xi = Tp{1} / x;
      const auto xi2 = Tp{2} * xi;
      const auto x2 = x / Tp{2};
      const auto pimu = s_pi * mu;
      const auto fact = (std::abs(pimu) < s_eps
			? Tp{1}
			: pimu / std::sin(pimu));
      auto d = -std::log(x2);
      auto e = mu * d;
      const auto fact2 = (std::abs(e) < s_eps
			 ? Tp{1}
			 : std::sinh(e) / e);
      const auto gamt = gamma_temme(mu);
      auto ff = (Tp{2} / s_pi) * fact
		* (gamt.gamma_1_value * std::cosh(e)
		 + gamt.gamma_2_value * fact2 * d);
      e = std::exp(e);
      auto p = e / (s_pi * gamt.gamma_plus_value);
      auto q = Tp{1} / (e * s_pi * gamt.gamma_minus_value);
      const auto pimu2 = pimu / Tp{2};
      const auto fact3 = (std::abs(pimu2) < s_eps
			 ? Tp{1} : std::sin(pimu2) / pimu2 );
      const auto r = s_pi * pimu2 * fact3 * fact3;
      auto c = Tp{1};
      d = -x2 * x2;
      auto sum = ff + r * q;
      auto sum1 = p;
      int i;
      for (i = 1; i <= max_iter; ++i)
	{
	  ff = (i * ff + p + q) / (i * i - mu2);
	  c *= d / Tp(i);
	  p /= Tp(i) - mu;
	  q /= Tp(i) + mu;
	  const auto del = c * (ff + r * q);
	  sum += del;
	  const auto del1 = c * p - Tp(i) * del;
	  sum1 += del1;
	  if (std::abs(del) < s_eps * (Tp{1} + std::abs(sum)))
	    break;
	}
      if (i > max_iter)
	throw std::runtime_error("cyl_bessel_nk_series: Series failed to converge");

      auto _Nmu = -sum;
      auto _Nnu1 = -sum1 * xi2;

      return {_Nmu, _Nnu1};
    }

  template<typename Tp>
    bessel_nk_series_t<Tp>
    old_k(Tp nu, Tp x, int max_iter = 10000)
    {
      const auto s_eps = emsr::epsilon<Tp>();
      const auto s_pi = emsr::pi_v<Tp>;
      const int n = std::nearbyint(nu);
      const auto mu = nu - Tp(n);
      const auto mu2 = mu * mu;
      const auto xi = Tp{1} / x;
      const auto xi2 = Tp{2} * xi;
      const auto x2 = x / Tp{2};
      const auto pimu = s_pi * mu;
      const auto fact = (std::abs(pimu) < s_eps
			? Tp{1}
			: pimu / std::sin(pimu));
      auto d = -std::log(x2);
      auto e = mu * d;
      const auto fact2 = (std::abs(e) < s_eps
			 ? Tp{1}
			 : std::sinh(e) / e);
      const auto gamt = gamma_temme(mu);
      auto ff = fact
		* (gamt.gamma_1_value * std::cosh(e)
		 + gamt.gamma_2_value * fact2 * d);
      auto sum = ff;
      e = std::exp(e);
      auto p = e / (Tp{2} * gamt.gamma_plus_value);
      auto q = Tp{1} / (Tp{2} * e * gamt.gamma_minus_value);
      auto c = Tp{1};
      d = x2 * x2;
      auto sum1 = p;
      int i;
      for (i = 1; i <= max_iter; ++i)
	{
	  ff = (i * ff + p + q) / (i * i - mu2);
	  c *= d / Tp(i);
	  p /= Tp(i) - mu;
	  q /= Tp(i) + mu;
	  const auto del = c * ff;
	  sum += del;
	  const auto del1 = c * (p - Tp(i) * ff);
	  sum1 += del1;
	  if (std::abs(del) < s_eps * std::abs(sum))
	    break;
	}
      if (i > max_iter)
	throw std::runtime_error("cyl_bessel_ik_steed: K-series failed to converge");
      auto _Kmu = sum;
      auto _Kmu1 = sum1 * xi2;

      return {_Kmu, _Kmu1};
    }

  /**
   * This routine computes the dominant cylindrical bessel function solutions
   * by series summation for order @f$ |\mu| < 1/2 @f$.
   *
   * @param mu The order of the Bessel functions @f$ |\mu| < 1/2 @f$.
   * @param x  The argument of the Bessel functions.
   * @param modified If true solve for the modified Bessel function
   *                   @f$ K_\mu @f$, otherwise solve for the Neumann function
   *                   @f$ N_\mu @f$,
   * @return A structure containing Z_{\mu} and Z_{\mu+1}.
   */
  template<typename Tp>
    bessel_nk_series_t<Tp>
    cyl_bessel_nk_series(Tp mu, Tp x, bool modified = false,
			   int max_iter = 100)
    {
      const auto s_eps = emsr::epsilon<Tp>();
      const auto s_pi = emsr::pi_v<Tp>;
      const auto xi = Tp{1} / x;
      const auto x2 = x / Tp{2};

      const auto fact = Tp{1} / emsr::detail::sinc_pi(mu);
      const auto lx2 = -std::log(x2);
      const auto arg = mu * lx2;
      const auto fact2 = emsr::detail::sinhc(arg);
      const auto gamt = emsr::detail::gamma_temme(mu);
      const auto norm = modified ? Tp{-1} : Tp{2} / s_pi;
      auto ff = norm * fact
		* (gamt.gamma_1_value * std::cosh(arg)
		 + gamt.gamma_2_value * fact2 * lx2);
      const auto e = std::exp(arg);
      auto p = norm * e / (Tp{2} * gamt.gamma_plus_value);
      auto q = norm / (e * Tp{2} * gamt.gamma_minus_value);
      const auto fact3 = modified
			 ? Tp{0}
			 : emsr::detail::sinc_pi(mu / Tp{2});
      const auto r = modified
		     ? Tp{0}
		     : fact3 * fact3 * s_pi * s_pi * mu / Tp{2};
      auto c = Tp{1};
      const auto d = modified ? x2 * x2 : -x2 * x2;
      auto sum_mu = ff + r * q;
      auto sum_mup1 = p;
      int i;
      for (i = 1; i <= max_iter; ++i)
	{
	  ff = (i * ff + p + q)
	       / ((Tp(i) - mu) * (Tp(i) + mu));
	  c *= d / Tp(i);
	  p /= Tp(i) - mu;
	  q /= Tp(i) + mu;
	  const auto del_mu = c * (ff + r * q);
	  sum_mu += del_mu;
	  const auto del_mup1 = c * p - Tp(i) * del_mu;
	  sum_mup1 += del_mup1;
	  if (std::abs(del_mu) < s_eps * std::abs(sum_mu))
	    break;
	}
      if (i > max_iter)
	throw std::runtime_error("cyl_bessel_nk_series: Series failed to converge");
      auto _N_mu = -sum_mu;
      auto _N_mup1 = -Tp{2} * xi * sum_mup1;

      return {_N_mu, _N_mup1};
    }

template<typename Tp>
  void
  test_cyl_bessel_nk_series()
  {
    const auto p = std::numeric_limits<Tp>::digits10;
    std::cout.precision(p);
    const auto w = 8 + std::cout.precision();

    std::cout << "\ncyl_neumann\n";
    std::cout << ' ' << std::setw(w) << "x"
	      << ' ' << std::setw(w) << "new N_mu"
	      << ' ' << std::setw(w) << "new N_mup1"
	      << ' ' << std::setw(w) << "old N_mu"
	      << ' ' << std::setw(w) << "old N_mup1"
	      << ' ' << std::setw(w) << "new N_mup1 - std"
	      << ' ' << std::setw(w) << "N_mu (new-lab)/lab"
	      << ' ' << std::setw(w) << "N'_mu (new-lab)/lab"
	      << '\n';
    for (auto nu : {Tp{0}, Tp{1}/Tp{3}, Tp{1}/Tp{2}})
      {
	for (int i = 1; i < 20; ++i)
	  {
	    const auto x = i * 0.1;
	    const auto jn = emsr::detail::cyl_bessel_jn(nu, x);
	    const auto np1 = emsr::cyl_neumann(nu + 1, x);
	    const auto nn = cyl_bessel_nk_series(nu, x);
	    const auto Np_mu = nu * nn._Z_mu / x - nn._Z_mup1;
	    const auto no = old_n(nu, x);
	    std::cout << ' ' << std::setw(w) << x
		      << ' ' << std::setw(w) << nn._Z_mu
		      << ' ' << std::setw(w) << nn._Z_mup1
		      << ' ' << std::setw(w) << no._Z_mu
		      << ' ' << std::setw(w) << no._Z_mup1
		      << ' ' << std::setw(w) << nn._Z_mup1 - np1
		      << ' ' << std::setw(w) << (nn._Z_mu - jn.N_value) / jn.N_value
		      << ' ' << std::setw(w) << (Np_mu - jn.N_deriv) / jn.N_deriv
		      << '\n';
	  }
      }

    std::cout << "\ncyl_bessel_k\n";
    std::cout << ' ' << std::setw(w) << "x"
	      << ' ' << std::setw(w) << "new K_mu"
	      << ' ' << std::setw(w) << "new K_mup1"
	      << ' ' << std::setw(w) << "old K_mu"
	      << ' ' << std::setw(w) << "old K_mup1"
	      << ' ' << std::setw(w) << "new K_mup1 - std"
	      << ' ' << std::setw(w) << "K_mu (new-lab)/lab"
	      << ' ' << std::setw(w) << "K'_mu (new-lab)/lab"
	      << '\n';
    for (auto nu : {Tp{0}, Tp{1}/Tp{3}, Tp{1}/Tp{2}})
      {
	for (int i = 1; i < 20; ++i)
	  {
	    const auto x = i * 0.1;
	    const auto ik = emsr::detail::cyl_bessel_ik(nu, x);
	    const auto kp1 = emsr::cyl_bessel_k(nu + 1, x);
	    const auto kn = cyl_bessel_nk_series(nu, x, true);
	    const auto Np_mu = nu * kn._Z_mu / x - kn._Z_mup1;
	    const auto ko = old_k(nu, x);
	    std::cout << ' ' << std::setw(w) << x
		      << ' ' << std::setw(w) << kn._Z_mu
		      << ' ' << std::setw(w) << kn._Z_mup1
		      << ' ' << std::setw(w) << ko._Z_mu
		      << ' ' << std::setw(w) << ko._Z_mup1
		      << ' ' << std::setw(w) << kn._Z_mup1 - kp1
		      << ' ' << std::setw(w) << (kn._Z_mu - ik.K_value) / ik.K_value
		      << ' ' << std::setw(w) << (Np_mu - ik.K_deriv) / ik.K_deriv
		      << '\n';
	  }
      }
  }


int
main()
{
  test_cyl_bessel_nk_series<double>();

  return 0;
}
