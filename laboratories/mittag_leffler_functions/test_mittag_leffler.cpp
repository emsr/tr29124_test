/**
 *
 */

#include <iostream>
#include <iomanip>
#include <algorithm> // max({...})
#include <cmath>

#include <emsr/numeric_limits.h>
#include <emsr/float128_io.h>
#include <emsr/special_functions.h>

#include <emsr/integration.h>

  /* Monotone integrand for the Mittag-Leffler function. */
  template<typename Tp>
    std::complex<Tp>
    mittag_leffler_K(Tp alpha, Tp beta, Tp chi,
		       const std::complex<Tp>& z)
    {
      const auto s_pi = emsr::pi_v<Tp>;
      const auto chip1 = std::pow(chi, Tp{1} / alpha);
      const auto chip2 = std::pow(chi, (Tp{1} - beta) / alpha);
      return chip2
	   * std::exp(-chip1)
	   * (chi * std::sin(s_pi * (Tp{1} - beta))
	      - z * std::sin(s_pi * (Tp{1} - beta + alpha)))
	   / (chi * chi - Tp{2} * chi * z * std::cos(alpha * s_pi)
		 + z * z)
	   / s_pi / alpha;
    }

  /* Monotone integral for the Mittag-Leffler function. */
  template<typename Tp>
    std::complex<Tp>
    mittag_leffler_K_integral(Tp alpha, Tp beta,
				Tp chi_min, Tp chi_max,
				const std::complex<Tp>& z)
    {
      const auto s_eps = emsr::epsilon(chi_min);
      auto func = [alpha, beta, z](Tp chi)
		    -> std::complex<Tp>
		    { return mittag_leffler_K(alpha, beta, chi, z); };

      const auto epsabs = Tp{100} * s_eps;
      const auto epsrel = Tp{0};
      auto ws = emsr::cquad_workspace<Tp, std::complex<Tp>>();

      auto quad
	= emsr::cquad_integrate(ws, func, chi_min, chi_max,
				     epsabs, epsrel);

      return quad.result;
    }

  /* Oscillatory integrand for the Mittag-Leffler function. */
  template<typename Tp>
    std::complex<Tp>
    mittag_leffler_P(Tp alpha, Tp beta, Tp epsilon, Tp phi,
		       const std::complex<Tp>& z)
    {
      const auto s_i = std::complex<Tp>{0, 1};
      const auto s_pi = emsr::pi_v<Tp>;
      const auto epsp1 = std::pow(epsilon, Tp{1} / alpha);
      const auto rat = Tp{1} + (Tp{1} - beta) / alpha;
      const auto epsp2 = std::pow(epsilon, rat);
      const auto omega = phi * rat + epsp1 * std::sin(phi / alpha);
      return epsp2
	   * std::exp(epsp1 * std::cos(phi / alpha))
	   * std::polar(Tp{1}, omega)
	   / (epsilon * s_i - z)
	   / Tp{2} / s_pi / alpha;
    }

  /* Oscillatory integral for the Mittag-Leffler function. */
  template<typename Tp>
    std::complex<Tp>
    mittag_leffler_P_integral(Tp alpha, Tp beta, Tp epsilon,
				Tp phi_min, Tp phi_max,
				const std::complex<Tp>& z)
    {
      const auto s_eps = emsr::epsilon(phi_min);
      auto func = [alpha, beta, epsilon, z](Tp phi)
		    -> std::complex<Tp>
		    {
		      return mittag_leffler_P(alpha, beta,
						epsilon, phi, z);
		    };

      const auto epsabs = Tp{100} * s_eps;
      const auto epsrel = Tp{0};
      auto ws = emsr::cquad_workspace<Tp, std::complex<Tp>>();

      auto quad
	= emsr::cquad_integrate(ws, func, phi_min, phi_max,
				     epsabs, epsrel);

      return quad.result;
    }


  /**
   * Compute the Mittag-Leffer function:
   * @f[
   *   E_{\alpha,\beta}(z) = \sum_{k=0}^{\infty}\frac{z^k}{\beta + \alpha k},
   *   \mbox{  } \alpha > 0, \beta \elem \complex, z \elem \complex
   * @f]
   *
   * @see COMPUTATION OF THE MITTAG-LEFFLER FUNCTION @f$ E_{\alpha,\beta}(z) @f$
   * AND ITS DERIVATIVE, Rudolf Gorenflo, Joulia Loutchko & Yuri Luchko
   */
  template<typename Tp>
    std::complex<Tp>
    mittag_leffler(Tp alpha, Tp beta, const std::complex<Tp>& z)
    {
      using Cmplx = std::complex<Tp>;
      const auto s_eps = emsr::epsilon(alpha);
      const auto s_2pi = emsr::tau_v<Tp>;
      const auto s_pi = emsr::pi_v<Tp>;

      const auto az = std::abs(z);
      if (alpha > Tp{1})
	{
          unsigned int k0 = Tp{1} + std::floor(alpha);
	  const auto alpha0 = alpha / k0;
	  const auto rho0 = std::pow(z, Tp{1} / Tp(k0));
	  const auto lamb = s_2pi / Tp(k0);

	  auto E = Cmplx{0};
	  for (auto k = 0u; k < k0; ++k)
	    {
	      auto zk = rho0 * std::polar(Tp{1}, lamb * Tp(k));
	      E += mittag_leffler(alpha0, beta, zk);
	    }
	  return E / Tp(k0);
	}
      else if (az < s_eps)
	return emsr::detail::gamma_reciprocal(beta);
      else if (az < Tp{1})
	{
	  unsigned int k0 = std::max(std::ceil((Tp{1} - beta) / alpha),
				std::ceil(std::log(s_eps * (Tp{1} - az))
					    / std::log(az)));
	  auto E = Cmplx{0};
	  auto zk = Cmplx{1};
	  for (auto k = 0u; k <= k0; ++k)
	    {
	      const auto arg = beta + alpha * k;
	      const auto term = zk
			* emsr::detail::gamma_reciprocal(arg);
	      E += term;
	      if (std::abs(term) < s_eps)
		break;
	      zk *= z;
	    }
	  return E;
	}
      else if (az > std::floor(Tp{10} + Tp{5} * alpha))
	{
	  unsigned int k0 = std::floor(-std::log(s_eps) / std::log(az));
	  auto E = Cmplx{0};
	  auto zk = Cmplx{1};
	  for (auto k = 1u; k <= k0; ++k)
	    {
	      zk /= z;
	      E += zk * emsr::detail::gamma_reciprocal(beta - alpha * k);
	    }
	  if (std::arg(z)
	      < s_pi * (alpha / Tp{4} + std::min(Tp{1}, alpha) / Tp{2}))
	    {
	      const auto zp1 = std::pow(z, Tp{1} / alpha);
	      const auto zp2 = std::pow(z, (Tp{1} - beta) / alpha);
	      const auto extra = zp2 * std::exp(zp1) / alpha;
	      return extra - E;
	    }
	  else
	    return -E;
	}
      else
	{
	  auto chi0 = Tp{0};
	  if (beta >= Tp{0})
	    chi0 = std::max({Tp{1}, Tp{2} * az,
			std::pow(-std::log(s_pi * s_eps / Tp{6}), alpha)});
	  else
	    {
	      const auto abeta = std::abs(beta);
	      chi0 = std::max({std::pow(Tp{1} + abeta, alpha),
			 Tp{2} * az, 
		std::pow(-Tp{2}
			  * std::log(s_pi * s_eps
			  / (Tp{6} * (abeta + Tp{2})
			   * std::pow(Tp{2} * abeta, abeta))),
			 alpha)});
	    }

	  const auto absarz = std::abs(std::arg(z));
	  if (absarz > alpha * s_pi + s_eps)
	    {
	      if (beta <= Tp{1})
		return mittag_leffler_K_integral(alpha, beta,
						   Tp{0}, chi0, z);
	      else
		{
		  const auto api = s_pi * alpha;
		  return mittag_leffler_K_integral(alpha, beta,
						     Tp{1}, chi0, z)
		       + mittag_leffler_P_integral(alpha, beta, Tp{1},
						     -api, api, z);
		}
	    }
	  else if (absarz < alpha * s_pi - s_eps)
	    {
	      const auto zp1 = std::pow(z, Tp{1} / alpha);
	      const auto zp2 = std::pow(z, (Tp{1} - beta) / alpha);
	      const auto extra = zp2 * std::exp(zp1) / alpha;
	      if (beta <= Tp{1})
		return mittag_leffler_K_integral(alpha, beta,
						   Tp{0}, chi0, z)
			+ extra;
	      else
		{
		  const auto lo = az / Tp{2};
		  const auto api = s_pi * alpha;
		  return mittag_leffler_K_integral(alpha, beta,
						     lo, chi0, z)
		       + mittag_leffler_P_integral(alpha, beta, lo,
							-api, api, z)
		       + extra;
		}
	    }
	  else
	    {
	      const auto lo = (az + Tp{1}) / Tp{2};
	      const auto api = s_pi * alpha;
	      return mittag_leffler_K_integral(alpha, beta,
						 lo, chi0, z)
		   + mittag_leffler_P_integral(alpha, beta, lo,
						 -api, api, z);
	    }
	}
    }

  /**
   * Compute the derivative of the Mittag-Leffer function:
   * @f[
   *   E_{\alpha,\beta}(z) = \sum_{k=0}^{\infty}
   *                       \frac{z^k}{\Gamma(\beta + \alpha k)},
   *   \mbox{  } \alpha > 0, \beta \elem \complex, z \elem \complex
   * @f]
   *
   * @see COMPUTATION OF THE MITTAG-LEFFLER FUNCTION @f$ E_{\alpha,\beta}(z) @f$
   * AND ITS DERIVATIVE, Rudolf Gorenflo, Joulia Loutchko & Yuri Luchko
   */
  template<typename Tp>
    std::complex<Tp>
    mittag_leffler_deriv(Tp alpha, Tp beta,
			   const std::complex<Tp>& z)
    {
      using Cmplx = std::complex<Tp>;
      const auto s_eps = emsr::epsilon(alpha);

      const auto az = std::abs(z);
      if (az < Tp{1})
	{
	  auto k1 = Tp{0};
	  if (alpha > Tp{1})
	    k1 = Tp{1} + (Tp{2} - alpha - beta) / (alpha - Tp{1});
	  else
	    {
	      const auto D = Tp{1}
			     + alpha * (alpha - Tp{4} * beta + Tp{6});
	      const auto omega = alpha + beta - Tp{3} / Tp{2};
	      const auto rat = Tp{1} + (Tp{3} - alpha - beta) / alpha;
	      if (D <= Tp{0})
		k1 = rat;
	      else
		k1 = std::max(rat,
			Tp{1}
			+ (Tp{1} - Tp{2} * omega * alpha + std::sqrt(D))
				  / (2 * alpha * alpha));
	    }
	  k1 = std::ceil(k1);
	  unsigned int k0 = std::max(k1,
				 std::ceil(std::log(s_eps * (Tp{1} - az))
					 / std::log(az)));
	  auto Ep = Cmplx{0};
	  auto zk = Cmplx{1};
	  for (auto k = 0u; k <= k0; ++k)
	    {
	      Ep += Tp(k + 1) * zk
		    * emsr::detail::gamma_reciprocal(beta
						      + alpha * Tp(k + 1));
	      zk *= z;
	    }
	  return Ep;
	}
      else
	return (mittag_leffler(alpha, beta - Tp{1}, z)
	      - (beta - Tp{1}) * mittag_leffler(alpha, beta, z))
	     / alpha / z;
    }

template<typename Tp>
  void
  test_mittag_leffler(Tp proto = Tp{})
  {
    using namespace std::complex_literals;
    using Cmplx = std::complex<Tp>;

    std::cout.precision(emsr::digits10(proto));
    auto width = std::cout.precision() + 8;
    std::cout << std::showpoint << std::scientific;

    // Figure 1
    {
      const auto alpha = Tp{1} / Tp{4};
      const auto beta = Tp{1};
      std::cout << '\n';
      std::cout << '\n';
      const auto del = Tp{1} / Tp{10};
      for (int i = 0; i <= 100; ++i)
	{
	  auto t = i * del;
	  auto ml_val = mittag_leffler(alpha, beta, Cmplx(-t, 0));
	  auto ml_der = -mittag_leffler_deriv(alpha, beta, Cmplx(-t, 0));
	  std::cout << std::setw(width) << t
		    << std::setw(width) << std::real(ml_val)
		    << std::setw(width) << std::real(ml_der)
		    << '\n';
	}
      std::cout << std::flush;
    }

    // Figure 2
    {
      const auto alpha = Tp{7} / Tp{4};
      const auto beta = Tp{1};
      std::cout << '\n';
      std::cout << '\n';
      const auto del = Tp{1} / Tp{10};
      for (int i = 0; i <= 500; ++i)
	{
	  auto t = i * del;
	  auto ml_val = mittag_leffler(alpha, beta, Cmplx(-t, 0));
	  auto ml_der = -mittag_leffler_deriv(alpha, beta, Cmplx(-t, 0));
	  std::cout << std::setw(width) << t
		    << std::setw(width) << std::real(ml_val)
		    << std::setw(width) << std::real(ml_der)
		    << '\n';
	}
      std::cout << std::flush;
    }

    // Figure 3
    {
      const auto alpha = Tp{9} / Tp{4};
      const auto beta = Tp{1};
      std::cout << '\n';
      std::cout << '\n';
      const auto del = Tp{1} / Tp{10};
      for (int i = 0; i <= 1000; ++i)
	{
	  auto t = i * del;
	  auto ml_val = mittag_leffler(alpha, beta, Cmplx(-t, 0));
	  auto ml_der = -mittag_leffler_deriv(alpha, beta, Cmplx(-t, 0));
	  std::cout << std::setw(width) << t
		    << std::setw(width) << std::real(ml_val)
		    << std::setw(width) << std::real(ml_der)
		    << '\n';
	}
      std::cout << std::flush;
    }

    // Figure 4
    {
      const auto alpha = Tp{3} / Tp{4};
      const auto beta = Tp{1};
      const auto s_pi = emsr::pi_v<Tp>;
      const auto phase = alpha * s_pi / Tp{4};
      std::cout << '\n';
      std::cout << '\n';
      const auto del = Tp{1} / Tp{10};
      for (int i = 0; i <= 50; ++i)
	{
	  auto z = std::polar(i * del, phase);
	  auto ml_val = mittag_leffler(alpha, beta, z);
	  std::cout << std::setw(width) << std::abs(z)
		    << std::setw(width) << std::abs(ml_val)
		    << std::setw(width) << std::real(ml_val)
		    << std::setw(width) << std::imag(ml_val)
		    << '\n';
	}
      std::cout << std::flush;
    }

    // Figure 5
    {
      const auto alpha = Tp{3} / Tp{4};
      const auto beta = Tp{1};
      const auto s_pi = emsr::pi_v<Tp>;
      const auto phase = alpha * s_pi / Tp{2};
      std::cout << '\n';
      std::cout << '\n';
      const auto del = Tp{1} / Tp{10};
      for (int i = 0; i <= 500; ++i)
	{
	  auto z = std::polar(i * del, phase);
	  auto ml_val = mittag_leffler(alpha, beta, z);
	  std::cout << std::setw(width) << std::abs(z)
		    << std::setw(width) << std::abs(ml_val)
		    << std::setw(width) << std::real(ml_val)
		    << std::setw(width) << std::imag(ml_val)
		    << '\n';
	}
      std::cout << std::flush;
    }

    // Figure 6
    {
      const auto alpha = Tp{3} / Tp{4};
      const auto beta = Tp{1};
      const auto s_pi = emsr::pi_v<Tp>;
      const auto phase = Tp{3} * alpha * s_pi / Tp{4};
      std::cout << '\n';
      std::cout << '\n';
      const auto del = Tp{1} / Tp{10};
      for (int i = 0; i <= 500; ++i)
	{
	  auto z = std::polar(i * del, phase);
	  auto ml_val = mittag_leffler(alpha, beta, z);
	  std::cout << std::setw(width) << std::abs(z)
		    << std::setw(width) << std::abs(ml_val)
		    << std::setw(width) << std::real(ml_val)
		    << std::setw(width) << std::imag(ml_val)
		    << '\n';
	}
      std::cout << std::flush;
    }

    // Figure 7
    {
      const auto alpha = Tp{3} / Tp{4};
      const auto beta = Tp{1};
      const auto s_pi = emsr::pi_v<Tp>;
      const auto phase = s_pi;
      std::cout << '\n';
      std::cout << '\n';
      const auto del = Tp{1} / Tp{10};
      for (int i = 0; i <= 200; ++i)
	{
	  auto z = std::polar(i * del, phase);
	  auto ml_val = mittag_leffler(alpha, beta, z);
	  std::cout << std::setw(width) << std::abs(z)
		    << std::setw(width) << std::abs(ml_val)
		    << std::setw(width) << std::real(ml_val)
		    << std::setw(width) << std::imag(ml_val)
		    << '\n';
	}
      std::cout << std::flush;
    }

    // Figure 8
    {
      const auto alpha = Tp{5} / Tp{4};
      const auto beta = Tp{1};
      const auto s_pi = emsr::pi_v<Tp>;
      const auto phase = alpha * s_pi / Tp{4};
      std::cout << '\n';
      std::cout << '\n';
      const auto del = Tp{1} / Tp{10};
      for (int i = 0; i <= 100; ++i)
	{
	  auto z = std::polar(i * del, phase);
	  auto ml_val = mittag_leffler(alpha, beta, z);
	  std::cout << std::setw(width) << std::abs(z)
		    << std::setw(width) << std::abs(ml_val)
		    << std::setw(width) << std::real(ml_val)
		    << std::setw(width) << std::imag(ml_val)
		    << '\n';
	}
      std::cout << std::flush;
    }

    // Figure 9
    {
      const auto alpha = Tp{5} / Tp{4};
      const auto beta = Tp{1};
      const auto s_pi = emsr::pi_v<Tp>;
      const auto phase = alpha * s_pi / Tp{2};
      std::cout << '\n';
      std::cout << '\n';
      const auto del = Tp{1} / Tp{10};
      for (int i = 0; i <= 500; ++i)
	{
	  auto z = std::polar(i * del, phase);
	  auto ml_val = mittag_leffler(alpha, beta, z);
	  std::cout << std::setw(width) << std::abs(z)
		    << std::setw(width) << std::abs(ml_val)
		    << std::setw(width) << std::real(ml_val)
		    << std::setw(width) << std::imag(ml_val)
		    << '\n';
	}
      std::cout << std::flush;
    }

    // Figure 10
    {
      const auto alpha = Tp{5} / Tp{4};
      const auto beta = Tp{1};
      const auto s_pi = emsr::pi_v<Tp>;
      const auto phase = Tp{3} * alpha * s_pi / Tp{4};
      std::cout << '\n';
      std::cout << '\n';
      const auto del = Tp{1} / Tp{10};
      for (int i = 0; i <= 500; ++i)
	{
	  auto z = std::polar(i * del, phase);
	  auto ml_val = mittag_leffler(alpha, beta, z);
	  std::cout << std::setw(width) << std::abs(z)
		    << std::setw(width) << std::abs(ml_val)
		    << std::setw(width) << std::real(ml_val)
		    << std::setw(width) << std::imag(ml_val)
		    << '\n';
	}
      std::cout << std::flush;
    }

    // Figure 11
    {
      const auto alpha = Tp{5} / Tp{4};
      const auto beta = Tp{1};
      const auto s_pi = emsr::pi_v<Tp>;
      const auto phase =  s_pi;
      std::cout << '\n';
      std::cout << '\n';
      const auto del = Tp{1} / Tp{10};
      for (int i = 0; i <= 1000; ++i)
	{
	  auto z = std::polar(i * del, phase);
	  auto ml_val = mittag_leffler(alpha, beta, z);
	  std::cout << std::setw(width) << std::abs(z)
		    << std::setw(width) << std::abs(ml_val)
		    << std::setw(width) << std::real(ml_val)
		    << std::setw(width) << std::imag(ml_val)
		    << '\n';
	}
      std::cout << std::flush;
    }
  }

int
main()
{
  test_mittag_leffler(1.0);
}
