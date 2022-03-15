
// Copyright (C) 2019 Free Software Foundation, Inc.
// Copyright (C) 2020-2022 Edward M. Smith-Rowland
//
// This program is free software; you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation; either version 3 of the License, or (at
// your option) any later version.
//
// This library is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// Under Section 7 of GPL version 3, you are granted additional
// permissions described in the GCC Runtime Library Exception, version
// 3.1, as published by the Free Software Foundation.

// You should have received a copy of the GNU General Public License and
// a copy of the GCC Runtime Library Exception along with this program;
// see the files COPYING3 and COPYING.RUNTIME respectively.  If not, see
// <http://www.gnu.org/licenses/>.

/** @file bits/sf_mittag_leffler.tcc
 *  This is an internal header file, included by other library headers.
 *  Do not attempt to use it directly. @headername{cmath}
 */

#ifndef SF_MITTAG_LEFFLER_TCC
#define SF_MITTAG_LEFFLER_TCC 1

#include <complex>

#include <emsr/integration.h>
#include <emsr/sf_gamma.h>

namespace emsr
{
namespace detail
{

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
		      return mittag_leffler_P(alpha, beta, epsilon, phi, z);
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
   *   E_{\alpha,\beta}(z)
   *     = \sum_{k=0}^\infty \frac{z^k}{\Gamma(\beta + \alpha k)},
   *   \mbox{  } \alpha > 0, \beta \in \$\mathbb{C}\$, z \in \$\mathbb{C}\$
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
      if (alpha == Tp{0})
	{
	  if (az < Tp{1})
	    return gamma_reciprocal(beta) / (Tp{1} - z);
	  else
	    {
	      const auto s_inf = emsr::infinity(alpha);
	      return std::complex<Tp>{s_inf, s_inf};
	    }
	}
      else if (alpha < Tp{0})
	return gamma_reciprocal(beta)
	     - mittag_leffler(-alpha, beta, Tp{1} / z);
      else if (alpha > Tp{1})
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
	return gamma_reciprocal(beta);
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
	      const auto term = zk * gamma_reciprocal(arg);
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
	      E += zk * gamma_reciprocal(beta - alpha * k);
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
		std::pow(-Tp{2} * std::log(s_pi * s_eps
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
		    * gamma_reciprocal(beta + alpha * Tp(k + 1));
	      zk *= z;
	    }
	  return Ep;
	}
      else
	return (mittag_leffler(alpha, beta - Tp{1}, z)
	      - (beta - Tp{1}) * mittag_leffler(alpha, beta, z))
	     / alpha / z;
    }

} // namespace detail
} // namespace emsr

#endif // SF_MITTAG_LEFFLER_TCC
