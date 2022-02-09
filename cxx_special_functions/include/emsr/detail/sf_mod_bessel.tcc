
// Copyright (C) 2006-2019 Free Software Foundation, Inc.
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

/** @file bits/sf_mod_bessel.tcc
 *  This is an internal header file, included by other library headers.
 *  Do not attempt to use it directly. @headername{cmath}
 */

//
// ISO C++ 14882 TR29124: Mathematical Special Functions
//

// Written by Edward Smith-Rowland.
//
// References:
// (1) Handbook of Mathematical Functions,
//     Ed. Milton Abramowitz and Irene A. Stegun,
//     Dover Publications,
//     Section 9, pp. 355-434, Section 10 pp. 435-478
// (2) The Gnu Scientific Library, http://www.gnu.org/software/gsl
// (3) Numerical Recipes in C, by W. H. Press, S. A. Teukolsky,
//     W. T. Vetterling, B. P. Flannery, Cambridge University Press (1992),
//     2nd ed, pp. 246-249.

#ifndef SF_MOD_BESSEL_TCC
#define SF_MOD_BESSEL_TCC 1

#include <stdexcept>
#include <utility> // For exchange.
#include <complex>

#include <emsr/fp_type_util.h>
#include <emsr/sf_trig.h> // sin_pi, cos_pi, polar_pi
#include <emsr/sf_bessel.h> // detail::cyl_bessel_ij_series, cyl_bessel_jn

namespace emsr
{
namespace detail
{

  /**
   * @brief This routine computes the asymptotic modified cylindrical
   * 	    Bessel and functions of order nu: @f$ I_{\nu}(x) @f$,
   * 	    @f$ N_{\nu}(x) @f$.  Use this for @f$ x >> \nu^2 + 1 @f$.
   *
   * References:
   *  (1) Handbook of Mathematical Functions,
   * 	  ed. Milton Abramowitz and Irene A. Stegun,
   * 	  Dover Publications,
   * 	  Section 9 p. 364, Equations 9.2.5-9.2.10
   *
   * @param  nu  The order of the Bessel functions.
   * @param  x   The argument of the Bessel functions.
   * @return A struct containing the modified cylindrical Bessel functions
   *         of the first and second kinds and their derivatives.
   */
  template<typename Tnu, typename Tp>
    constexpr emsr::cyl_mod_bessel_t<Tnu, Tp, Tp>
    cyl_bessel_ik_scaled_asymp(Tnu nu, Tp x)
    {
      // FIXME: This will promote float to double if Tnu is integral.
      using Val = emsr::fp_promote_t<Tnu, Tp>;
      using Real = emsr::num_traits_t<Val>;
      using bess_t = emsr::cyl_mod_bessel_t<Tnu, Tp, Tp>;
      const auto s_pi = emsr::pi_v<Real>;

      const auto sums = cyl_bessel_asymp_sums(nu, x, +1);

      const auto coef = std::sqrt(Real{1} / (Real{2} * s_pi * x));
      return bess_t{nu, x,
		      coef * (sums.Psum - sums.Qsum),
		      coef * (sums.Rsum - sums.Ssum),
		      coef * s_pi * (sums.Psum + sums.Qsum),
		     -coef * s_pi * (sums.Rsum + sums.Ssum)};
    }

  /**
   * @param  nu  The order of the Bessel functions.
   * @param  x   The argument of the Bessel functions.
   * @param  do_scaled  If true, scale I, I' by exp(-x) and K, K' by exp(+x).
   */
  template<typename Tnu, typename Tp>
    constexpr emsr::cyl_mod_bessel_t<Tnu, Tp, Tp>
    cyl_bessel_ik_asymp(Tnu nu, Tp x, bool do_scaled = false)
    {
      using bess_t = emsr::cyl_mod_bessel_t<Tnu, Tp, Tp>;

      auto ik = cyl_bessel_ik_scaled_asymp(nu, x);

      if (do_scaled)
	return ik;
      else
	{
	  const auto exp = std::exp(x);
	  const auto iexp = Tp{1} / exp;
	  // @todo Check for over/under-flow in cyl_bessel_ik_asymp.
	  return bess_t{ik.nu_arg, ik.x_arg,
			  exp * ik.I_value, exp * ik.I_deriv,
			  iexp * ik.K_value, iexp * ik.K_deriv};
	}
    }

  /**
   * @brief  Compute the modified Bessel functions @f$ I_\nu(x) @f$ and
   * 	     @f$ K_\nu(x) @f$ and their first derivatives
   * 	     @f$ I'_\nu(x) @f$ and @f$ K'_\nu(x) @f$ respectively.
   * 	     These four functions are computed together for numerical
   * 	     stability.
   *
   * @param  nu  The order of the Bessel functions.
   * @param  x   The argument of the Bessel functions.
   * @param  do_scaled  If true, scale I, I' by exp(-x) and K, K' by exp(+x).
   * @return A struct containing the modified cylindrical Bessel functions
   *         of the first and second kinds and their derivatives.
   */
  template<typename Tp>
    emsr::cyl_mod_bessel_t<Tp, Tp, Tp>
    cyl_bessel_ik_steed(Tp nu, Tp x, bool do_scaled = false)
    {
      using bess_t = emsr::cyl_mod_bessel_t<Tp, Tp, Tp>;
      const auto s_eps = emsr::epsilon(x);
      const auto s_pi = emsr::pi_v<Tp>;
      const auto s_fp_min = Tp{10} * s_eps;
      constexpr int s_max_iter = 15000;
      const auto s_x_min = Tp{2};

      const int n = std::nearbyint(nu);

      const auto mu = nu - Tp(n);
      const auto mu2 = mu * mu;
      const auto xi = Tp{1} / x;
      const auto xi2 = Tp{2} * xi;
      auto h = std::max(s_fp_min, nu * xi);
      auto b = xi2 * nu;
      auto d = Tp{0};
      auto c = h;
      int i;
      for (i = 1; i <= s_max_iter; ++i)
	{
	  b += xi2;
	  d = Tp{1} / (b + d);
	  c = b + Tp{1} / c;
	  const auto del = c * d;
	  h *= del;
	  if (std::abs(del - Tp{1}) < s_eps)
	    break;
	}
      if (i > s_max_iter)
	return cyl_bessel_ik_asymp(nu, x, do_scaled);

      auto Inul = s_fp_min;
      auto Ipnul = h * Inul;
      auto Inul1 = Inul;
      auto Ipnu1 = Ipnul;
      auto fact = nu * xi;
      for (int l = n; l >= 1; --l)
	{
	  const auto Inutemp = fact * Inul + Ipnul;
	  fact -= xi;
	  Ipnul = fact * Inutemp + Inul;
	  Inul = Inutemp;
	}

      const auto f = Ipnul / Inul;
      bool scaled = false;
      Tp Kmu, Knu1;
      if (x < s_x_min)
	{
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
	  for (i = 1; i <= s_max_iter; ++i)
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
	  if (i > s_max_iter)
	    throw std::runtime_error("cyl_bessel_ik_steed: K-series failed to converge");
	  Kmu = sum;
	  Knu1 = sum1 * xi2;
	}
      else
	{
	  scaled = true;
	  auto b = Tp{2} * (Tp{1} + x);
	  auto d = Tp{1} / b;
	  auto delh = d;
	  auto h = delh;
	  auto q1 = Tp{0};
	  auto q2 = Tp{1};
	  const auto a1 = Tp{0.25L} - mu2;
	  auto q = c = a1;
	  auto a = -a1;
	  auto s = Tp{1} + q * delh;
	  int i;
	  for (i = 2; i <= s_max_iter; ++i)
	    {
	      a -= Tp(2 * (i - 1));
	      c = -a * c / i;
	      const auto qnew = (q1 - b * q2) / a;
	      q1 = q2;
	      q2 = qnew;
	      q += c * qnew;
	      b += Tp{2};
	      d = Tp{1} / (b + a * d);
	      delh = (b * d - Tp{1}) * delh;
	      h += delh;
	      const auto dels = q * delh;
	      s += dels;
	      if (std::abs(dels / s) < s_eps)
		break;
	    }
	  if (i > s_max_iter)
	    throw std::runtime_error("cyl_bessel_ik_steed: Steed's method failed");
	  h = a1 * h;
	  // We are scaling this branch to prevent under/overflow. Removing...
	  // * std::exp(-x)
	  Kmu = std::sqrt(s_pi / (Tp{2} * x)) / s;
	  Knu1 = Kmu * (mu + x + Tp{0.5L} - h) * xi;
	}

      auto Kpmu = mu * xi * Kmu - Knu1;
      auto Inumu = xi / (f * Kmu - Kpmu);
      auto Inu = Inumu * Inul1 / Inul;
      auto Ipnu = Inumu * Ipnu1 / Inul;
      for (int i = 1; i <= n; ++i)
	Kmu = std::exchange(Knu1, (mu + Tp(i)) * xi2 * Knu1 + Kmu);
      auto Knu = Kmu;
      auto Kpnu = nu * xi * Kmu - Knu1;

      if (do_scaled && !scaled)
	{
	  const auto exp = std::exp(x);
	  const auto iexp = Tp{1} / exp;
	  Inu *= iexp;
	  Ipnu *= iexp;
	  Knu *= exp;
	  Kpnu *= exp;
	}
      else if (!do_scaled && scaled)
	{
	  const auto exp = std::exp(x);
	  const auto iexp = Tp{1} / exp;
	  /// @todo Check for over/underflow for large-argument modified Bessel.
	  Inu *= exp;
	  Ipnu *= exp;
	  Knu *= iexp;
	  Kpnu *= iexp;
	}

      return bess_t{nu, x, Inu, Ipnu, Knu, Kpnu};
    }

  /**
   * @brief  Return the modified cylindrical Bessel functions
   *         and their derivatives of order @f$ \nu @f$ by various means.
   *
   * @param  nu  The order of the Bessel functions.
   * @param  x   The argument of the Bessel functions.
   * @param  do_scaled  If true, scale I, I' by exp(-x) and K, K' by exp(+x).
   * @return A struct containing the modified cylindrical Bessel functions
   *         of the first and second kinds and their derivatives.
   */
  template<typename Tp>
    emsr::cyl_mod_bessel_t<Tp, Tp, Tp>
    cyl_bessel_ik(Tp nu, Tp x, bool do_scaled = false)
    {
      using bess_t = emsr::cyl_mod_bessel_t<Tp, Tp, Tp>;
      const auto s_eps = emsr::epsilon(x);
      const auto s_inf = emsr::infinity(x);
      const auto s_pi = emsr::pi_v<Tp>;
      if (nu < Tp{0})
	{
	  const auto Bessm = cyl_bessel_ik(-nu, x, do_scaled);
	  const auto sinnupi = emsr::sin_pi(-nu);
	  if (std::abs(sinnupi) < s_eps) // Carefully preserve +-inf.
	    return bess_t{nu, x, Bessm.I_value, Bessm.I_deriv,
					Bessm.K_value, Bessm.K_deriv};
	  else
	    return bess_t{nu, x,
	      Bessm.I_value + Tp{2} * sinnupi * Bessm.K_value / s_pi,
	      Bessm.I_deriv + Tp{2} * sinnupi * Bessm.K_deriv / s_pi,
	      Bessm.K_value, Bessm.K_deriv};
	}
      else if (x == Tp{0})
	{
	  Tp Inu, Ipnu;
	  if (nu == Tp{0})
	    {
	      Inu = Tp{1};
	      Ipnu = Tp{0};
	    }
	  else if (nu == Tp{1})
	    {
	      Inu = Tp{0};
	      Ipnu = Tp{0.5L};
	    }
	  else
	    {
	      Inu = Tp{0};
	      Ipnu = Tp{0};
	    }
	  return bess_t{nu, x, Inu, Ipnu, s_inf, -s_inf};
	}
      else if (x > Tp{1000})
	return cyl_bessel_ik_asymp(nu, x, do_scaled);
      else
	return cyl_bessel_ik_steed(nu, x, do_scaled);
    }

  /**
   * @brief  Return the regular modified Bessel function of order
   * 	     @f$ \nu @f$: @f$ I_{\nu}(x) @f$.
   *
   * The regular modified cylindrical Bessel function is:
   * @f[
   *  I_{\nu}(x) = \sum_{k=0}^{\infty}
   * 		\frac{(x/2)^{\nu + 2k}}{k!\Gamma(\nu+k+1)}
   * @f]
   *
   * @param  nu  The order of the regular modified Bessel function.
   * @param  x   The argument of the regular modified Bessel function.
   * @return  The output regular modified Bessel function.
   */
  template<typename Tp>
    Tp
    cyl_bessel_i(Tp nu, Tp x)
    {
      if (x < Tp{0})
	throw std::domain_error("cyl_bessel_i: Argument < 0");
      else if (std::isnan(nu) || std::isnan(x))
	return emsr::quiet_NaN(x);
      else if (nu >= Tp{0} && x * x < Tp{10} * (nu + Tp{1}))
	return detail::cyl_bessel_ij_series(nu, x, +1, 200);
      else
	return cyl_bessel_ik(nu, x).I_value;
    }

  template<typename Tp>
    Tp
    cyl_bessel_i_scaled(Tp nu, Tp x)
    {
      if (x < Tp{0})
	throw std::domain_error("cyl_bessel_i: Argument < 0");
      else if (std::isnan(nu) || std::isnan(x))
	return emsr::quiet_NaN(x);
      else
	return cyl_bessel_ik(nu, x, true).I_value;
    }

  /**
   * @brief  Return the irregular modified Bessel function
   * 	     @f$ K_{\nu}(x) @f$ of order @f$ \nu @f$.
   *
   * The irregular modified Bessel function is defined by:
   * @f[
   * 	K_{\nu}(x) = \frac{\pi}{2}
   * 		     \frac{I_{-\nu}(x) - I_{\nu}(x)}{\sin \nu\pi}
   * @f]
   * where for integral @f$ \nu = n @f$ a limit is taken:
   * @f$ lim_{\nu \to n} @f$.
   * For negative argument we have simply:
   * @f[
   * 	K_{-\nu}(x) = K_{\nu}(x)
   * @f]
   *
   * @param  nu  The order of the irregular modified Bessel function.
   * @param  x   The argument of the irregular modified Bessel function.
   * @return  The output irregular modified Bessel function.
   */
  template<typename Tp>
    Tp
    cyl_bessel_k(Tp nu, Tp x)
    {
      if (x < Tp{0})
	throw std::domain_error("cyl_bessel_k: Argument < 0");
      else if (std::isnan(nu) || std::isnan(x))
	return emsr::quiet_NaN(x);
      else
	return cyl_bessel_ik(nu, x).K_value;
    }

  template<typename Tp>
    Tp
    cyl_bessel_k_scaled(Tp nu, Tp x)
    {
      if (x < Tp{0})
	throw std::domain_error("cyl_bessel_k_scaled: Argument < 0");
      else if (std::isnan(nu) || std::isnan(x))
	return emsr::quiet_NaN(x);
      else
	return cyl_bessel_ik(nu, x, true).K_value;
    }

  /**
   * @brief  Compute the spherical modified Bessel functions
   * 	     @f$ i_n(x) @f$ and @f$ k_n(x) @f$ and their first
   * 	     derivatives @f$ i'_n(x) @f$ and @f$ k'_n(x) @f$
   * 	     respectively.
   *
   * @param  n  The order of the modified spherical Bessel function.
   * @param  x  The argument of the modified spherical Bessel function.
   * @return A struct containing the modified spherical Bessel functions
   *         of the first and second kinds and their derivatives.
   */
  template<typename Tp>
    emsr::sph_mod_bessel_t<unsigned int, Tp, Tp>
    sph_bessel_ik(unsigned int n, Tp x)
    {
      using sph_t = emsr::sph_mod_bessel_t<unsigned int, Tp, Tp>;
      const auto s_NaN = emsr::quiet_NaN(x);

      if (std::isnan(x))
	return sph_t{n, x, s_NaN, s_NaN, s_NaN, s_NaN};
      else if (x == Tp{0})
	{
	  const auto s_inf = emsr::infinity(x);
	  if (n == 0)
	    return sph_t{n, x, Tp{1}, Tp{0}, s_inf, -s_inf};
	  else
	    return sph_t{n, x, Tp{0}, Tp{0}, s_inf, -s_inf};
	}
      else
	{
	  const auto nu = Tp(n + 0.5L);
	  auto Bess = cyl_bessel_ik(nu, x);

	  const auto factor = (emsr::sqrtpi_v<Tp> / emsr::sqrt2_v<Tp>)
			      / std::sqrt(x);

	  const auto i_n = factor * Bess.I_value;
	  const auto ip_n = factor * Bess.I_deriv
			    - i_n / (Tp{2} * x);
	  const auto k_n = factor * Bess.K_value;
	  const auto kp_n = factor * Bess.K_deriv
			    - k_n / (Tp{2} * x);

	  return sph_t{n, x, i_n, ip_n, k_n, kp_n};
	}
    }


  /**
   * @brief  Compute the Airy functions
   * 	     @f$ Ai(x) @f$ and @f$ Bi(x) @f$ and their first
   * 	     derivatives @f$ Ai'(x) @f$ and @f$ Bi(x) @f$
   * 	     respectively.
   *
   * @param  z  The argument of the Airy functions.
   * @param  do_scaled  If true, scale @f$ Ai @f$, @f$ Ai' @f$ by @f$ \exp(+\xi) @f$
   *                    and @f$ Bi @f$, @f$ Bi' @f$ by @f$ \exp(-\xi) @f$
   *                    where @f$ \xi = 2 x^{3/2} / 3 @f$.
   * @return A struct containing the Airy functions
   *         of the first and second kinds and their derivatives.
   */
  template<typename Tp>
    emsr::airy_t<Tp, Tp>
    airy(Tp z, bool do_scaled = false)
    {
      using ai_t = emsr::airy_t<Tp, Tp>;
      const auto s_NaN = emsr::quiet_NaN(z);
      const auto s_inf = emsr::infinity(z);
      const auto s_pi = emsr::pi_v<Tp>;
      const auto s_sqrt3 = emsr::sqrt3_v<Tp>;
      const auto absz = std::abs(z);
      const auto rootz = std::sqrt(absz);
      const auto xi = Tp{2} * absz * rootz / Tp{3};

      if (std::isnan(z))
	return ai_t{z, s_NaN, s_NaN, s_NaN, s_NaN};
      else if (z == s_inf)
	return ai_t{z, Tp{0}, Tp{0}, s_inf, s_inf};
      else if (z == -s_inf)
	return ai_t{z, Tp{0}, Tp{0}, Tp{0}, Tp{0}};
      else if (z > Tp{0})
	{
          // I and K have the opposite scaling factor so propagating scaling can't work.
	  const auto Bess13 = cyl_bessel_ik(Tp{1} / Tp{3}, xi, false/*do_scaled*/);
	  const auto Ai = rootz * Bess13.K_value / (s_sqrt3 * s_pi);
	  const auto Bi = rootz * (Bess13.K_value / s_pi
				    + Tp{2} * Bess13.I_value / s_sqrt3);

	  const auto Bess23 = cyl_bessel_ik(Tp{2} / Tp{3}, xi, do_scaled);
	  const auto Aip = -z * Bess23.K_value / (s_sqrt3 * s_pi);
	  const auto Bip = z * (Bess23.K_value / s_pi
				 + Tp{2} * Bess23.I_value / s_sqrt3);

          auto exp = Tp{1};
          if (do_scaled)
            exp = std::exp(xi);

	  return ai_t{z, Ai * exp, Aip * exp, Bi / exp, Bip / exp};
	}
      else if (z < Tp{0})
	{
	  const auto Bess13 = emsr::detail::cyl_bessel_jn(Tp{1} / Tp{3}, xi);
	  const auto Ai = +rootz * (Bess13.J_value
				     - Bess13.N_value / s_sqrt3) / Tp{2};
	  const auto Bi = -rootz * (Bess13.N_value
				     + Bess13.J_value / s_sqrt3) / Tp{2};

	  const auto Bess23 = emsr::detail::cyl_bessel_jn(Tp{2} / Tp{3}, xi);
	  const auto Aip = absz * (Bess23.N_value / s_sqrt3
				    + Bess23.J_value) / Tp{2};
	  const auto Bip = absz * (Bess23.J_value / s_sqrt3
				    - Bess23.N_value) / Tp{2};

          auto exp = Tp{1};
          if (do_scaled)
            exp = std::exp(xi);

	  return ai_t{z, Ai * exp, Aip * exp, Bi / exp, Bip / exp};
	}
      else
	{
	  // Reference:
	  //  Abramowitz & Stegun, page 446 section 10.4.4 on Airy functions.
	  // The number is Ai(0) = 3^{-2/3}/\Gamma(2/3).
	  const auto Ai
	    = Tp{0.3550280538878172392600631860041831763979791741991772L};
	  const auto Bi = Ai * s_sqrt3;

	  // Reference:
	  //  Abramowitz & Stegun, page 446 section 10.4.5 on Airy functions.
	  // The number is Ai'(0) = -3^{-1/3}/\Gamma(1/3).
	  const auto Aip
	    = -Tp{0.25881940379280679840518356018920396347909113835493L};
	  const auto Bip = -Aip * s_sqrt3;

	  return ai_t{z, Ai, Aip, Bi, Bip};
	}
    }

  /**
   * @brief  Compute the Fock-type Airy functions
   * 	     @f$ w_1(x) @f$ and @f$ w_2(x) @f$ and their first
   * 	     derivatives @f$ w_1'(x) @f$ and @f$ w_2'(x) @f$
   * 	     respectively.
   * @f[
   *   w_1(x) = \sqrt{\pi}(Ai(x) + iBi(x))
   * @f]
   * @f[
   *   w_2(x) = \sqrt{\pi}(Ai(x) - iBi(x))
   * @f]
   *
   * @param  x   The argument of the Airy functions.
   * @return A struct containing the Fock-type Airy functions
   *         of the first and second kinds and their derivatives.
   */
  template<typename Tp>
    emsr::fock_airy_t<Tp, std::complex<Tp>>
    fock_airy(Tp x)
    {
      using Cmplx = std::complex<Tp>;
      using fock_t = emsr::fock_airy_t<Tp, Cmplx>;
      const auto s_sqrtpi = emsr::sqrtpi_v<Tp>;

      const auto Ai = airy(x);

      const auto w1 = s_sqrtpi * Cmplx(Ai.Ai_value, Ai.Bi_value);
      const auto w1p = s_sqrtpi * Cmplx(Ai.Ai_deriv, Ai.Bi_deriv);
      const auto w2 = s_sqrtpi * Cmplx(Ai.Ai_value, -Ai.Bi_value);
      const auto w2p = s_sqrtpi * Cmplx(Ai.Ai_deriv, -Ai.Bi_deriv);

      return fock_t{x, w1, w1p, w2, w2p};
    }

} // namespace detail
} // namespace std

#endif // SF_MOD_BESSEL_TCC
