
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

/** @file bits/sf_bessel.tcc
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
//     ed. Milton Abramowitz and Irene A. Stegun,
//     Dover Publications,
//     Section 9, pp. 355-434, Section 10 pp. 435-478
// (2) The Gnu Scientific Library, http://www.gnu.org/software/gsl
// (3) Numerical Recipes in C, by W. H. Press, S. A. Teukolsky,
//     W. T. Vetterling, B. P. Flannery, Cambridge University Press (1992),
//     2nd ed, pp. 240-245

#ifndef SF_BESSEL_TCC
#define SF_BESSEL_TCC 1

#include <stdexcept>
#include <complex>
#include <utility> // For exchange.

#include <emsr/fp_type_util.h>
#include <emsr/sf_trig.h> // sin_pi, cos_pi, polar_pi
#include <emsr/sf_gamma.h> // log_gamma, gamma_reciprocal_series
#include <emsr/specfun_state.h>

namespace emsr
{
namespace detail
{

  /**
   * @brief This routine returns the cylindrical Bessel functions
   * 	    of order @f$ \nu @f$: @f$ J_{\nu}(z) @f$ or @f$ I_{\nu}(z) @f$
   * 	    by series expansion.
   *
   * The modified cylindrical Bessel function is:
   * @f[
   *  Z_{\nu}(x) = \sum_{k=0}^{\infty}
   * 		\frac{\sigma^k (x/2)^{\nu + 2k}}{k!\Gamma(\nu+k+1)}
   * @f]
   * where @f$ \sigma = +1 @f$ or @f$ -1 @f$ for
   * @f$ Z = I @f$ or @f$ J @f$ respectively.
   *
   * See Abramowitz & Stegun, 9.1.10
   * 	 Abramowitz & Stegun, 9.6.7
   *  (1) Handbook of Mathematical Functions,
   * 	  ed. Milton Abramowitz and Irene A. Stegun,
   * 	  Dover Publications,
   * 	  Equation 9.1.10 p. 360 and Equation 9.6.10 p. 375
   *
   * @param  nu  The order of the Bessel function.
   * @param  x   The argument of the Bessel function.
   * @param  sgn  The sign of the alternate terms
   * 		    -1 for the Bessel function of the first kind.
   * 		    +1 for the modified Bessel function of the first kind.
   * @param  max_iter  The maximum number of iterations for sum.
   * @return  The output Bessel function.
   */
  template<typename Tnu, typename Tp>
    constexpr Tp
    cyl_bessel_ij_series(Tnu nu, Tp x, int sgn,
			   unsigned int max_iter)
    {
      // FIXME: This will promote float to double if Tnu is integral.
      using Val = emsr::fp_promote_t<Tnu, Tp>;
      using Real = emsr::num_traits_t<Val>;
      const auto s_eps = emsr::epsilon<Real>();
      if (std::abs(x) < s_eps)
	{
	  if (nu == Tnu{0})
	    return Tp{1};
	  else
	    return Tp{0};
	}
      else
	{
	  const auto x2 = x / Real{2};

	  const auto xx4 = Val(sgn) * x2 * x2;
	  auto Jn = Val{1};
	  auto term = Val{1};
	  for (unsigned int i = 1; i < max_iter; ++i)
	    {
	      term *= xx4 / (Val(i) * (Val(nu) + Val(i)));
	      Jn += term;
	      if (std::abs(term / Jn) < s_eps)
		break;
	    }

	  auto fact = Val(nu) * std::log(x2);
	  fact -= emsr::detail::log_gamma(Real{1} + nu);
	  fact = std::exp(fact);
	  return fact * Jn;
	}
    }

  /**
   * A type for Bessel asymptotic sums.
   */
  template<typename Tnu, typename Tp>
    struct cyl_bessel_asymp_sums_t
    {
      // FIXME: This will promote float to double if Tnu is integral.
      using Val = emsr::fp_promote_t<Tnu, Tp>;
      Val Psum;
      Val Qsum;
      Val Rsum;
      Val Ssum;
    };

  /**
   * @brief This routine computes the asymptotic cylindrical Bessel
   * 	    and Neumann functions of order nu: @f$ J_{\nu}(z) @f$,
   * 	    @f$ N_{\nu}(z) @f$.  Use this for @f$ z >> \nu^2 + 1 @f$.
   *
   * @f[
   *   J_{\nu}(z) = \left(\frac{2}{\pi z}\right)^{1/2} \left(
   *    \cos(\omega)\sum_{k=0}^{\infty}(-1)^k\frac{a_{2k}(\nu)}{z^{2k}}
   *  - \sin(\omega)\sum_{k=0}^{\infty}(-1)^k\frac{a_{2k+1}(\nu)}{z^{2k+1}}
   *    \right)
   * @f]
   * and
   * @f[
   *   N_{\nu}(z) = \left(\frac{2}{\pi z}\right)^{1/2} \left(
   *    \sin(\omega)\sum_{k=0}^{\infty}(-1)^k\frac{a_{2k}(\nu)}{z^{2k}}
   *  + \cos(\omega)\sum_{k=0}^{\infty}(-1)^k\frac{a_{2k+1}(\nu)}{z^{2k+1}}
   *    \right)
   * @f]
   * where @f$ \omega = z - \nu\pi/2 - \pi/4 @f$ and
   * @f[
   *   a_{k}(\nu) = \frac{(4\nu^2 - 1^2)(4\nu^2 - 3^2)...(4\nu^2 - (2k-1)^2)}
   *                     {8^k k!}
   * @f]
   * The derivatives are also computed:
   * @f[
   *   J'_{\nu}(z) = -\left(\frac{2}{\pi z}\right)^{1/2} \left(
   *    \sin(\omega)\sum_{k=0}^{\infty}(-1)^k\frac{b_{2k}(\nu)}{z^{2k}}
   *  + \cos(\omega)\sum_{k=0}^{\infty}(-1)^k\frac{b_{2k+1}(\nu)}{z^{2k+1}}
   *    \right)
   * @f]
   * and
   * @f[
   *   N'_{\nu}(z) = \left(\frac{2}{\pi z}\right)^{1/2} \left(
   *    \cos(\omega)\sum_{k=0}^{\infty}(-1)^k\frac{b_{2k}(\nu)}{z^{2k}}
   *  - \sin(\omega)\sum_{k=0}^{\infty}(-1)^k\frac{b_{2k+1}(\nu)}{z^{2k+1}}
   *    \right)
   * @f]
   * where
   * @f[
   *   b_{k}(\nu) = \frac{(4\nu^2 - 1^2)(4\nu^2 - 3^2)...(4\nu^2 - (2k-3)^2)}
   *                     {8^k k!}
   * @f]
   *
   * These sums work everywhere but on the negative real axis:
   * @f$ |ph(z)| < \pi - \delta @f$.
   *
   * References:
   *  (1) Handbook of Mathematical Functions,
   * 	  ed. Milton Abramowitz and Irene A. Stegun,
   * 	  Dover Publications,
   * 	  Section 9 p. 364, Equations 9.2.5-9.2.10
   *  (2) Digital Library of Mathematical Fundtions
   *      https://dlmf.nist.gov/10.17
   *
   * @param  nu  The order of the Bessel functions.
   * @param  x   The argument of the Bessel functions.
   * @return A struct containing the cylindrical Bessel functions
   *         of the first and second kinds and their derivatives.
   */
  template<typename Tnu, typename Tp>
    constexpr cyl_bessel_asymp_sums_t<Tnu, Tp>
    cyl_bessel_asymp_sums(Tnu nu, Tp x, int sgn)
    {
      // FIXME: This will promote float to double if Tnu is integral.
      using Val = emsr::fp_promote_t<Tnu, Tp>;
      using Real = emsr::num_traits_t<Val>;
      using bess_t = cyl_bessel_asymp_sums_t<Tnu, Tp>;
      const auto s_eps = emsr::epsilon<Real>();
      const auto __2nu = Real{2} * nu;
      const auto __4nu2 = __2nu * __2nu;
      const auto r8x = Tp{1} / (Real{8} * x);
      const auto nu_min = std::real(nu / Real{2});
      const auto nu_max = std::abs(Real{100} * (nu + Tnu{1}));
      auto k = 0;
      auto bk_xk = Val{1};
      auto Rsum = bk_xk;
      auto ak_xk = Val{1};
      auto Psum = ak_xk;
      auto convP = false;
      ++k;
      auto __2km1 = 1;
      bk_xk *= (__4nu2 + 3) * r8x;
      auto Ssum = bk_xk;
      ak_xk *= (__2nu - __2km1) * (__2nu + __2km1) * r8x;
      auto Qsum = ak_xk;
      auto convQ = false;
      auto ak_xk_prev = std::abs(ak_xk);
      do
	{
	  ++k;
	  auto rk8x = r8x / Real(k);
	  __2km1 += 2;
	  bk_xk = sgn * (__4nu2 + __2km1 * (__2km1 + 2)) * ak_xk * rk8x;
	  Rsum += bk_xk;
	  ak_xk *= sgn * (__2nu - __2km1) * (__2nu + __2km1) * rk8x;
	  if (k > nu_min && std::abs(ak_xk) > ak_xk_prev)
	    break;
	  Psum += ak_xk;
	  ak_xk_prev = std::abs(ak_xk);
	  convP = std::abs(ak_xk) < s_eps * std::abs(Psum);

	  ++k;
	  rk8x = r8x / Real(k);
	  __2km1 += 2;
	  bk_xk = (__4nu2 + __2km1 * (__2km1 + 2)) * ak_xk * rk8x;
	  Ssum += bk_xk;
	  ak_xk *= (__2nu - __2km1) * (__2nu + __2km1) * rk8x;
	  if (k > nu_min && std::abs(ak_xk) > ak_xk_prev)
	    break;
	  Qsum += ak_xk;
	  ak_xk_prev = std::abs(ak_xk);
	  convQ = std::abs(ak_xk) < s_eps * std::abs(Qsum);

	  if (convP && convQ)
	    break;
	}
      while (k < nu_max);

      return bess_t{Psum, Qsum, Rsum, Ssum};
    }

  /**
   *
   */
  template<typename Tnu, typename Tp>
    constexpr emsr::cyl_bessel_t<Tnu, Tp, Tp>
    cyl_bessel_jn_asymp(Tnu nu, Tp x)
    {
      // FIXME: This will promote float to double if Tnu is integral.
      using Val = emsr::fp_promote_t<Tnu, Tp>;
      using Real = emsr::num_traits_t<Val>;
      using bess_t = emsr::cyl_bessel_t<Tnu, Tp, Tp>;
      const auto s_pi = emsr::pi_v<Real>;
      const auto s_pi_2 = emsr::pi_v<Real> / Real{2};

      const auto sums = cyl_bessel_asymp_sums(nu, x, -1);

      const auto omega = x - (nu + Real{0.5L}) * s_pi_2;
      const auto c = std::cos(omega);
      const auto s = std::sin(omega);

      const auto coef = std::sqrt(Real{2} / (s_pi * x));
      return bess_t{nu, x,
		 coef * (c * sums.Psum - s * sums.Qsum),
		-coef * (s * sums.Rsum + c * sums.Ssum),
		 coef * (s * sums.Psum + c * sums.Qsum),
		 coef * (c * sums.Rsum - s * sums.Ssum)};
    }

  /**
   * @brief Compute the gamma functions required by the Temme series
   * 	    expansions of @f$ N_\nu(x) @f$ and @f$ K_\nu(x) @f$.
   * @f[
   *   \Gamma_1 = \frac{1}{2\mu}
   * 	 \left[\frac{1}{\Gamma(1 - \mu)} - \frac{1}{\Gamma(1 + \mu)}\right]
   * @f]
   * and
   * @f[
   *   \Gamma_2 = \frac{1}{2}
   *     \left[\frac{1}{\Gamma(1 - \mu)} + \frac{1}{\Gamma(1 + \mu)}\right]
   * @f]
   * where @f$ -1/2 <= \mu <= 1/2 @f$ is @f$ \mu = \nu - N @f$ and @f$ N @f$.
   * is the nearest integer to @f$ \nu @f$.
   * The values of @f$ \Gamma(1 + \mu) @f$ and @f$ \Gamma(1 - \mu) @f$
   * are returned as well.
   *
   * The accuracy requirements on this are exquisite.
   *
   * @param mu The input parameter of the gamma functions.
   * @return  An output structure containing four gamma functions.
   */
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

  /**
   * @brief  Compute the Bessel @f$ J_\nu(x) @f$ and Neumann
   * 	     @f$ N_\nu(x) @f$ functions and their first derivatives
   * 	     @f$ J'_\nu(x) @f$ and @f$ N'_\nu(x) @f$ respectively.
   * 	     These four functions are computed together for numerical
   * 	     stability.
   *
   * @param nu The order of the Bessel functions.
   * @param x  The argument of the Bessel functions.
   * @return A struct containing the cylindrical Bessel functions
   *         of the first and second kinds and their derivatives.
   */
  template<typename Tp>
    emsr::cyl_bessel_t<Tp, Tp, Tp>
    cyl_bessel_jn_steed(Tp nu, Tp x)
    {
      using bess_t = emsr::cyl_bessel_t<Tp, Tp, Tp>;
      const auto s_inf = emsr::infinity(x);
      const auto s_eps = emsr::epsilon(x);
      const auto s_tiny = emsr::lim_min(x);
      const auto s_pi = emsr::pi_v<Tp>;
      // When the multiplier is N i.e.
      // fp_min = N * min()
      // Then J_0 and N_0 tank at x = 8 * N (J_0 = 0 and N_0 = nan)!
      //const Tp s_fp_min = Tp{20} * emsr::lim_min(nu);
      constexpr int s_max_iter = 15000;
      const auto s_x_min = Tp{2};
      const auto s_fp_min = emsr::sqrt_min(nu);

      const int n = (x < s_x_min
		    ? std::nearbyint(nu)
		    : std::max(0,
			       static_cast<int>(nu - x + Tp{1.5L})));

      const auto mu = nu - Tp(n);
      const auto mu2 = mu * mu;
      const auto xi = Tp{1} / x;
      const auto xi2 = Tp{2} * xi;
      const auto Wronski = xi2 / s_pi;
      int isign = 1;
      auto h = std::max(s_fp_min, nu * xi);
      auto b = xi2 * nu;
      auto d = Tp{0};
      auto c = h;
      int i;
      for (i = 1; i <= s_max_iter; ++i)
	{
	  b += xi2;
	  d = b - d;
	  if (std::abs(d) < s_fp_min)
	    d = s_fp_min;
	  d = Tp{1} / d;
	  c = b - Tp{1} / c;
	  if (std::abs(c) < s_fp_min)
	    c = s_fp_min;
	  const auto del = c * d;
	  h *= del;
	  if (d < Tp{0})
	    isign = -isign;
	  if (std::abs(del - Tp{1}) < s_eps)
	    break;
	}
      if (i > s_max_iter)
	return cyl_bessel_jn_asymp(nu, x);

      auto Jnul = isign * s_fp_min;
      auto Jpnul = h * Jnul;
      auto Jnul1 = Jnul;
      auto Jpnu1 = Jpnul;
      auto fact = nu * xi;
      for (int l = n; l >= 1; --l)
	{
	  const auto Jnutemp = fact * Jnul + Jpnul;
	  fact -= xi;
	  Jpnul = fact * Jnutemp - Jnul;
	  Jnul = Jnutemp;
	}
      if (Jnul == Tp{0})
	Jnul = s_eps;

      const auto f = Jpnul / Jnul;
      Tp Nmu, Nnu1, Npmu, Jmu;
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
	  for (i = 1; i <= s_max_iter; ++i)
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
	  if (i > s_max_iter)
	    throw std::runtime_error("cyl_bessel_jn_steed: Y-series failed to converge");
	  Nmu = -sum;
	  Nnu1 = -sum1 * xi2;
	  Npmu = mu * xi * Nmu - Nnu1;
	  Jmu = Wronski / (Npmu - f * Nmu);
	}
      else
	{
	  const auto s_i = std::complex<Tp>{0, 1};
	  auto a = Tp{0.25L} - mu2;
	  auto pq = std::complex<Tp>(-xi / Tp{2}, Tp{1});
	  auto b = std::complex<Tp>(Tp{2} * x, Tp{2});
	  auto fact = a * xi / std::norm(pq);
	  auto c = b + s_i * fact * std::conj(pq);
	  auto d = std::conj(b) / std::norm(b);
	  auto dl = c * d;
	  pq *= dl;
	  int i;
	  for (i = 2; i <= s_max_iter; ++i)
	    {
	      a += Tp(2 * (i - 1));
	      b += s_i * Tp{2};
	      d = a * d + b;
	      if (std::abs(d) < s_fp_min)
		d = s_fp_min;
	      fact = a / std::norm(c);
	      c = b + fact * std::conj(c);
	      if (std::abs(c) < s_fp_min)
		c = s_fp_min;
	      d = std::conj(d) / std::norm(d);
	      dl = c * d;
	      pq *= dl;
	      if (std::abs(dl - Tp{1}) < s_eps)
		break;
	    }
	  if (i > s_max_iter)
	    throw std::runtime_error("cyl_bessel_jn_steed: Lentz's method failed");
	  //const auto [p, q] = pq; // This should be a thing.
	  const auto [p, q] = reinterpret_cast<Tp(&)[2]>(pq);
	  const auto gam = (p - f) / q;
	  Jmu = std::sqrt(Wronski / ((p - f) * gam + q));
	  Jmu = std::copysign(Jmu, Jnul);
	  Nmu = gam * Jmu;
	  Npmu = (p + q / gam) * Nmu;
	  Nnu1 = mu * xi * Nmu - Npmu;
        }
      fact = Jmu / Jnul;
      const auto Jnu = fact * Jnul1;
      const auto Jpnu = fact * Jpnu1;
      if (std::abs(s_pi * x * Jnu / Tp{2}) > s_tiny)
	{
	  for (int i = 1; i <= n; ++i)
	    Nmu = std::exchange(Nnu1, (mu + i) * xi2 * Nnu1 - Nmu);
	  const auto Nnu = Nmu;
	  const auto Npnu = nu * xi * Nmu - Nnu1;
	  return bess_t{nu, x, Jnu, Jpnu, Nnu, Npnu};
	}
      else
	return bess_t{nu, x, Jnu, Jpnu, -s_inf, s_inf};
    }

  /**
   * @brief  Return the cylindrical Bessel functions and their derivatives
   * of order @f$ \nu @f$ by various means.
   */
  template<typename Tp>
    emsr::cyl_bessel_t<Tp, Tp, Tp>
    cyl_bessel_jn(Tp nu, Tp x)
    {
      using bess_t = emsr::cyl_bessel_t<Tp, Tp, Tp>;
      const auto s_eps = emsr::epsilon(x);
      const auto s_inf = emsr::infinity(x);
      if (nu < Tp{0})
	{
	  const auto Bess = cyl_bessel_jn(-nu, x);
	  const auto sinnupi = emsr::sin_pi(-nu);
	  const auto cosnupi = emsr::cos_pi(-nu);
	  if (std::abs(sinnupi) < s_eps)
	    { // Carefully preserve +-inf.
	      const auto sign = std::copysign(Tp{1}, cosnupi);
	      return bess_t{nu, x,
			sign * Bess.J_value, sign * Bess.J_deriv,
			sign * Bess.N_value, sign * Bess.N_deriv};
	    }
	  else if (std::abs(cosnupi) < s_eps)
	    { // Carefully preserve +-inf.
	      const auto sign = std::copysign(Tp{1}, sinnupi);
	      return bess_t{nu, x,
			-sign * Bess.N_value, -sign * Bess.N_deriv,
			 sign * Bess.J_value,  sign * Bess.J_deriv};
	    }
	  else
	    {
	      return bess_t{nu, x,
		cosnupi * Bess.J_value - sinnupi * Bess.N_value,
		cosnupi * Bess.J_deriv - sinnupi * Bess.N_deriv,
		sinnupi * Bess.J_value + cosnupi * Bess.N_value,
		sinnupi * Bess.J_deriv + cosnupi * Bess.N_deriv};
	    }
	}
      else if (x == Tp{0})
	{
	  Tp Jnu, Jpnu;
	  if (nu == Tp{0})
	    {
	      Jnu = Tp{1};
	      Jpnu = Tp{0};
	    }
	  else if (nu == Tp{1})
	    {
	      Jnu = Tp{0};
	      Jpnu = Tp{0.5L};
	    }
	  else
	    {
	      Jnu = Tp{0};
	      Jpnu = Tp{0};
	    }
	  return bess_t{nu, x, Jnu, Jpnu, -s_inf, s_inf};
	}
      else if (x > Tp{1000})
	return cyl_bessel_jn_asymp(nu, x);
      else
	return cyl_bessel_jn_steed(nu, x);
    }

  /**
   * @brief  Return the cylindrical Bessel functions and their derivatives
   *         of real order @f$ \nu @f$ and argument @f$ x < 0 @f$.
   */
  template<typename Tp>
    emsr::cyl_bessel_t<Tp, Tp, std::complex<Tp>>
    cyl_bessel_jn_neg_arg(Tp nu, Tp x)
    {
      using Cmplx = std::complex<Tp>;
      using bess_t = emsr::cyl_bessel_t<Tp, Tp, Cmplx>;
      constexpr Cmplx s_i{0, 1};
      if (x >= Tp{0})
	throw std::domain_error("cyl_bessel_jn_neg_arg: non-negative argument");
      else
	{
	  const auto Bess = cyl_bessel_jn(nu, -x);
	  const auto phm = emsr::polar_pi(Tp{1}, -nu);
	  const auto php = emsr::polar_pi(Tp{1}, nu);
	  const auto cosp = emsr::cos_pi(nu);
	  return bess_t{nu, x,
			  php * Bess.J_value,
			  -php * Bess.J_deriv,
			  phm * Bess.N_value
				+ s_i * Tp{2} * cosp * Bess.J_value,
			  -phm * Bess.N_deriv
				- s_i * Tp{2} * cosp * Bess.J_deriv};
	}
    }


  /**
   * @brief  Return the Bessel function of order @f$ \nu @f$:
   * 	     @f$ J_{\nu}(x) @f$.
   *
   * The cylindrical Bessel function is:
   * @f[
   *  J_{\nu}(x) = \sum_{k=0}^{\infty}
   * 		\frac{(-1)^k (x/2)^{\nu + 2k}}{k!\Gamma(\nu+k+1)}
   * @f]
   *
   * @param  nu  The order of the Bessel function.
   * @param  x   The argument of the Bessel function.
   * @return  The output Bessel function.
   */
  template<typename Tp>
    Tp
    cyl_bessel_j(Tp nu, Tp x)
    {
      if (x < Tp{0})
	throw std::domain_error("cyl_bessel_j: bad argument");
      else if (std::isnan(nu) || std::isnan(x))
	return emsr::quiet_NaN(x);
      else if (nu >= Tp{0} && x * x < Tp{10} * (nu + Tp{1}))
	return cyl_bessel_ij_series(nu, x, -1, 200);
      else
	return cyl_bessel_jn(nu, x).J_value;
    }


  /**
   * @brief  Return the Neumann function of order @f$ \nu @f$:
   * 	     @f$ N_{\nu}(x) @f$.
   *
   * The Neumann function is defined by:
   * @f[
   * 	N_{\nu}(x) = \frac{J_{\nu}(x) \cos \nu\pi - J_{-\nu}(x)}
   * 			  {\sin \nu\pi}
   * @f]
   * where for integral @f$ \nu = n @f$ a limit is taken:
   * @f$ lim_{\nu \to n} @f$.
   *
   * @param  nu  The order of the Neumann function.
   * @param  x   The argument of the Neumann function.
   * @return  The output Neumann function.
   */
  template<typename Tp>
    Tp
    cyl_neumann_n(Tp nu, Tp x)
    {
      if (x < Tp{0})
	throw std::domain_error("cyl_neumann_n: bad argument");
      else if (std::isnan(nu) || std::isnan(x))
	return emsr::quiet_NaN(x);
      else
	return cyl_bessel_jn(nu, x).N_value;
    }

  /**
   * @brief  Return the cylindrical Hankel functions of the first
   *         and second kinds and their derivatives.
   *
   */
  template<typename Tp>
    emsr::cyl_hankel_t<Tp, Tp, std::complex<Tp>>
    cyl_hankel_h1h2(Tp nu, Tp x)
    {
      using Cmplx = std::complex<Tp>;
      constexpr Cmplx s_i{0, 1};

      Cmplx ph1 = Tp{1}, ph2 = Tp{1};
      if (nu < Tp{0})
	{
	  ph1 = emsr::polar_pi(Tp{1}, -nu);
	  ph2 = emsr::polar_pi(Tp{1}, +nu);
	  nu = -nu;
	}

      // The two Bess types are different.
      // We might still be able to assign the real output to the complex one.
      if (x < Tp{0})
	{
	  const auto Bess = cyl_bessel_jn_neg_arg(nu, x);
	  const auto H1 = ph1 * (Bess.J_value + s_i * Bess.N_value);
	  const auto H1p = ph1 * (Bess.J_deriv + s_i * Bess.N_deriv);
	  const auto H2 = ph2 * (Bess.J_value - s_i * Bess.N_value);
	  const auto H2p = ph2 * (Bess.J_deriv - s_i * Bess.N_deriv);
	  return {nu, x, H1, H1p, H2, H2p};
	}
      else
	{
	  const auto Bess = cyl_bessel_jn(nu, x);
	  const auto H1 = ph1 * Cmplx{Bess.J_value, Bess.N_value};
	  const auto H1p = ph1 * Cmplx{Bess.J_deriv, Bess.N_deriv};
	  const auto H2 = ph2 * Cmplx{Bess.J_value, -Bess.N_value};
	  const auto H2p = ph2 * Cmplx{Bess.J_deriv, -Bess.N_deriv};
	  return {nu, x, H1, H1p, H2, H2p};
	}
    }

  /**
   * @brief  Return the cylindrical Hankel function of the first kind
   * 	     @f$ H^{(1)}_\nu(x) @f$.
   *
   * The cylindrical Hankel function of the first kind is defined by:
   * @f[
   *   H^{(1)}_\nu(x) = J_\nu(x) + i N_\nu(x)
   * @f]
   *
   * @param  nu  The order of the spherical Neumann function.
   * @param  x  The argument of the spherical Neumann function.
   * @return  The output spherical Neumann function.
   */
  template<typename Tp>
    std::complex<Tp>
    cyl_hankel_1(Tp nu, Tp x)
    {
      using Cmplx = std::complex<Tp>;
      const auto s_nan = emsr::quiet_NaN(x);
      constexpr Cmplx s_i{0, 1};
      if (nu < Tp{0})
	return emsr::polar_pi(Tp{1}, -nu)
	     * cyl_hankel_1(-nu, x);
      else if (std::isnan(x))
	return Cmplx{s_nan, s_nan};
      else if (x < Tp{0})
	{
	  const auto Bess = cyl_bessel_jn_neg_arg(nu, x);
	  return Bess.J_value + s_i * Bess.N_value;
	}
      else
	{
	  const auto Bess = cyl_bessel_jn(nu, x);
	  return Cmplx{Bess.J_value, Bess.N_value};
	}
    }


  /**
   *   @brief  Return the cylindrical Hankel function of the second kind
   *           @f$ H^{(2)}_nu(x) @f$.
   *
   *   The cylindrical Hankel function of the second kind is defined by:
   *   @f[
   *     H^{(2)}_\nu(x) = J_\nu(x) - i N_\nu(x)
   *   @f]
   *
   *   @param  nu  The order of the spherical Neumann function.
   *   @param  x  The argument of the spherical Neumann function.
   *   @return  The output spherical Neumann function.
   */
  template<typename Tp>
    std::complex<Tp>
    cyl_hankel_2(Tp nu, Tp x)
    {
      using Cmplx = std::complex<Tp>;
      const auto s_nan = emsr::quiet_NaN(x);
      constexpr Cmplx s_i{0, 1};
      if (nu < Tp{0})
	return emsr::polar_pi(Tp{1}, nu)
	     * cyl_hankel_2(-nu, x);
      else if (std::isnan(x))
	return Cmplx{s_nan, s_nan};
      else if (x < Tp{0})
	{
	  const auto Bess = cyl_bessel_jn_neg_arg(nu, x);
	  return Bess.J_value - s_i * Bess.N_value;
	}
      else
	{
	  const auto Bess = cyl_bessel_jn(nu, x);
	  return Cmplx{Bess.J_value, -Bess.N_value};
	}
    }


  /**
   * @brief  Compute the spherical Bessel @f$ j_n(x) @f$
   * 	     and Neumann @f$ n_n(x) @f$ functions and their first
   * 	     derivatives @f$ j_n(x) @f$ and @f$ n'_n(x) @f$
   * 	     respectively.
   *
   * @param  n  The order of the spherical Bessel function.
   * @param  x  The argument of the spherical Bessel function.
   * @return  The output derivative of the spherical Neumann function.
   */
  template<typename Tp>
    emsr::sph_bessel_t<unsigned int, Tp, Tp>
    sph_bessel_jn(unsigned int n, Tp x)
    {
      using bess_t = emsr::sph_bessel_t<unsigned int, Tp, Tp>;
      const auto nu = Tp(n + 0.5L);

      const auto Bess = cyl_bessel_jn(nu, x);

      const auto factor = (emsr::sqrtpi_v<Tp> / emsr::sqrt2_v<Tp>)
			  / std::sqrt(x);

      const auto j_n = factor * Bess.J_value;
      const auto jp_n = factor * Bess.J_deriv - j_n / (Tp{2} * x);
      const auto n_n = factor * Bess.N_value;
      const auto np_n = factor * Bess.N_deriv - n_n / (Tp{2} * x);

      return bess_t{n, x, j_n, jp_n, n_n, np_n};
    }

  /**
   * Return the spherical Bessel functions and their derivatives
   * of order @f$ \nu @f$ and argument @f$ x < 0 @f$.
   */
  template<typename Tp>
    emsr::sph_bessel_t<unsigned int, Tp, std::complex<Tp>>
    sph_bessel_jn_neg_arg(unsigned int n, Tp x)
    {
      using Cmplx = std::complex<Tp>;
      using bess_t = emsr::sph_bessel_t<unsigned int, Tp, Cmplx>;
      if (x >= Tp{0})
	throw std::domain_error("sph_bessel_jn_neg_arg: non-negative argument");
      else
	{
	  const auto nu = Tp(n + 0.5L);
	  const auto Bess = cyl_bessel_jn_neg_arg(nu, x);

	  const auto factor
	    = (emsr::sqrtpi_v<Tp> / emsr::sqrt2_v<Tp>)
	      / std::sqrt(Cmplx(x));

	  const auto j_n = factor * Bess.J_value;
	  const auto jp_n = factor * Bess.J_deriv
			    - j_n / (Tp{2} * x);
	  const auto n_n = factor * Bess.N_value;
	  const auto np_n = factor * Bess.N_deriv
			    - n_n / (Tp{2} * x);

	  return bess_t{n, x, j_n, jp_n, n_n, np_n};
	}
    }


  /**
   * @brief  Return the spherical Bessel function @f$ j_n(x) @f$ of order n
   * and non-negative real argument @c x.
   *
   * The spherical Bessel function is defined by:
   * @f[
   *   j_n(x) = \left(\frac{\pi}{2x} \right) ^{1/2} J_{n+1/2}(x)
   * @f]
   *
   * @param  n  The non-negative integral order
   * @param  x  The non-negative real argument
   * @return  The output spherical Bessel function.
   */
  template<typename Tp>
    Tp
    sph_bessel(unsigned int n, Tp x)
    {
      if (x < Tp{0})
	throw std::domain_error("sph_bessel: bad argument");
      else if (std::isnan(x))
	return emsr::quiet_NaN(x);
      else if (x == Tp{0})
	{
	  if (n == 0)
	    return Tp{1};
	  else
	    return Tp{0};
	}
      else
	return sph_bessel_jn(n, x).j_value;
    }


  /**
   * @brief  Return the spherical Neumann function @f$ n_n(x) @f$ of order n
   * and non-negative real argument @c x.
   *
   * The spherical Neumann function is defined by:
   * @f[
   *  n_n(x) = \left(\frac{\pi}{2x} \right) ^{1/2} N_{n+1/2}(x)
   * @f]
   *
   * @param  n  The order of the spherical Neumann function.
   * @param  x  The argument of the spherical Neumann function.
   * @return  The output spherical Neumann function.
   */
  template<typename Tp>
    Tp
    sph_neumann(unsigned int n, Tp x)
    {
      if (x < Tp{0})
	throw std::domain_error("sph_neumann: bad argument");
      else if (std::isnan(x))
	return emsr::quiet_NaN(x);
      else if (x == Tp{0})
	return -emsr::infinity(x);
      else
	return sph_bessel_jn(n, x).n_value;
    }


  /**
   * @brief  Return the spherical Hankel function of the first kind
   * 	     @f$ h^{(1)}_n(x) @f$.
   *
   * The spherical Hankel function of the first kind is defined by:
   * @f[
   *   h^{(1)}_n(x) = j_n(x) + i n_n(x)
   * @f]
   *
   * @param  n  The order of the spherical Neumann function.
   * @param  x  The argument of the spherical Neumann function.
   * @return  The output spherical Neumann function.
   */
  template<typename Tp>
    std::complex<Tp>
    sph_hankel_1(unsigned int n, Tp x)
    {
      using Cmplx = std::complex<Tp>;
      constexpr Cmplx s_i{0, 1};
      const auto s_nan = emsr::quiet_NaN(x);
      if (std::isnan(x))
	return Cmplx{s_nan, s_nan};
      else if (x < Tp{0})
	{
	  const auto Bess = sph_bessel_jn_neg_arg(n, x);
	  return Bess.j_value + s_i * Bess.n_value;
	}
      else
	{
	  const auto Bess = sph_bessel_jn(n, x);
	  return Cmplx{Bess.j_value, Bess.n_value};
	}
    }


  /**
   * @brief  Return the spherical Hankel function of the second kind
   * 	     @f$ h^{(2)}_n(x) @f$.
   *
   * The spherical Hankel function of the second kind is defined by:
   * @f[
   *   h^{(2)}_n(x) = j_n(x) - i n_n(x)
   * @f]
   *
   * @param  n  The non-negative integral order
   * @param  x  The non-negative real argument
   * @return  The output spherical Neumann function.
   */
  template<typename Tp>
    std::complex<Tp>
    sph_hankel_2(unsigned int n, Tp x)
    {
      using Cmplx = std::complex<Tp>;
      constexpr Cmplx s_i{0, 1};
      const auto s_nan = emsr::quiet_NaN(x);
      if (std::isnan(x))
	return Cmplx{s_nan, s_nan};
      else if (x < Tp{0})
	{
	  const auto Bess = sph_bessel_jn_neg_arg(n, x);
	  return Bess.j_value - s_i * Bess.n_value;
	}
      else
	{
	  const auto Bess = sph_bessel_jn(n, x);
	  return Cmplx{Bess.j_value, -Bess.n_value};
	}
    }

} // namespace detail
} // namespace emsr

#endif // SF_BESSEL_TCC
