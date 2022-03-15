
// Copyright (C) 2016-2019 Free Software Foundation, Inc.
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

/** @file bits/sf_coulomb.tcc
 *  This is an internal header file, included by other library headers.
 *  You should not attempt to use it directly.
 */

#ifndef SF_COULOMB_TCC
#define SF_COULOMB_TCC 1

#include <stdexcept>
#include <complex>
#include <utility>

#include <emsr/numeric_limits.h>
#include <emsr/math_constants.h>
#include <emsr/sf_gamma.h>
#include <emsr/sf_mod_bessel.h> // airy
#include <emsr/notsospecfun.h>

namespace emsr
{
namespace detail
{

/**
 * Return the digamma function
 * @f[
 *   \psi(1+iy)
 * @f[
 * For some reason GSL goes ham on this (gsl_sf_psi_1piy_e)
 * A&S has lots of formulea on this too.
 */
template<typename Tp>
  Tp
  re_digamma_1pi(Tp y)
  {
    // FIXME: Revisit this.
    return std::real(emsr::digamma(std::complex(Tp{1}, y)));
  }

/**
 * 
 */
template<typename Tp>
  Tp
  coulomb_norm(Tp lambda, Tp eta)
  {
    using Cmplx = std::complex<Tp>;
    const auto s_pi = emsr::pi_v<Tp>;
    const auto s_ln2 = emsr::ln2_v<Tp>;

    auto lgam_numer = log_gamma(Cmplx(lambda + Tp{1}, eta));
    auto lgam_denom = log_gamma(Cmplx(Tp{2} * (lambda + Tp{1})));

    auto exp_arg = lambda * s_ln2 - Tp{0.5} * eta * s_pi + lgam_numer - lgam_denom;

    return std::abs(std::exp(exp_arg));
  }

/* Determine the connection phase, phi_lambda.
 * See coulomb_FG_series() below. We have
 * to be careful about sin(phi)->0. Note that
 * there is an underflow condition for large 
 * positive eta in any case.
 */
template<typename Tp>
  sincos_t<Tp>
  coulomb_connection(Tp lam, Tp eta)
  {
    const auto s_pi = emsr::pi_v<Tp>;
    if (eta > -emsr::log_min<Tp>() / Tp{2} * s_pi - Tp{1})
      {
	const Tp cos_phi = Tp{1};
	const Tp sin_phi = Tp{0};
	return {sin_phi, cos_phi};
      }
    else if (eta > -emsr::log_eps<Tp>() / (Tp{4} * s_pi))
      {
	const Tp eps = Tp{2} * std::exp(-Tp{2} * s_pi * eta);
	const Tp tpl = std::tan(s_pi * lam);
	const Tp dth = eps * tpl / (tpl * tpl + Tp{1});
	const Tp cos_phi = -Tp{1} + Tp{0.5} * dth * dth;
	const Tp sin_phi = -dth;
	return {sin_phi, cos_phi};
      }
    else
      {
	Tp X   = std::tanh(s_pi * eta) / std::tan(s_pi * lam);
	Tp phi = -std::atan(X) - (lam + Tp{0.5}) * s_pi;
	return emsr::sincos(phi);
      }
  }

/**
 * A structure for the WKB evaluation of F, G.
 * The F, G functions are scaled such that:
 *
 *   F = f_value * exp(-exponent)
 *   G = g_value * exp(+exponent)
 */
template<typename Tp>
  struct Coulomb
  {
    Tp F_value;
    Tp G_value;
  };

/* Evaluate the Frobenius series for F_lam(eta,rho) and G_lam(eta,rho).
 * Homegrown algebra. Evaluates the series for F_{lam} and
 * F_{-lam-1}, then uses
 *    G_{lam} = (F_{lam} cos(phi) - F_{-lam-1}) / sin(phi)
 * where
 *    phi = Arg[Gamma[1+lam+I eta]] - Arg[Gamma[-lam + I eta]] - (lam+1/2)Pi
 *        = Arg[Sin[Pi(-lam+I eta)] - (lam+1/2)Pi
 *        = atan2(-cos(lam Pi)sinh(eta Pi), -sin(lam Pi)cosh(eta Pi)) - (lam+1/2)Pi
 *
 *        = -atan(X) - (lam+1/2) Pi,  X = tanh(eta Pi)/tan(lam Pi)
 *
 * Not appropriate for lam <= -1/2, lam = 0, or lam >= 1/2.
 */
template<typename Tp>
  Coulomb<Tp>
  coulomb_FG_series(Tp lam, Tp eta, Tp rho)
  {
    const int max_iter = 800;
    const auto s_eps = emsr::epsilon<Tp>();
    auto ClamA = coulomb_norm(lam, eta);
    auto ClamB = coulomb_norm(-lam - Tp{1}, eta);
    const auto tlp1 = Tp{2} * lam + Tp{1};
    const auto pow_x = std::pow(rho, lam);

    auto uA_mm2 = Tp{1};    // uA sum is for F_{lam}
    auto uA_mm1 = rho * eta / (lam + Tp{1});
    auto uB_mm2 = Tp{1};    // uB sum is for F_{-lam-1}
    auto uB_mm1 = -rho * eta / lam;
    auto A_sum = uA_mm2 + uA_mm1;
    auto B_sum = uB_mm2 + uB_mm1;
    auto A_abs_del_prev = std::abs(A_sum);
    auto B_abs_del_prev = std::abs(B_sum);
    int m = 2;

    auto [cos_phi_lam, sin_phi_lam] = coulomb_connection(lam, eta);

    while (m < max_iter)
      {
	auto uA_m = rho * (Tp{2} * eta * uA_mm1 - rho * uA_mm2) / (m * (m + tlp1));
	auto uB_m = rho * (Tp{2} * eta * uB_mm1 - rho * uB_mm2) / (m * (m - tlp1));
	A_sum += uA_m;
	B_sum += uB_m;
	auto abs_dA = std::abs(uA_m);
	auto abs_dB = std::abs(uB_m);
	if (m > 15)
	  {
	    /* Don't bother checking until we have gone out a little ways;
	     * a minor optimization. Also make sure to check both the
	     * current and the previous increment because the odd and even
	     * terms of the sum can have very different behaviour, depending
	     * on the value of eta.
	     */
	    auto max_abs_dA = std::max(abs_dA, A_abs_del_prev);
	    auto max_abs_dB = std::max(abs_dB, B_abs_del_prev);
	    auto abs_A = std::abs(A_sum);
	    auto abs_B = std::abs(B_sum);
	    if (   max_abs_dA / (max_abs_dA + abs_A) < Tp{4} * s_eps
                && max_abs_dB / (max_abs_dB + abs_B) < Tp{4} * s_eps)
              break;
	  }
	A_abs_del_prev = abs_dA;
	B_abs_del_prev = abs_dB;
	uA_mm2 = uA_mm1;
	uA_mm1 = uA_m;
	uB_mm2 = uB_mm1;
	uB_mm1 = uB_m;
	++m;
      }

    auto FA = A_sum * ClamA * pow_x * rho;
    auto FB = B_sum * ClamB / pow_x;

    auto F = FA;
    auto G = (FA * cos_phi_lam - FB) / sin_phi_lam;

    if (m >= max_iter)
      throw std::logic_error("coulomb_FG_series: Exceeded maximum iterations");
    else
      return {F, G};
  }


/* Evaluate the Frobenius series for F_0(eta,rho) and G_0(eta,rho).
 * See [Bardin et al., CPC 3, 73 (1972), (14)-(17)];
 * note the misprint in (17): nu_0=1 is correct, not nu_0=0.
 */
template<typename Tp>
  Coulomb<Tp>
  coulomb_FG0_series(Tp eta, Tp rho)
  {
    const int max_iter = 800;
    const auto s_eps = emsr::epsilon<Tp>();
    const auto s_gamma_e = emsr::egamma_v<Tp>;
    const auto rho2  = rho * rho;
    const auto tex = Tp{2} * eta * rho;
    auto C0 = coulomb_norm(Tp{0}, eta);
    auto r1pie = re_digamma_1pi(eta);
    auto u_mm2 = Tp{0};   /* u_0 */
    auto u_mm1 = rho;       /* u_1 */
    auto v_mm2 = Tp{1};   /* nu_0 */
    auto v_mm1 = tex * (Tp{2} * s_gamma_e - Tp{1} + r1pie);  /* nu_1 */
    auto u_sum = u_mm2 + u_mm1;
    auto v_sum = v_mm2 + v_mm1;
    auto u_abs_del_prev = std::abs(u_sum);
    auto v_abs_del_prev = std::abs(v_sum);
    int m = 2;
    const auto ln2x = std::log(Tp{2} * rho);

    while (m < max_iter)
      {
	Tp m_mm1 = m * (m - Tp{1});
	auto u_m = (tex * u_mm1 - rho2 * u_mm2) / m_mm1;
	auto v_m = (tex * v_mm1 - rho2 * v_mm2 - Tp{2} * eta * (2 * m - 1) * u_m) / m_mm1;
	u_sum += u_m;
	v_sum += v_m;
	auto abs_du = std::abs(u_m);
	auto abs_dv = std::abs(v_m);
	if (m > 15)
	  {
	    /* Don't bother checking until we have gone out a little ways;
	     * a minor optimization. Also make sure to check both the
	     * current and the previous increment because the odd and even
	     * terms of the sum can have very different behaviour, depending
	     * on the value of eta.
	     */
	    auto max_abs_du = std::max(abs_du, u_abs_del_prev);
	    auto max_abs_dv = std::max(abs_dv, v_abs_del_prev);
	    auto abs_u = std::abs(u_sum);
	    auto abs_v = std::abs(v_sum);
	    if (   max_abs_du / (max_abs_du + abs_u) < Tp{40} * s_eps
        	&& max_abs_dv / (max_abs_dv + abs_v) < Tp{40} * s_eps)
              break;
	  }
	u_abs_del_prev = abs_du;
	v_abs_del_prev = abs_dv;
	u_mm2 = u_mm1;
	u_mm1 = u_m;
	v_mm2 = v_mm1;
	v_mm1 = v_m;
	++m;
      }

    auto F  = C0 * u_sum;
    auto G  = (v_sum + Tp{2} * eta * u_sum * ln2x) / C0;

    if (m == max_iter)
      throw std::logic_error("coulomb_FG0_series: Exceeded maximum iterations");
    else
      return {F, G};
  }


/* Evaluate the Frobenius series for F_{-1/2}(eta,rho) and G_{-1/2}(eta,rho).
 * Homegrown algebra.
 */
template<typename Tp>
  Coulomb<Tp>
  coulomb_FGmhalf_series(Tp eta, Tp rho)
  {
    const int max_iter = 800;
    const auto s_eps = emsr::epsilon<Tp>();
    const auto s_gamma_e = emsr::egamma_v<Tp>;
    const auto s_ln2 = emsr::ln2_v<Tp>;
    const Tp rx  = std::sqrt(rho);
    const Tp rho2  = rho * rho;
    const Tp tex = Tp{2} * eta * rho;
    auto Cmhalf = coulomb_norm(-0.5, eta);
    Tp u_mm2 = Tp{1};         /* u_0 */
    Tp u_mm1 = tex * u_mm2; /* u_1 */
    int m = 2;

    // Re[\psi(1 + i \eta)]
    auto rpsi_1pe = re_digamma_1pi(eta);
    // Re[\psi(1 + i 2\eta)]
    auto rpsi_1p2e = re_digamma_1pi(Tp{2} * eta);

    Tp v_mm2 = Tp{2} * s_gamma_e - s_ln2 - rpsi_1pe + Tp{2} * rpsi_1p2e;
    Tp v_mm1 = tex * (v_mm2 - Tp{2} * u_mm2);

    Tp f_sum = u_mm2 + u_mm1;
    Tp g_sum = v_mm2 + v_mm1;

    while (m < max_iter)
      {
	Tp m2 = m * m;
	auto u_m = (tex * u_mm1 - rho2 * u_mm2) / m2;
	auto v_m = (tex * v_mm1 - rho2 * v_mm2 - Tp{2} * m * u_m) / m2;
	f_sum += u_m;
	g_sum += v_m;
	if (f_sum != Tp{0} && g_sum != Tp{0}
	   && (std::abs(u_m / f_sum) + std::abs(v_m / g_sum) < Tp{10} * s_eps))
          break;
	u_mm2 = u_mm1;
	u_mm1 = u_m;
	v_mm2 = v_mm1;
	v_mm1 = v_m;
	++m;
      }

    auto F = Cmhalf * rx * f_sum;
    auto G = -rx * (f_sum * std::log(rho) + g_sum) / Cmhalf;

    if (m == max_iter)
      throw std::logic_error("coulomb_FGmhalf_series: Exceeded maximum iterations");
    else
      return {F, G};
  }

  /**
   * Evolve the backwards recurrence for F, F'.
   * @f[
   *    F_{l-1}  = (S_l F_l + F_l') / R_l
   *    F_{l-1}' = (S_l F_{l-1} - R_l F_l)
   * @f]
   * where
   * @f[
   *    R_l = \sqrt{1 + (\eta / l)^2}
   *    S_l = l / \rho + \eta / l
   * @f]
   */
  template<typename Tp>
    std::pair<Tp, Tp>
    coulomb_f_recur(Tp lam_min, unsigned int k_max,
		    Tp eta, Tp rho, Tp F_lam_max, Tp Fp_lam_max)
    {
      const auto x_inv = Tp{1} / rho;
      auto fcl = F_lam_max;
      auto fpl = Fp_lam_max;
      const auto lam_max = lam_min + Tp(k_max);
      auto lam = lam_max;

      for (int k = k_max - 1; k >= 0; --k)
	{
	  const auto el = eta / lam;
	  const auto rl = std::hypot(Tp{1}, el);
	  const auto sl = el + lam * x_inv;
	  const auto fc_lm1 = (fcl * sl + fpl) / rl;
	  fpl = fc_lm1 * sl - fcl * rl;
	  fcl = fc_lm1;
	  lam -= Tp{1};
	}

      return std::make_pair(fcl, fpl);
    }

  /**
   * Evolve the forward recurrence for G, G'.
   * @f[
   *   G_{l+1}  = (S_l G_l - G_l')/R_l
   *   G_{l+1}' = R_{l+1} G_l - S_l G_{l+1}
   * @f]
   * where
   * @f[
   *    R_l = \sqrt{1 + (\eta / l)^2}
   *    S_l = l / \rho + \eta / l
   * @f]
   */
  template<typename Tp>
    std::pair<Tp, Tp>
    coulomb_g_recur(Tp lam_min, unsigned int k_max,
		    Tp eta, Tp rho, Tp G_lam_min, Tp Gp_lam_min)
    {
      const auto x_inv = Tp{1} / rho;
      auto gcl = G_lam_min;
      auto gpl = Gp_lam_min;
      auto lam = lam_min + Tp{1};

      for (unsigned int k = 1; k <= k_max; ++k)
	{
	  const auto el = eta / lam;
	  const auto rl = std::hypot(Tp{1}, el);
	  const auto sl = el + lam * x_inv;
	  const auto gc_lm1 = (sl * gcl - gpl) / rl;
	  gpl = rl * gcl - sl * gc_lm1;
	  gcl = gc_lm1;
	  lam += Tp{1};
	}

      return std::make_pair(gcl, gpl);
    }

  /**
   * Evaluate the first continued fraction, giving the ratio F'/F
   * at the upper l value.
   * We also determine the sign of F at that point, since it is the sign
   * of the last denominator in the continued fraction.
   */
  template<typename Tp>
    std::pair<Tp, Tp>
    coulomb_CF1(Tp lambda, Tp eta, Tp rho)
    {
      const auto CF1_small = Tp{1.0e-30};
      const auto CF1_abort = Tp{1.0e+05};
      const auto CF1_acc = Tp{2} * std::numeric_limits<Tp>::epsilon();
      const auto x_inv = Tp{1} / rho;
      const auto px = lambda + Tp{1} + CF1_abort;

      auto pk = lambda + Tp{1};
      auto F = eta / pk + pk * x_inv;

      auto fcl_sign = Tp{1};

      if (std::abs(F) < CF1_small)
	F = CF1_small;
      auto D = Tp{0};
      auto C = F;

      Tp df;
      do
	{
	  const auto pk1 = pk + Tp{1};
	  const auto ek  = eta / pk;
	  const auto rk2 = Tp{1} + ek * ek;
	  const auto tk  = (pk + pk1) * (x_inv + ek / pk1);
	  D = tk - rk2 * D;
	  C = tk - rk2 / C;
	  if (std::abs(C) < CF1_small)
	    C = CF1_small;
	  if (std::abs(D) < CF1_small)
	    D = CF1_small;
	  D = Tp{1} / D;
	  df = D * C;
	  F *= df;
	  if (D < Tp{0})
	    {
	      // sign of result depends on sign of denominator.
	      fcl_sign = -fcl_sign;
	    }
	  pk = pk1;
	  if (pk > px)
	    throw std::runtime_error("coulomb_CF1: Too many iterations.");
	}
      while (std::abs(df - Tp{1}) > CF1_acc);

      return std::make_pair(fcl_sign, F);
    }

  /**
   * Evaluate the second continued fraction to obtain the ratio
   * @f[
   *    (G' + i F') / (G + i F) := P + i Q
   * @f]
   * at the specified l value.
   */
  template<typename Tp>
    std::complex<Tp>
    coulomb_CF2(Tp lambda, Tp eta, Tp rho)
    {
      const auto s_i = std::complex<Tp>{0, 1};
      const auto CF2_acc = Tp{4} * std::numeric_limits<Tp>::epsilon();
      const auto CF2_abort = Tp{2.0e+05};

      const auto wi = Tp{2} * eta;
      const auto x_inv = Tp{1} / rho;
      const auto e2mm1 = eta * eta + Tp(lambda * (lambda + 1));

      auto a = std::complex<Tp>(-e2mm1, eta);
      auto b = Tp{2} * std::complex<Tp>(rho - eta, Tp{1});
      auto d = std::conj(b) / std::norm(b);

      auto dpq = s_i * x_inv * a * d;

      auto pk = Tp{0};
      auto PQ = std::complex<Tp>(Tp{0}, Tp{1} - eta * x_inv);

      do
        {
	  PQ += dpq;
	  pk += Tp{2};
	  a += std::complex<Tp>(pk, wi);
	  b += Tp{2} * s_i;
	  d = b + a * d;
	  d = std::conj(d) / std::norm(d);
	  dpq *= b * d - Tp{1};
	  if (pk > CF2_abort)
	    throw std::runtime_error("coulomb_CF2: Too many iterations.");
        }
      while (std::abs(dpq) > (std::abs(PQ)) * CF2_acc);

      //if (Q < CF2_abort * std::numeric_limits<Tp>::epsilon() * std::abs(P))
	//status = GSL_ELOSS;

      return PQ;
    }

  /**
   * A structure for the WKB evaluation of F, G.
   * The F, G functions are scaled such that:
   *
   *   F = f_value * exp(-exponent)
   *   G = g_value * exp(+exponent)
   */
  template<typename Tp>
    struct coulomb_jwkb_t
    {
      Tp F_value;
      Tp G_value;
      Tp exponent;
    };

  /**
   * WKB evaluation of F, G. Assumes  0 < rho < turning point.
   * Overflows are trapped, GSL_EOVRFLW is signalled,
   * and an exponent is returned such that:
   *
   *   result_F = fjwkb * exp(-exponent)
   *   result_G = gjwkb * exp( exponent)
   *
   * See [Biedenharn et al. Phys. Rev. 97, 542-554 (1955), Section IV]
   *
   * Unfortunately, this is not very accurate in general. The
   * test cases typically have 3-4 digits of precision. One could
   * argue that this is ok for general use because, for instance,
   * F is exponentially small in this region and so the absolute
   * accuracy is still roughly acceptable. But it would be better
   * to have a systematic method for improving the precision. See
   * the Abad & Sesma method discussion below.
   */
  template<typename Tp>
    coulomb_jwkb_t<Tp>
    coulomb_jwkb(Tp lambda, Tp eta, Tp rho)
    {
      const auto s_pi = emsr::pi_v<Tp>;
      const auto llp1 = lambda * (lambda + Tp{1}) + Tp{6} / Tp{35};
      const auto llp1_eff = std::max(llp1, Tp{0});
      const auto rho_ghalf = std::sqrt(rho * (Tp{2} * eta - rho) + llp1_eff);
      const auto sinh_arg = std::sqrt(llp1_eff / (eta * eta + llp1_eff)) * rho_ghalf / rho;
      const auto sinh_inv = std::log(sinh_arg + std::hypot(Tp{1}, sinh_arg));

      const auto phi = std::abs(rho_ghalf
                              - eta * std::atan2(rho_ghalf, rho - eta)
                              - std::sqrt(llp1_eff) * sinh_inv);

      const auto zeta_half = std::pow(Tp{3} * phi / Tp{2}, Tp{1} / Tp{3});
      const auto prefactor = std::sqrt(s_pi * phi * rho / (Tp{6} * rho_ghalf));

      auto F = prefactor * Tp{3} / zeta_half;
      auto G = prefactor * Tp{3} / zeta_half; // Note the sqrt(3) from Bi normalization.

      const auto airy_scale_exp = phi;
      const auto AiBi = emsr::detail::airy(zeta_half * zeta_half, true); // scaled
      F *= AiBi.Ai_value;
      G *= AiBi.Bi_value;
      auto F_exp = std::log(F) - airy_scale_exp;
      auto G_exp = std::log(G) + airy_scale_exp;

      if (G_exp >= emsr::log_max<Tp>())
	{
	  return {F, G, airy_scale_exp};
	}
      else
	{
	  return {std::exp(F_exp), std::exp(G_exp), Tp{0}};
	}
    }

/**
 * This is the main deal. It computes the two Coulomg functions ad their derivatives.
 *
 * @param  lam_F  The lambda or degree parameter for the F function.
 * @param  eta  The Sommerfeld parameter.
 * @param  rho  The radial distance.
 * @param  k_lam_G  The integral difference between the F lambda and the G lambda: lam_G = lam_F - k_lam_G
 */
template<typename Tp>
  emsr::Coulomb<Tp>
  coulomb_wave_FG(Tp lam_F, Tp eta, Tp rho, int k_lam_G)
  {
    const auto s_pi = emsr::pi_v<Tp>;
    const auto s_inf = emsr::infinity<Tp>();
    const Tp lam_G = lam_F - k_lam_G;

    Tp F, Fp, G, Gp, exp_F, exp_G;

    if (rho < Tp{0} || lam_F <= -Tp{0.5} || lam_G <= -Tp{0.5})
      {
	throw std::domain_error("coulomb_wave_FG: domain error");
      }
    else if (rho == Tp{0})
      {
	const auto C0 = coulomb_norm(Tp{0}, eta);
	F = Tp{0};
	Fp = Tp{0};
	G = s_inf;
	Gp = s_inf;
	exp_F = Tp{0};
	exp_G = Tp{0};
	if (lam_F == Tp{0})
	  Fp = C0;
	if (lam_G == Tp{0})
	  Gp = Tp{1} / C0;
	return {F, Fp, exp_F, G, Gp, exp_G};
      }
    else if (rho < 1.2 && Tp{2} * s_pi * eta < 0.9 * (-emsr::log_min<Tp>()) && std::abs(eta*rho) < Tp{10})
      {
	/* Reduce to a small lambda value and use the series
	 * representations for F and G. We cannot allow eta to
	 * be large and positive because the connection formula
	 * for G_lam is badly behaved due to an underflow in sin(phi_lam) 
	 * [see coulomb_FG_series() and coulomb_connection() above].
	 * Note that large negative eta is ok however.
	 */
	const Tp SMALL = emsr::sqrt_eps<Tp>();
	const int N    = static_cast<int>(lam_F + Tp{0.5});
	const int span = std::max(k_lam_G, N);
	const Tp lam_min = lam_F - N;    /* -1/2 <= lam_min < 1/2 */
	Tp G_lam_G = Tp{0}, Gp_lam_G = Tp{0};
	Tp Fp_over_F_lam_min;
	Tp F_lam_min, G_lam_min, Gp_lam_min;
	Tp F_scale;

	// Determine F'/F at lam_F.
	auto F_CF1 = coulomb_CF1(lam_F, eta, rho);
        auto Fp_over_F_lam_F = F_CF1.second;

	// Recurse down with unnormalized F, F' values.
	auto F_lam_F  = SMALL;
	auto Fp_lam_F = Fp_over_F_lam_F * F_lam_F;
	Tp F_lam_min_unnorm, Fp_lam_min_unnorm;
	if (span != 0)
          {
	    auto Fr = coulomb_f_recur(lam_min, span, eta, rho, F_lam_F, Fp_lam_F);
            F_lam_min_unnorm = Fr.first;
            Fp_lam_min_unnorm = Fr.second;
	  }
	else
          {
	    F_lam_min_unnorm = F_lam_F;
	    Fp_lam_min_unnorm = Fp_lam_F;
	  }

	// Determine F and G at lam_min.
	if (lam_min == -0.5)
          {
	    auto FG_ser = coulomb_FGmhalf_series(eta, rho);
            F_lam_min = FG_ser.F_value;
            G_lam_min = FG_ser.G_value;
	  }
	else if (lam_min == Tp{0})
          {
	    auto FG_ser = coulomb_FG0_series(eta, rho);
            F_lam_min = FG_ser.F_value;
            G_lam_min = FG_ser.G_value;
	  }
	else if (lam_min == 0.5)
          {
	    // This cannot happen.
	    F = F_lam_F;
	    Fp = Fp_lam_F;
	    exp_F = Tp{0};
	    G = G_lam_G;
	    Gp = Gp_lam_G;
	    exp_G = Tp{0};
	    //GSL_ERROR ("error", GSL_ESANITY);
            return {F, Fp, exp_F, G, Gp, exp_G};
	  }
	else
          {
	    auto FG_ser = coulomb_FG_series(lam_min, eta, rho);
            F_lam_min = FG_ser.F_value;
            G_lam_min = FG_ser.G_value;
	  }

	// Determine remaining quantities.
	Fp_over_F_lam_min = Fp_lam_min_unnorm / F_lam_min_unnorm;
	Gp_lam_min  = Fp_over_F_lam_min * G_lam_min - Tp{1} / F_lam_min;
	F_scale = F_lam_min / F_lam_min_unnorm;

	// Apply scale to the original F,F' values.
	F_lam_F *= F_scale;
	Fp_lam_F *= F_scale;

	// Recurse up to get the required G, G' values.
	auto Gr = coulomb_g_recur(lam_min, std::max(N - k_lam_G, 0), eta, rho, G_lam_min, Gp_lam_min);
        G_lam_min = Gr.first;
        Gp_lam_min = Gr.second;

	F = F_lam_F;
	Fp = Fp_lam_F;
	G = Gr.first;
	Gp = Gr.second;
	exp_F = Tp{0};
	exp_G = Tp{0};

	return {F, Fp, exp_F, G, Gp, exp_G};
      }
    else if (rho < Tp{2} * eta)
      {
	/* Use WKB approximation to obtain F and G at the two
	 * lambda values, and use the Wronskian and the
	 * continued fractions for F'/F to obtain F' and G'.
	 */

	auto jwkb = coulomb_jwkb(lam_F, eta, rho);
        auto F_lam_F = jwkb.F_value;
        auto exp_lam_F = jwkb.exponent;
        Tp F_lam_G, G_lam_G, exp_lam_G;
	if (k_lam_G == 0)
          {
	    F_lam_G = jwkb.F_value;
	    G_lam_G = jwkb.G_value;
	    exp_lam_G = exp_lam_F;
	  }
	else
          {
	    auto jwkb = coulomb_jwkb(lam_G, eta, rho);
            F_lam_G = jwkb.F_value;
	    G_lam_G = jwkb.G_value;
            exp_lam_G = jwkb.exponent;
          }

	auto F_CF1 = coulomb_CF1(lam_F, eta, rho);
        auto Fp_over_F_lam_F = F_CF1.second;
	Tp Fp_over_F_lam_G;
	if (k_lam_G == 0)
          {
	    Fp_over_F_lam_G = Fp_over_F_lam_F;
	  }
	else
          {
	    auto G_CF1 = coulomb_CF1(lam_G, eta, rho);
	    Fp_over_F_lam_G = G_CF1.second;
          }
	F = F_lam_F;
	G = G_lam_G;
	Fp = Fp_over_F_lam_F * F_lam_F;
	Gp = Fp_over_F_lam_G * G_lam_G - Tp{1} / F_lam_G;
	exp_F = exp_lam_F;
	exp_G = exp_lam_G;

	return  {F, Fp, exp_F, G, Gp, exp_G};
      }
    else
      {
	/* rho > 2 eta, so we know that we can find a lambda value such
	 * that rho is above the turning point. We do this, evaluate
	 * using Steed's method at that oscillatory point, then
	 * use recursion on F and G to obtain the required values.
	 *
	 * lam_0   = a value of lambda such that rho is below the turning point
	 * lam_min = minimum of lam_0 and the requested lam_G, since
	 *           we must go at least as low as lam_G
	 */
	const Tp SMALL = emsr::sqrt_eps<Tp>();
	const Tp C = std::sqrt(Tp{1} + Tp{4} * rho * (rho - Tp{2} * eta));
	const int N = std::ceil(lam_F - C + 0.5);
	const Tp lam_0 = lam_F - std::max(N, 0);
	const Tp lam_min = std::min(lam_0, lam_G);

	auto F_CF1 = coulomb_CF1(lam_F, eta, rho);
	auto F_sign_lam_F = F_CF1.first;
	auto Fp_over_F_lam_F = F_CF1.second;

	auto F_lam_F  = F_sign_lam_F * SMALL; // unnormalized
	auto Fp_lam_F = Fp_over_F_lam_F * F_lam_F;

	// Backward recurrence to get F,Fp at lam_min
        int F_recur_count = std::max(k_lam_G, N);
	auto Fr = coulomb_f_recur(lam_min, F_recur_count, eta, rho, F_lam_F, Fp_lam_F);
        auto F_lam_min_unnorm = Fr.first;
	auto Fp_over_F_lam_min = Fr.second / Fr.first;

	// Steed evaluation to complete evaluation of F, Fp, G, Gp at lam_min
	auto CF2 = coulomb_CF2(lam_min, eta, rho);
	auto P_lam_min = std::real(CF2);
	auto Q_lam_min = std::imag(CF2);
	auto alpha = Fp_over_F_lam_min - P_lam_min;
	auto gamma = alpha / Q_lam_min;

	auto F_sign_lam_min = emsr::sign(F_lam_min_unnorm);

	auto F_lam_min  = F_sign_lam_min / std::sqrt(alpha * alpha / Q_lam_min + Q_lam_min);
	auto G_lam_min  = gamma * F_lam_min;
	auto Gp_lam_min = (P_lam_min * gamma - Q_lam_min) * F_lam_min;

	// Apply scale to values of F,Fp at lam_F (the top).
	auto F_scale = F_lam_min / F_lam_min_unnorm;    
	F_lam_F *= F_scale;
	Fp_lam_F *= F_scale;

	// Forward recurrence to get G,Gp at lam_G (the top).
        int G_recur_count = std::max(N - k_lam_G, 0);
	auto Gr = coulomb_g_recur(lam_min, G_recur_count, eta, rho, G_lam_min, Gp_lam_min);

	F  = F_lam_F;
	Fp = Fp_lam_F;
	G  = Gr.first;
	Gp = Gr.second;
	exp_F = Tp{0};
	exp_G = Tp{0};

	return  {F, Fp, exp_F, G, Gp, exp_G};
      }
  }

  /**
   * Return the bound-state Coulomb wave-function.
   */
  template <typename Tp>
    std::complex<Tp>
    hydrogen(unsigned int n,
             unsigned int l, unsigned int m,
             Tp Z, Tp r, Tp theta, Tp phi)
    {
      const auto s_NaN = emsr::quiet_NaN(r);

      if (std::isnan(Z) || std::isnan(r)
	 || std::isnan(theta) || std::isnan(phi))
	return std::complex<Tp>{s_NaN, s_NaN};
      else if(n < 1)
	throw std::domain_error("hydrogen: level number less than one");
      else if(l > n - 1)
	throw std::domain_error("hydrogen: angular momentum number too large");
      else if(Z <= Tp(0))
	throw std::domain_error("hydrogen: non-positive charge");
      else if(r < Tp(0))
	throw std::domain_error("hydrogen: negative radius");
      else
	{
	  const auto A = Tp(2) * Z / n;

	  const auto pre = std::sqrt(A * A * A / (Tp(2) * n));
	  const auto ln_a = emsr::lgamma(n + l + 1);
	  const auto ln_b = emsr::lgamma(n - l);
	  const auto ex = std::exp((ln_b - ln_a) / Tp(2));
	  const auto norm = pre * ex;

	  const auto rho = A * r;
	  const auto ea = std::exp(-rho / Tp(2));
	  const auto pp = std::pow(rho, l);
	  const auto lag = assoc_laguerre(n - l - 1, 2 * l + 1, rho);
	  const auto sphh = sph_legendre(l, m, theta) * std::polar(Tp(1), Tp(m) * phi);

	  const auto psi = norm * ea * pp * lag * sphh;

	  return psi;
	}
    }

} // namespace detail
} // namespace emsr

#endif // SF_COULOMB_TCC
