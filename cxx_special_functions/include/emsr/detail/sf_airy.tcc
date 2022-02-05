
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

/** @file bits/sf_airy.tcc
 * This is an internal header file, included by other library headers.
 * You should not attempt to use it directly.
 */

#ifndef SF_AIRY_TCC
#define SF_AIRY_TCC 1

#include <complex>

#include <emsr/summation.h>
#include <emsr/polynomial.h>
#include <emsr/math_constants.h>

namespace emsr
{
namespace detail
{

  /**
   * This struct defines the Airy function state with two presumably
   * numerically useful Airy functions and their derivatives.
   * The data mambers are directly accessible.
   * The lone method computes the Wronskian from the stored functions.
   * A static method returns the correct Wronskian.
   */
  template<typename Tp>
    struct AiryState
    {
      using Real = emsr::num_traits_t<Tp>;

      Tp z;
      Tp Ai_value;
      Tp Ai_deriv;
      Tp Bi_value;
      Tp Bi_deriv;

      Tp
      Wronskian() const
      { return Ai_value * Bi_deriv - Bi_value * Ai_deriv; }

      Real
      true_Wronskian()
      { return Real{1} / emsr::pi_v<Real>; }
    };


  /**
   * A structure containing three auxilliary Airy functions
   * and their derivatives.
   */
  template<typename Tp>
    struct AiryAuxilliaryState
    {
      using Val = emsr::num_traits_t<Tp>;

      Tp z;
      Tp fai_value;
      Tp fai_deriv;
      Tp gai_value;
      Tp gai_deriv;
      Tp hai_value;
      Tp hai_deriv;
    };


  /**
   * This class orgianizes series solutions of the Airy function.
   * @f[
   *    fai(x) = \sum_{k=0}^\infty \frac{(2k+1)!!!x^{3k}}{(2k+1)!}
   * @f]
   * @f[
   *    gai(x) = \sum_{k=0}^\infty \frac{(2k+2)!!!x^{3k+1}}{(2k+2)!}
   * @f]
   * @f[
   *    hai(x) = \sum_{k=0}^\infty \frac{(2k+3)!!!x^{3k+2}}{(2k+3)!}
   * @f]
   * This class contains tabulations of the factors appearing in the sums above.
   */
  template<typename Tp>
    class AirySeries
    {
    public: // FIXME!!!

      using Cmplx = std::complex<Tp>;

      static constexpr int N_FGH = 200;
    private: // FIXME!!!
      static constexpr Tp s_slope_F{-2.660L}, s_intercept_F{-0.778L};
      static constexpr Tp s_slope_Fp{-2.576L}, s_intercept_Fp{-0.301L};
      static constexpr Tp s_slope_G{-2.708L}, s_intercept_G{-1.079L};
      static constexpr Tp s_slope_Gp{-2.632L}, s_intercept_Gp{-0.477L};
      static constexpr Tp s_slope_H{-2.75L}, s_intercept_H{-1.25L};
      static constexpr Tp s_slope_Hp{-2.625L}, s_intercept_Hp{-0.6L};

    public:

      static constexpr Tp s_eps = emsr::epsilon(Tp{});
      static constexpr Tp s_pi = emsr::pi_v<Tp>;
      static constexpr Tp s_sqrt_pi
      		 = emsr::sqrtpi_v<Tp>;
      static constexpr Tp s_Ai0
      		 = Tp{3.550280538878172392600631860041831763980e-1L};
      static constexpr Tp s_Aip0
      		 = Tp{-2.588194037928067984051835601892039634793e-1L};
      static constexpr Tp s_Bi0
      		 = Tp{6.149266274460007351509223690936135535960e-1L};
      static constexpr Tp s_Bip0
      		 = Tp{4.482883573538263579148237103988283908668e-1L};
      static constexpr Tp s_Hi0
      		 = Tp{4.099510849640004901006149127290757023959e-1L};
      static constexpr Tp s_Hip0
      		 = Tp{2.988589049025509052765491402658855939102e-1L};
      static constexpr Tp s_Gi0
      		 = Tp{2.049755424820002450503074563645378511979e-1L};
      static constexpr Tp s_Gip0
      		 = Tp{1.494294524512754526382745701329427969551e-1L};
      static constexpr Cmplx s_i{Tp{0}, Tp{1}};

      static AiryState<Cmplx>
      s_Airy(Cmplx t);

      static AiryState<Cmplx>
      s_Fock(Cmplx t);

      static std::pair<Cmplx, Cmplx>
      s_Ai(Cmplx t);

      static std::pair<Cmplx, Cmplx>
      s_Bi(Cmplx t);

      static AiryAuxilliaryState<Cmplx>
      s_FGH(Cmplx t);

      static AiryState<Cmplx>
      s_Scorer(Cmplx t);
      static AiryState<Cmplx>
      s_Scorer2(Cmplx t);

    private:

      static AiryState<Cmplx>
      s_AiryHelp(Cmplx t, bool return_fock_airy = false);

      std::pair<Cmplx, Cmplx>
      static s_AiBi(Cmplx t, std::pair<Tp, Tp> Z0);
    };

  // Type-dependent limits for the arrays.
  // FIXME: Make these limits digits10-based.
  template<typename Tp>
    constexpr int max_FGH = AirySeries<Tp>::N_FGH;

  template<>
    constexpr int max_FGH<float> = 15;

  template<>
    constexpr int max_FGH<double> = 79;

  template<typename Tp>
    constexpr Tp
    AirySeries<Tp>::s_eps;

  template<typename Tp>
    constexpr Tp
    AirySeries<Tp>::s_pi;

  template<typename Tp>
    constexpr Tp
    AirySeries<Tp>::s_sqrt_pi;

  template<typename Tp>
    constexpr Tp
    AirySeries<Tp>::s_Ai0;

  template<typename Tp>
    constexpr Tp
    AirySeries<Tp>::s_Aip0;

  template<typename Tp>
    constexpr Tp
    AirySeries<Tp>::s_Bi0;

  template<typename Tp>
    constexpr Tp
    AirySeries<Tp>::s_Bip0;

  template<typename Tp>
    constexpr Tp
    AirySeries<Tp>::s_Hi0;

  template<typename Tp>
    constexpr Tp
    AirySeries<Tp>::s_Hip0;

  template<typename Tp>
    constexpr Tp
    AirySeries<Tp>::s_Gi0;

  template<typename Tp>
    constexpr Tp
    AirySeries<Tp>::s_Gip0;

  template<typename Tp>
    constexpr std::complex<Tp>
    AirySeries<Tp>::s_i;

  /**
   * Return the Airy functions by using the series expansions of
   * the auxilliary Airy functions:
   * @f[
   *    fai(x) = \sum_{k=0}^\infty \frac{(2k+1)!!!x^{3k}}{(2k+1)!}
   * @f]
   * @f[
   *    gai(x) = \sum_{k=0}^\infty \frac{(2k+2)!!!x^{3k+1}}{(2k+2)!}
   * @f]
   * The Airy functions are then defined by:
   * @f[
   *    Ai(x) = Ai(0)fai(x) + Ai'(0)gai(x)
   * @f]
   * @f[
   *    Bi(x) = Bi(0)fai(x) + Bi'(0)gai(x)
   * @f]
   * where @f$ Ai(0) = 3^{-2/3}/\Gamma(2/3) @f$, @f$ Ai'(0) = -3{-1/2}Bi'(0) @f$
   * and @f$ Bi(0) = 3^{1/2}Ai(0) @f$, @f$ Bi'(0) = 3^{1/6}/\Gamma(1/3) @f$
   *
   * @tparam Tp A real type
   */
  template<typename Tp>
    AiryState<std::complex<Tp>>
    AirySeries<Tp>::s_AiryHelp(std::complex<Tp> t,
				   bool return_fock_airy)
    {
      const Tp s_log10min = emsr::log10_min(Tp{});
      const auto s_min = std::numeric_limits<Tp>::min();
      const auto log10t = std::log10(std::abs(t));
      const auto ttt = t * t * t;

      auto Fai = Tp{1};
      auto Gai = Tp{1};
      auto Faip = Tp{0};
      auto Gaip = Tp{1};
      auto term = Cmplx{Tp{1}};
      auto F = Cmplx{Tp{1}};
      auto G = t;
      auto Fp = Cmplx{Tp{0}};
      auto Gp = Cmplx{Tp{1}};
      for (int k = 1; k < max_FGH<Tp>; ++k)
	{
	  if (std::abs(t) < s_eps)
	    break;

	  auto xx = log10t * (3 * k + 1)
		    + s_slope_G * k + s_intercept_G;
	  if (xx < s_log10min)
	    break;

	  Fai /= (3ULL * k - 1ULL) * (3ULL * k);
	  if (Fai < Tp{10} * s_min)
	    break;
	  Faip = (3ULL * k) * Fai;
	  Gai /= (3ULL * k) * (3ULL * k + 1ULL);
	  Gaip = (3ULL * k + 1ULL) * Gai;

	  term *= ttt;
	  F += Fai * term;
	  G += Gai * term * t;
	  Fp += Faip * term / t;
	  Gp += Gaip * term;
	}
      auto UU = s_sqrt_pi * (s_Bi0 * F + s_Bip0 * G);
      auto VV = s_sqrt_pi * (s_Ai0 * F + s_Aip0 * G);
      auto UUp = s_sqrt_pi * (s_Bi0 * Fp + s_Bip0 * Gp);
      auto VVp = s_sqrt_pi * (s_Ai0 * Fp + s_Aip0 * Gp);

      if (!return_fock_airy)
	{
	  auto Bi = UU / s_sqrt_pi;
	  auto Ai = VV / s_sqrt_pi;
	  auto Bip = UUp / s_sqrt_pi;
	  auto Aip = VVp / s_sqrt_pi;
	  return AiryState<std::complex<Tp>>{t, Ai, Aip, Bi, Bip};
	}
      else
	{
	  auto w1 = UU - s_i * VV;
	  auto w2 = UU + s_i * VV;
	  auto w1p = UUp - s_i * VVp;
	  auto w2p = UUp + s_i * VVp;
	  return AiryState<std::complex<Tp>>{t, w1, w1p, w2, w2p};
	}
    }

  /**
   * Return the Airy function of the first kind and its derivative
   * by using the series expansions of the auxilliary Airy functions:
   * @f[
   *    fai(x) = \sum_{k=0}^\infty \frac{(2k+1)!!!x^{3k}}{(2k+1)!}
   * @f]
   * @f[
   *    gai(x) = \sum_{k=0}^\infty \frac{(2k+2)!!!x^{3k+1}}{(2k+2)!}
   * @f]
   * The Airy function of the first kind is then defined by:
   * @f[
   *    Ai(x) = Ai(0)fai(x) + Ai'(0)gai(x)
   * @f]
   * where @f$ Ai(0) = 3^{-2/3}/\Gamma(2/3) @f$, @f$ Ai'(0) = -3{-1/2}Bi'(0) @f$
   * and @f$ Bi(0) = 3^{1/2}Ai(0) @f$, @f$ Bi'(0) = 3^{1/6}/\Gamma(1/3) @f$
   *
   * @tparam Tp A real type
   */
  template<typename Tp>
    std::pair<std::complex<Tp>, std::complex<Tp>>
    AirySeries<Tp>::s_Ai(std::complex<Tp> t)
    { return s_AiBi(t, std::make_pair(s_Ai0, s_Aip0)); }

  /**
   * Return the Airy function of the second kind and its derivative
   * by using the series expansions of the auxilliary Airy functions:
   * @f[
   *    fai(x) = \sum_{k=0}^\infty \frac{(2k+1)!!!x^{3k}}{(2k+1)!}
   * @f]
   * @f[
   *    gai(x) = \sum_{k=0}^\infty \frac{(2k+2)!!!x^{3k+1}}{(2k+2)!}
   * @f]
   * The Airy function of the second kind is then defined by:
   * @f[
   *    Bi(x) = Bi(0)fai(x) + Bi'(0)gai(x)
   * @f]
   * where @f$ Ai(0) = 3^{-2/3}/\Gamma(2/3) @f$, @f$ Ai'(0) = -3{-1/2}Bi'(0) @f$
   * and @f$ Bi(0) = 3^{1/2}Ai(0) @f$, @f$ Bi'(0) = 3^{1/6}/\Gamma(1/3) @f$
   *
   * @tparam Tp A real type
   */
  template<typename Tp>
    std::pair<std::complex<Tp>, std::complex<Tp>>
    AirySeries<Tp>::s_Bi(std::complex<Tp> t)
    { return s_AiBi(t, std::make_pair(s_Bi0, s_Bip0)); }

  /**
   * Return the auxilliary Airy functions:
   * @f[
   *    fai(x) = \sum_{k=0}^\infty \frac{(2k+1)!!!x^{3k}}{(2k+1)!}
   * @f]
   * @f[
   *    gai(x) = \sum_{k=0}^\infty \frac{(2k+2)!!!x^{3k+1}}{(2k+2)!}
   * @f]
   * @f[
   *    hai(x) = \sum_{k=0}^\infty \frac{(2k+3)!!!x^{3k+2}}{(2k+3)!}
   * @f]
   *
   * @tparam Tp A real type
   */
  template<typename Tp>
    AiryAuxilliaryState<std::complex<Tp>>
    AirySeries<Tp>::s_FGH(std::complex<Tp> t)
    {
      const Tp s_log10min = emsr::log10_min(t);
      const auto s_min = emsr::lim_min(t);
      const auto log10t = std::log10(std::abs(t));
      const auto tt = t * t;
      const auto ttt = t * tt;

      auto Fai = Tp{1};
      auto Gai = Tp{1};
      auto Faip = Tp{0};
      auto Gaip = Tp{1};
      auto Hai = Tp{1} / Tp{2};
      auto Haip = Tp{1};
      auto term = Cmplx{Tp{1}};
      auto F = Cmplx{Tp{1}};
      auto G = t;
      auto H = t * t / Tp{2};
      auto Fp = Cmplx{Tp{0}};
      auto Gp = Cmplx{Tp{1}};
      auto Hp = t;
      for (int k = 1; k < max_FGH<Tp>; ++k)
	{
	  if (std::abs(t) < s_eps)
	    break;

	  auto xx = log10t * (3 * k + 2)
		    + s_slope_H * k + s_intercept_H;
	  if (xx < s_log10min)
	    break;

	  Fai /= (3ULL * k - 1ULL) * (3ULL * k);
	  if (Fai < Tp{10} * s_min)
	    break;
	  Faip = (3ULL * k) * Fai;
	  Gai /= (3ULL * k) * (3ULL * k + 1ULL);
	  Gaip = (3ULL * k + 1ULL) * Gai;
	  Hai /= (3ULL * k + 1ULL) * (3ULL * k + 2ULL);
	  Haip = (3ULL * k + 2ULL) * Hai;

	  term *= ttt;
	  F += Fai * term;
	  G += Gai * term * t;
	  H += Hai * term * tt;
	  Fp += Faip * term / t;
	  Gp += Gaip * term;
	  Hp += Haip * term * t;
	}

      return AiryAuxilliaryState<Cmplx>{t, F, G, H, Fp, Gp, Hp};
    }

  /**
   * Return the Scorer functions by using the series expansions of
   * the auxilliary Airy functions:
   * @f[
   *    fai(x) = \sum_{k=0}^\infty \frac{(2k+1)!!!x^{3k}}{(2k+1)!}
   * @f]
   * @f[
   *    gai(x) = \sum_{k=0}^\infty \frac{(2k+2)!!!x^{3k+1}}{(2k+2)!}
   * @f]
   * @f[
   *    hai(x) = \sum_{k=0}^\infty \frac{(2k+3)!!!x^{3k+2}}{(2k+3)!}
   * @f]
   * The Scorer function is then defined by:
   * @f[
   *    Hi(x) = Hi(0)\left(fai(x) + gai(x) + hai(x)\right)
   * @f]
   * where @f$ Hi(0) = 2/(3^{7/6}\Gamma(2/3)) @f$
   *   and @f$ Hi'(0) = 2/(3^{5/6}\Gamma(1/3)) @f$.
   * The other Scorer function is found from the identity
   * @f[
   *    Gi(x) + Hi(x) = Bi(x)
   * @f]
   *
   * @todo Find out what is wrong with the Hi = fai + gai + hai scorer function.
   *
   * @tparam Tp A real type
   */
  template<typename Tp>
    AiryState<std::complex<Tp>>
    AirySeries<Tp>::s_Scorer(std::complex<Tp> t)
    {
      auto aux = FGH(t);

      auto Hi = s_Hi0 * (aux.fai + aux.gai + aux.hai_value);
      auto Hip = s_Hip0 * (aux.faip + aux.gaip + aux.haip);
      auto Bi = s_Bi0 * aux.fai + s_Bip0 * aux.gai;
      auto Bip = s_Bi0 * aux.faip + s_Bip0 * aux.gaip;
      auto Gi = Bi - Hi;
      auto Gip = Bip - Hip;

      return AiryState<Cmplx>{t, Gi, Gip, Hi, Hip};
    }

  /**
   * Return the Scorer functions by using the series expansions:
   * @f[
   *    Hi(x) = \frac{3^{-2/3}}{\pi} \sum_{k=0}^\infty
   *                     \Gamma\left(\frac{k+1}{3}\right) \frac{3^{1/3}x}{k!}
   * @f]
   * @f[
   *    Hi'(x) = \frac{3^{-1/3}}{\pi} \sum_{k=0}^\infty
   *                     \Gamma\left(\frac{k+2}{3}\right) \frac{3^{1/3}x}{k!}
   * @f]
   * @f[
   *    Gi(x) = \frac{3^{-2/3}}{\pi} \sum_{k=0}^\infty
   *                     \cos\left(\frac{2k-1}{3}\pi\right)
   *                     \Gamma\left(\frac{k+1}{3}\right) \frac{3^{1/3}x}{k!}
   * @f]
   * @f[
   *    Gi'(x) = \frac{3^{-1/3}}{\pi} \sum_{k=0}^\infty
   *                     \cos\left(\frac{2k+1}{3}\pi\right)
   *                     \Gamma\left(\frac{k+2}{3}\right) \frac{3^{1/3}x}{k!}
   * @f]
   */
  template<typename Tp>
    AiryState<std::complex<Tp>>
    AirySeries<Tp>::s_Scorer2(std::complex<Tp> t)
    {
      const auto s_cbrt3 = emsr::cbrt3_v<Tp>;
      const auto s_1d3 = Tp{1} / Tp{3};
      const auto s_2d3 = Tp{2} / Tp{3};
      const auto s = s_cbrt3 * t;
      const std::array<Tp, 3>
	cos{ Tp{1} / Tp{2}, Tp{-1}, Tp{1} / Tp{2} };
      auto Hi = Cmplx{gamma(s_1d3)};
      auto Hip = Cmplx{gamma(s_2d3)};
      auto Gi = Cmplx{cos[2] * gamma(s_1d3)};
      auto Gip = Cmplx{cos[0] * gamma(s_2d3)};
      auto term = Cmplx(Tp{1});
      auto termp = Cmplx(Tp{1});
      for (int k = 1; k < max_FGH<Tp>; ++k)
	{
	  term *= s / Tp(k);
	  termp *= s / Tp(k);
	  if (std::abs(term) < s_eps)
	    break;

	  const auto gam = gamma(Tp(k + 1) /Tp{3});
	  const auto gamp = gamma(Tp(k + 2) /Tp{3});
	  Hi += gam * term;
	  Hip += gamp * termp;
	  Gi += cos[(k + 2) % 3] * gam * term;
	  Gip += cos[k % 3] * gamp * termp;
	}

      const auto fact = Tp{1} / (s_cbrt3 * s_cbrt3 * s_pi);
      const auto factp = Tp{1} / (s_cbrt3 * s_pi);
      Gi *= fact;
      Gip *= factp;
      Hi *= fact;
      Hip *= factp;

      return AiryState<Cmplx>{t, Gi, Gip, Hi, Hip};
    }

  /**
   * Return the Airy function of either the first or second kind and it's
   * derivative by series expansion as a pair of complex numbers.  The type
   * of function is determined by the input value of the function and it's
   * derivative at the origin.
   *
   * @tparam Tp A real type
   */
  template<typename Tp>
    std::pair<std::complex<Tp>, std::complex<Tp>>
    AirySeries<Tp>::s_AiBi(std::complex<Tp> t, std::pair<Tp, Tp> Z0)
    {
      const Tp s_log10min = emsr::log10_min(t.real());
      const auto log10t = std::log10(std::abs(t));
      const auto ttt = t * t * t;

      auto Fai = Tp{1};
      auto Gai = Tp{1};
      auto termF = Z0.first * Cmplx{Tp{1}};
      auto termG = Z0.second * t;
      auto Ai = termF + termG;
      if (std::abs(t) >= s_eps)
	for (int k = 1; k < max_FGH<Tp>; ++k)
	  {
	    auto xx = log10t * (3 * k + 1)
		      + s_slope_G * k + s_intercept_G;
	    if (xx < s_log10min)
	      break;

	    termF *= ttt;
	    termG *= ttt;
	    Ai += Fai * termF + Gai * termG;
	  }

      auto Faip = Tp{0};
      auto Gaip = Tp{1};
      termF = Z0.first * Cmplx{Tp{1}};
      termG = Z0.second * Cmplx{Tp{1}};
      auto Aip = termG;
      if (std::abs(t) >= s_eps)
	{
	  termF *= t * t;
	  termG *= ttt;
	  Aip += Faip * termF + Gaip * termG;
	  for (int k = 2; k < max_FGH<Tp>; ++k)
	    {
	      auto xx = log10t * 3 * k
			+ s_slope_Gp * k + s_intercept_Gp;
	      if (xx < s_log10min)
		break;

	      termF *= ttt;
	      termG *= ttt;
	      Aip += Faip * termF + Gaip * termG;
	    }
	}

      return std::make_pair(Ai, Aip);
    }

  /**
   * Return the Fock-type Airy functions @f$ Ai(t) @f$, and @f$ Bi(t) @f$
   * and their derivatives of complex argument.
   *
   * @tparam Tp A real type
   * @param t The complex argument
   */
  template<typename Tp>
    AiryState<std::complex<Tp>>
    AirySeries<Tp>::s_Airy(std::complex<Tp> t)
    { return s_AiryHelp(t, false); }

  /**
   * Return the Fock-type Airy functions @f$ w_1(t) @f$, and @f$ w_2(t) @f$
   * and their derivatives of complex argument.
   *
   * @tparam Tp A real type
   * @param t The complex argument
   */
  template<typename Tp>
    AiryState<std::complex<Tp>>
    AirySeries<Tp>::s_Fock(std::complex<Tp> t)
    { return s_AiryHelp(t, true); }


  /**
   * A class encapsulating data for the asymptotic expansions of Airy functions
   * and their derivatives.
   *
   * @tparam Tp A real type
   */
  template<typename Tp>
    struct Airy_asymp_data
    { };

  template<>
    struct Airy_asymp_data<float>
    {
      static constexpr int s_max_cd = 43;

      static constexpr float
      s_c[s_max_cd]
      {
	1.000000e+00F,
	6.944445e-02F,
	3.713349e-02F,
	3.799306e-02F,
	5.764920e-02F,
	1.160991e-01F,
	2.915914e-01F,
	8.776670e-01F,
	3.079453e+00F,
	1.234157e+01F,
	5.562279e+01F,
	2.784651e+02F,
	1.533170e+03F,
	9.207208e+03F,
	5.989252e+04F,
	4.195249e+05F,
	3.148258e+06F,
	2.519892e+07F,
	2.142880e+08F,
	1.929376e+09F,
	1.833577e+10F,
	1.834183e+11F,
	1.926471e+12F,
	2.119700e+13F,
	2.438269e+14F,
	2.926600e+15F,
	3.659031e+16F,
	4.757682e+17F,
	6.424051e+18F,
	8.995209e+19F,
	1.304514e+21F,
	1.957062e+22F,
	3.033872e+23F,
	4.854833e+24F,
	8.011465e+25F,
	1.362108e+27F,
	2.383952e+28F,
	4.291560e+29F,
	7.940171e+30F,
	1.508774e+32F,
	2.942371e+33F,
	5.885241e+34F,
	1.206572e+36F,
      };

      static constexpr float
      s_d[s_max_cd]
      {
	-1.000000e+00F,
	-9.722223e-02F,
	-4.388503e-02F,
	-4.246284e-02F,
	-6.266217e-02F,
	-1.241059e-01F,
	-3.082538e-01F,
	-9.204800e-01F,
	-3.210494e+00F,
	-1.280729e+01F,
	-5.750832e+01F,
	-2.870333e+02F,
	-1.576357e+03F,
	-9.446356e+03F,
	-6.133572e+04F,
	-4.289525e+05F,
	-3.214537e+06F,
	-2.569791e+07F,
	-2.182934e+08F,
	-1.963524e+09F,
	-1.864393e+10F,
	-1.863530e+11F,
	-1.955883e+12F,
	-2.150644e+13F,
	-2.472370e+14F,
	-2.965883e+15F,
	-3.706245e+16F,
	-4.816783e+17F,
	-6.500985e+18F,
	-9.099199e+19F,
	-1.319089e+21F,
	-1.978220e+22F,
	-3.065640e+23F,
	-4.904121e+24F,
	-8.090396e+25F,
	-1.375143e+27F,
	-2.406128e+28F,
	-4.330398e+29F,
	-8.010129e+30F,
	-1.521725e+32F,
	-2.966994e+33F,
	-5.933284e+34F,
	-1.216186e+36F,
      };
    };


  template<>
    struct Airy_asymp_data<double>
    { 
      static constexpr int s_max_cd = 198;

      static constexpr double
      s_c[s_max_cd]
      {
	1.000000000000000e+00,
	6.944444444444445e-02,
	3.713348765432099e-02,
	3.799305912780064e-02,
	5.764919041266972e-02,
	1.160990640255154e-01,
	2.915913992307505e-01,
	8.776669695100169e-01,
	3.079453030173167e+00,
	1.234157333234524e+01,
	5.562278536591707e+01,
	2.784650807776025e+02,
	1.533169432012795e+03,
	9.207206599726414e+03,
	5.989251356587907e+04,
	4.195248751165511e+05,
	3.148257417866826e+06,
	2.519891987160237e+07,
	2.142880369636803e+08,
	1.929375549182493e+09,
	1.833576693789056e+10,
	1.834183035288325e+11,
	1.926471158970446e+12,
	2.119699938864764e+13,
	2.438268268797160e+14,
	2.926599219297924e+15,
	3.659030701264312e+16,
	4.757681020363067e+17,
	6.424049357901937e+18,
	8.995207427058378e+19,
	1.304513299317610e+21,
	1.957062178658161e+22,
	3.033871086594339e+23,
	4.854832179436168e+24,
	8.011464687609595e+25,
	1.362107954526322e+27,
	2.383951672727106e+28,
	4.291560449285805e+29,
	7.940171107576634e+30,
	1.508773895252729e+32,
	2.942371035655193e+33,
	5.885240440388241e+34,
	1.206571599149305e+36,
	2.533995217967925e+37,
	5.448489654744985e+38,
	1.198751805674371e+40,
	2.697372533752491e+41,
	6.204355375581934e+42,
	1.458113275347628e+44,
	3.499678509541131e+45,
	8.574698414835430e+46,
	2.143791361584876e+48,
	5.466954268964723e+49,
	1.421479741931207e+51,
	3.767104119582454e+52,
	1.017165676733217e+54,
	2.797331747633005e+55,
	7.832869698897223e+56,
	2.232461648545130e+58,
	6.474401546982448e+59,
	1.910023391562912e+61,
	5.730287618153169e+62,
	1.747801906865773e+64,
	5.418378570224247e+65,
	1.706848042790888e+67,
	5.462096092490968e+68,
	1.775238701609359e+70,
	5.858471716005496e+71,
	1.962647854025608e+73,
	6.673200232657881e+74,
	2.302320282650229e+76,
	8.058346177096876e+77,
	2.860790616115698e+79,
	1.029911836324255e+81,
	3.759274853469072e+82,
	1.390966503884052e+84,
	5.216251488112698e+85,
	1.982222609597978e+87,
	7.631733526885406e+88,
	2.976443161750489e+90,
	1.175720886071667e+92,
	4.702984343402413e+93,
	1.904748487874923e+95,
	7.809628166793868e+96,
	3.241060252944379e+98,
	1.361271785487072e+100,
	5.785515010137358e+101,
	2.487817635034432e+103,
	1.082220303639244e+105,
	4.761853778920262e+106,
	2.119061674318428e+108,
	9.535939245488659e+109,
	4.338924336915075e+111,
	1.995937594356210e+113,
	9.281257267774871e+114,
	4.362258761302053e+116,
	2.072104467309746e+118,
	9.946249789626537e+119,
	4.824001628763865e+121,
	2.363794636489558e+123,
	1.170094760302862e+125,
	5.850554253574287e+126,
	2.954569730259901e+128,
	1.506850482670752e+130,
	7.760380603441518e+131,
	4.035449239058130e+133,
	2.118637288197074e+135,
	1.122891512986455e+137,
	6.007541796863912e+138,
	3.244110844655371e+140,
	1.768060890834934e+142,
	9.724445514012234e+143,
	5.397127555697884e+145,
	3.022424599378842e+147,
	1.707688310104938e+149,
	9.733926488872914e+150,
	5.597066004129277e+152,
	3.246331503347047e+154,
	1.899123034516305e+156,
	1.120493673015381e+158,
	6.667002197825375e+159,
	4.000239582022806e+161,
	2.420167717157848e+163,
	1.476315971466497e+165,
	9.079425903504815e+166,
	5.629294501428005e+168,
	3.518340089045620e+170,
	2.216573494616288e+172,
	1.407536194762195e+174,
	9.008307418237021e+175,
	5.810406406063190e+177,
	3.776794965501750e+179,
	2.473820571905779e+181,
	1.632734494231895e+183,
	1.085776900182112e+185,
	7.274761083941352e+186,
	4.910500878112026e+188,
	3.339165488138471e+190,
	2.287345162743855e+192,
	1.578279589877007e+194,
	1.096912143732327e+196,
	7.678439030562050e+197,
	5.413337067597843e+199,
	3.843495606538892e+201,
	2.748117894051496e+203,
	1.978658045201244e+205,
	1.434536494196155e+207,
	1.047218417674069e+209,
	7.697104507405235e+210,
	5.695893209382314e+212,
	4.243466810865792e+214,
	3.182619623725185e+216,
	2.402892356389597e+218,
	1.826209097230261e+220,
	1.397058194451032e+222,
	1.075741068948596e+224,
	8.337041171685528e+225,
	6.502928990423786e+227,
	5.104827839273243e+229,
	4.032836288744937e+231,
	3.206122353181950e+233,
	2.564911711575723e+235,
	2.064764922810362e+237,
	1.672468384191300e+239,
	1.363068815045043e+241,
	1.117722165158547e+243,
	9.221254621350063e+244,
	7.653679680924400e+246,
	6.390854170806005e+248,
	5.368343764383400e+250,
	4.536272410412814e+252,
	3.855849971009803e+254,
	3.296767293083097e+256,
	2.835233105703090e+258,
	2.452487952018684e+260,
	2.133674250351303e+262,
	1.866973387910966e+264,
	1.642943906272933e+266,
	1.454011766789008e+268,
	1.294076113394137e+270,
	1.158203114065351e+272,
	1.042387246347866e+274,
	9.433644353076327e+275,
	8.584652159889140e+277,
	7.854989126102911e+279,
	7.226619481709595e+281,
	6.684650001687565e+283,
	6.216749325730167e+285,
	5.812683583318513e+287,
	5.463943925916343e+289,
	5.163446980546229e+291,
	4.905293404959073e+293,
	4.684572943682547e+295,
	4.497206881767656e+297,
	4.339820739154935e+299,
	4.209641572182349e+301,
	4.104415447991071e+303,
	4.022341607500821e+305,
      };

      static constexpr double
      s_d[s_max_cd]
      {
	-1.000000000000000e+00,
	-9.722222222222224e-02,
	-4.388503086419752e-02,
	-4.246283078989483e-02,
	-6.266216349203231e-02,
	-1.241058960272751e-01,
	-3.082537649010791e-01,
	-9.204799924129445e-01,
	-3.210493584648621e+00,
	-1.280729308073562e+01,
	-5.750830351391426e+01,
	-2.870332371092210e+02,
	-1.576357303337099e+03,
	-9.446354823095931e+03,
	-6.133570666385206e+04,
	-4.289524004000691e+05,
	-3.214536521400865e+06,
	-2.569790838391133e+07,
	-2.182934208321604e+08,
	-1.963523788991033e+09,
	-1.864393108810722e+10,
	-1.863529963852939e+11,
	-1.955882932389843e+12,
	-2.150644463519725e+13,
	-2.472369922906211e+14,
	-2.965882430295212e+15,
	-3.706244000635465e+16,
	-4.816782647945217e+17,
	-6.500984080751062e+18,
	-9.099198264365413e+19,
	-1.319088866907751e+21,
	-1.978219607616628e+22,
	-3.065639370223599e+23,
	-4.904119815775620e+24,
	-8.090395374187028e+25,
	-1.375142480406956e+27,
	-2.406127967357125e+28,
	-4.330398100410563e+29,
	-8.010128562268939e+30,
	-1.521724744139019e+32,
	-2.966993387417997e+33,
	-5.933283219493451e+34,
	-1.216185715477188e+36,
	-2.553715025111644e+37,
	-5.489923036149890e+38,
	-1.207664458504664e+40,
	-2.716989788543418e+41,
	-6.248514488575399e+42,
	-1.468274343468517e+44,
	-3.523567100049944e+45,
	-8.632054257075131e+46,
	-2.157849009857564e+48,
	-5.502111531144559e+49,
	-1.430448068378723e+51,
	-3.790429841685131e+52,
	-1.023349054707279e+54,
	-2.814032235678575e+55,
	-7.878810283641486e+56,
	-2.245328862657782e+58,
	-6.511083708721726e+59,
	-1.920664190401703e+61,
	-5.761686454417021e+62,
	-1.757224019571249e+64,
	-5.447123284124641e+65,
	-1.715761087400762e+67,
	-5.490178848750562e+68,
	-1.784227251997255e+70,
	-5.887691026309764e+71,
	-1.972292315224751e+73,
	-6.705515972283344e+74,
	-2.313309878271472e+76,
	-8.096267806165569e+77,
	-2.874065746584913e+79,
	-1.034625391639241e+81,
	-3.776246748970062e+82,
	-1.397162345772177e+84,
	-5.239180066082425e+85,
	-1.990822273847861e+87,
	-7.664417610512326e+88,
	-2.989028545098271e+90,
	-1.180629950314137e+92,
	-4.722378093272115e+93,
	-1.912507137520035e+95,
	-7.841055241911753e+96,
	-3.253947172439189e+98,
	-1.366620594074447e+100,
	-5.807983029594204e+101,
	-2.497367798700592e+103,
	-1.086327401565769e+105,
	-4.779721898165742e+106,
	-2.126924611885473e+108,
	-9.570933517949168e+109,
	-4.354673608555420e+111,
	-2.003104336167184e+113,
	-9.314227986310485e+114,
	-4.377591832519284e+116,
	-2.079311787196041e+118,
	-9.980488171002189e+119,
	-4.840437750156587e+121,
	-2.371766962413637e+123,
	-1.174001587549283e+125,
	-5.869894928792715e+126,
	-2.964240989606087e+128,
	-1.511734925078113e+130,
	-7.785293542778409e+131,
	-4.048280556193450e+133,
	-2.125310161545726e+135,
	-1.126395074649439e+137,
	-6.026112250640926e+138,
	-3.254046865619094e+140,
	-1.773426781247180e+142,
	-9.753691966685956e+143,
	-5.413214374045717e+145,
	-3.031353475595618e+147,
	-1.712688861525451e+149,
	-9.762181718158467e+150,
	-5.613172668889362e+152,
	-3.255593504783130e+154,
	-1.904495376905319e+156,
	-1.123636712771385e+158,
	-6.685547405607922e+159,
	-4.011274725697352e+161,
	-2.426789243059785e+163,
	-1.480322256328008e+165,
	-9.103865811724210e+166,
	-5.644325995423808e+168,
	-3.527660195241767e+170,
	-2.222398917729603e+172,
	-1.411206432558185e+174,
	-9.031614811298435e+175,
	-5.825324009159629e+177,
	-3.786417373057168e+179,
	-2.480075491177349e+181,
	-1.636831694970243e+183,
	-1.088481201303362e+185,
	-7.292745660168647e+186,
	-4.922551187015368e+188,
	-3.347299874224070e+190,
	-2.292876831819414e+192,
	-1.582068976647419e+194,
	-1.099526952179841e+196,
	-7.696612850752728e+197,
	-5.426059363878918e+199,
	-3.852465257896041e+201,
	-2.754486649310016e+203,
	-1.983211918722996e+205,
	-1.437815434754317e+207,
	-1.049595758009311e+209,
	-7.714459872698144e+210,
	-5.708649969089776e+212,
	-4.252907226462823e+214,
	-3.189653037258280e+216,
	-2.408167641474974e+218,
	-1.830192105075910e+220,
	-1.400085406139984e+222,
	-1.078056980830078e+224,
	-8.354874414833517e+225,
	-6.516750306025111e+227,
	-5.115608890676142e+229,
	-4.041299743705576e+231,
	-3.212808739737073e+233,
	-2.570227590770699e+235,
	-2.069017785679178e+237,
	-1.675892065632837e+239,
	-1.365842098493761e+241,
	-1.119982472873731e+243,
	-9.239789806518604e+244,
	-7.668971748218954e+246,
	-6.403547029139283e+248,
	-5.378942666188698e+250,
	-4.545175791002436e+252,
	-3.863373580709335e+254,
	-3.303162573962017e+256,
	-2.840701250554976e+258,
	-2.457190709357166e+260,
	-2.137742266081238e+262,
	-1.870512673954399e+264,
	-1.646040878763552e+266,
	-1.456737187157872e+268,
	-1.296488184434946e+270,
	-1.160349922432478e+272,
	-1.044308697493207e+274,
	-9.450937926592507e+275,
	-8.600303303298236e+277,
	-7.869232080094031e+279,
	-7.239652158863265e+281,
	-6.696640405278035e+283,
	-6.227840760744850e+285,
	-5.822998904673116e+287,
	-5.473589016694659e+289,
	-5.172513612645520e+291,
	-4.913861603046338e+293,
	-4.692712948797547e+295,
	-4.504980791675638e+297,
	-4.347283887459587e+299,
	-4.216843696343482e+301,
	-4.111401687051481e+303,
	-4.029153362974997e+305,
      };
    };


  template<>
    struct Airy_asymp_data<long double>
    {
      static constexpr int s_max_cd = 201;

      static constexpr long double
      s_c[s_max_cd]
      {
	1.000000000000000000e+00L,
	6.944444444444444445e-02L,
	3.713348765432098766e-02L,
	3.799305912780064015e-02L,
	5.764919041266972134e-02L,
	1.160990640255154110e-01L,
	2.915913992307505115e-01L,
	8.776669695100169165e-01L,
	3.079453030173166994e+00L,
	1.234157333234523871e+01L,
	5.562278536591708279e+01L,
	2.784650807776025672e+02L,
	1.533169432012795616e+03L,
	9.207206599726414699e+03L,
	5.989251356587906863e+04L,
	4.195248751165510687e+05L,
	3.148257417866826379e+06L,
	2.519891987160236768e+07L,
	2.142880369636803196e+08L,
	1.929375549182493053e+09L,
	1.833576693789056766e+10L,
	1.834183035288325634e+11L,
	1.926471158970446564e+12L,
	2.119699938864764906e+13L,
	2.438268268797160419e+14L,
	2.926599219297925047e+15L,
	3.659030701264312806e+16L,
	4.757681020363067633e+17L,
	6.424049357901937700e+18L,
	8.995207427058378954e+19L,
	1.304513299317609818e+21L,
	1.957062178658161504e+22L,
	3.033871086594338300e+23L,
	4.854832179436167361e+24L,
	8.011464687609593663e+25L,
	1.362107954526321589e+27L,
	2.383951672727105667e+28L,
	4.291560449285803547e+29L,
	7.940171107576632359e+30L,
	1.508773895252729248e+32L,
	2.942371035655192299e+33L,
	5.885240440388239474e+34L,
	1.206571599149304506e+36L,
	2.533995217967924041e+37L,
	5.448489654744983645e+38L,
	1.198751805674370862e+40L,
	2.697372533752490593e+41L,
	6.204355375581932929e+42L,
	1.458113275347627820e+44L,
	3.499678509541130410e+45L,
	8.574698414835427995e+46L,
	2.143791361584876000e+48L,
	5.466954268964722378e+49L,
	1.421479741931207335e+51L,
	3.767104119582453966e+52L,
	1.017165676733216895e+54L,
	2.797331747633004845e+55L,
	7.832869698897222655e+56L,
	2.232461648545130119e+58L,
	6.474401546982448099e+59L,
	1.910023391562912264e+61L,
	5.730287618153167906e+62L,
	1.747801906865772587e+64L,
	5.418378570224246184e+65L,
	1.706848042790887377e+67L,
	5.462096092490966838e+68L,
	1.775238701609358950e+70L,
	5.858471716005495788e+71L,
	1.962647854025607485e+73L,
	6.673200232657881231e+74L,
	2.302320282650229520e+76L,
	8.058346177096876306e+77L,
	2.860790616115697965e+79L,
	1.029911836324255108e+81L,
	3.759274853469072083e+82L,
	1.390966503884051755e+84L,
	5.216251488112698105e+85L,
	1.982222609597977811e+87L,
	7.631733526885405277e+88L,
	2.976443161750488134e+90L,
	1.175720886071666341e+92L,
	4.702984343402412251e+93L,
	1.904748487874923120e+95L,
	7.809628166793867769e+96L,
	3.241060252944379015e+98L,
	1.361271785487071739e+100L,
	5.785515010137358431e+101L,
	2.487817635034432399e+103L,
	1.082220303639244464e+105L,
	4.761853778920262973e+106L,
	2.119061674318428445e+108L,
	9.535939245488660482e+109L,
	4.338924336915075186e+111L,
	1.995937594356210236e+113L,
	9.281257267774873180e+114L,
	4.362258761302054246e+116L,
	2.072104467309746760e+118L,
	9.946249789626540298e+119L,
	4.824001628763866564e+121L,
	2.363794636489558111e+123L,
	1.170094760302862443e+125L,
	5.850554253574289032e+126L,
	2.954569730259901407e+128L,
	1.506850482670752657e+130L,
	7.760380603441520329e+131L,
	4.035449239058131851e+133L,
	2.118637288197074896e+135L,
	1.122891512986455128e+137L,
	6.007541796863915236e+138L,
	3.244110844655372943e+140L,
	1.768060890834934916e+142L,
	9.724445514012239322e+143L,
	5.397127555697886649e+145L,
	3.022424599378843182e+147L,
	1.707688310104939300e+149L,
	9.733926488872918803e+150L,
	5.597066004129280281e+152L,
	3.246331503347048857e+154L,
	1.899123034516305991e+156L,
	1.120493673015382185e+158L,
	6.667002197825379055e+159L,
	4.000239582022809076e+161L,
	2.420167717157850167e+163L,
	1.476315971466498391e+165L,
	9.079425903504822323e+166L,
	5.629294501428009312e+168L,
	3.518340089045622156e+170L,
	2.216573494616289233e+172L,
	1.407536194762195531e+174L,
	9.008307418237028691e+175L,
	5.810406406063194601e+177L,
	3.776794965501753246e+179L,
	2.473820571905781698e+181L,
	1.632734494231896465e+183L,
	1.085776900182112434e+185L,
	7.274761083941356088e+186L,
	4.910500878112028622e+188L,
	3.339165488138473218e+190L,
	2.287345162743856462e+192L,
	1.578279589877007681e+194L,
	1.096912143732327268e+196L,
	7.678439030562053892e+197L,
	5.413337067597845030e+199L,
	3.843495606538893038e+201L,
	2.748117894051497464e+203L,
	1.978658045201245087e+205L,
	1.434536494196155738e+207L,
	1.047218417674069521e+209L,
	7.697104507405240282e+210L,
	5.695893209382317244e+212L,
	4.243466810865795709e+214L,
	3.182619623725187776e+216L,
	2.402892356389598849e+218L,
	1.826209097230261744e+220L,
	1.397058194451033127e+222L,
	1.075741068948596812e+224L,
	8.337041171685535040e+225L,
	6.502928990423792373e+227L,
	5.104827839273246737e+229L,
	4.032836288744939317e+231L,
	3.206122353181952212e+233L,
	2.564911711575725598e+235L,
	2.064764922810364009e+237L,
	1.672468384191301163e+239L,
	1.363068815045044049e+241L,
	1.117722165158548263e+243L,
	9.221254621350072972e+244L,
	7.653679680924408496e+246L,
	6.390854170806011899e+248L,
	5.368343764383406036e+250L,
	4.536272410412819536e+252L,
	3.855849971009808288e+254L,
	3.296767293083101081e+256L,
	2.835233105703092986e+258L,
	2.452487952018686700e+260L,
	2.133674250351305122e+262L,
	1.866973387910968178e+264L,
	1.642943906272935389e+266L,
	1.454011766789009870e+268L,
	1.294076113394138357e+270L,
	1.158203114065351801e+272L,
	1.042387246347866902e+274L,
	9.433644353076337553e+275L,
	8.584652159889149460e+277L,
	7.854989126102919164e+279L,
	7.226619481709603434e+281L,
	6.684650001687572355e+283L,
	6.216749325730173453e+285L,
	5.812683583318519397e+287L,
	5.463943925916348704e+289L,
	5.163446980546234308e+291L,
	4.905293404959078738e+293L,
	4.684572943682552715e+295L,
	4.497206881767661669e+297L,
	4.339820739154940732e+299L,
	4.209641572182355027e+301L,
	4.104415447991076276e+303L,
	4.022341607500826067e+305L,
	3.962020590927621690e+307L,
	3.922414211165019564e+309L,
	3.902815759602983234e+311L,
      };

      static constexpr long double
      s_d[s_max_cd]
      {
	-1.000000000000000000e+00L,
	-9.722222222222222222e-02L,
	-4.388503086419753087e-02L,
	-4.246283078989483311e-02L,
	-6.266216349203230580e-02L,
	-1.241058960272750945e-01L,
	-3.082537649010791121e-01L,
	-9.204799924129445710e-01L,
	-3.210493584648620908e+00L,
	-1.280729308073562507e+01L,
	-5.750830351391427204e+01L,
	-2.870332371092211078e+02L,
	-1.576357303337099718e+03L,
	-9.446354823095931964e+03L,
	-6.133570666385205824e+04L,
	-4.289524004000690703e+05L,
	-3.214536521400864830e+06L,
	-2.569790838391132545e+07L,
	-2.182934208321603255e+08L,
	-1.963523788991032753e+09L,
	-1.864393108810721585e+10L,
	-1.863529963852938844e+11L,
	-1.955882932389842695e+12L,
	-2.150644463519724977e+13L,
	-2.472369922906211613e+14L,
	-2.965882430295212631e+15L,
	-3.706244000635465229e+16L,
	-4.816782647945217541e+17L,
	-6.500984080751062702e+18L,
	-9.099198264365412236e+19L,
	-1.319088866907750710e+21L,
	-1.978219607616628114e+22L,
	-3.065639370223598387e+23L,
	-4.904119815775620837e+24L,
	-8.090395374187028084e+25L,
	-1.375142480406956246e+27L,
	-2.406127967357125255e+28L,
	-4.330398100410561950e+29L,
	-8.010128562268937490e+30L,
	-1.521724744139018769e+32L,
	-2.966993387417997255e+33L,
	-5.933283219493449593e+34L,
	-1.216185715477187411e+36L,
	-2.553715025111643294e+37L,
	-5.489923036149888465e+38L,
	-1.207664458504663582e+40L,
	-2.716989788543417798e+41L,
	-6.248514488575398644e+42L,
	-1.468274343468517212e+44L,
	-3.523567100049943587e+45L,
	-8.632054257075129856e+46L,
	-2.157849009857563711e+48L,
	-5.502111531144559821e+49L,
	-1.430448068378722839e+51L,
	-3.790429841685131699e+52L,
	-1.023349054707279003e+54L,
	-2.814032235678575023e+55L,
	-7.878810283641487892e+56L,
	-2.245328862657782166e+58L,
	-6.511083708721725425e+59L,
	-1.920664190401702861e+61L,
	-5.761686454417020881e+62L,
	-1.757224019571248450e+64L,
	-5.447123284124640063e+65L,
	-1.715761087400761462e+67L,
	-5.490178848750560497e+68L,
	-1.784227251997254439e+70L,
	-5.887691026309762600e+71L,
	-1.972292315224750519e+73L,
	-6.705515972283343125e+74L,
	-2.313309878271471665e+76L,
	-8.096267806165567488e+77L,
	-2.874065746584912340e+79L,
	-1.034625391639240257e+81L,
	-3.776246748970061121e+82L,
	-1.397162345772176707e+84L,
	-5.239180066082424250e+85L,
	-1.990822273847860578e+87L,
	-7.664417610512323500e+88L,
	-2.989028545098270324e+90L,
	-1.180629950314136764e+92L,
	-4.722378093272112919e+93L,
	-1.912507137520034823e+95L,
	-7.841055241911750536e+96L,
	-3.253947172439187679e+98L,
	-1.366620594074447266e+100L,
	-5.807983029594202542e+101L,
	-2.497367798700591449e+103L,
	-1.086327401565769111e+105L,
	-4.779721898165742383e+106L,
	-2.126924611885472706e+108L,
	-9.570933517949169328e+109L,
	-4.354673608555420287e+111L,
	-2.003104336167184061e+113L,
	-9.314227986310485517e+114L,
	-4.377591832519284665e+116L,
	-2.079311787196041531e+118L,
	-9.980488171002191038e+119L,
	-4.840437750156588426e+121L,
	-2.371766962413637565e+123L,
	-1.174001587549282684e+125L,
	-5.869894928792716434e+126L,
	-2.964240989606087664e+128L,
	-1.511734925078113281e+130L,
	-7.785293542778411245e+131L,
	-4.048280556193451825e+133L,
	-2.125310161545727101e+135L,
	-1.126395074649439387e+137L,
	-6.026112250640928884e+138L,
	-3.254046865619095371e+140L,
	-1.773426781247180546e+142L,
	-9.753691966685960342e+143L,
	-5.413214374045719396e+145L,
	-3.031353475595619676e+147L,
	-1.712688861525451568e+149L,
	-9.762181718158471542e+150L,
	-5.613172668889364540e+152L,
	-3.255593504783131735e+154L,
	-1.904495376905319587e+156L,
	-1.123636712771386062e+158L,
	-6.685547405607925311e+159L,
	-4.011274725697354756e+161L,
	-2.426789243059786829e+163L,
	-1.480322256328008563e+165L,
	-9.103865811724216192e+166L,
	-5.644325995423811740e+168L,
	-3.527660195241769499e+170L,
	-2.222398917729604053e+172L,
	-1.411206432558185611e+174L,
	-9.031614811298444031e+175L,
	-5.825324009159634125e+177L,
	-3.786417373057171726e+179L,
	-2.480075491177351310e+181L,
	-1.636831694970245013e+183L,
	-1.088481201303363025e+185L,
	-7.292745660168652396e+186L,
	-4.922551187015371024e+188L,
	-3.347299874224072422e+190L,
	-2.292876831819415970e+192L,
	-1.582068976647420664e+194L,
	-1.099526952179841755e+196L,
	-7.696612850752733309e+197L,
	-5.426059363878921047e+199L,
	-3.852465257896043314e+201L,
	-2.754486649310017736e+203L,
	-1.983211918722997090e+205L,
	-1.437815434754318380e+207L,
	-1.049595758009311450e+209L,
	-7.714459872698149505e+210L,
	-5.708649969089780440e+212L,
	-4.252907226462827513e+214L,
	-3.189653037258282114e+216L,
	-2.408167641474976673e+218L,
	-1.830192105075911170e+220L,
	-1.400085406139984445e+222L,
	-1.078056980830079259e+224L,
	-8.354874414833525489e+225L,
	-6.516750306025118181e+227L,
	-5.115608890676146942e+229L,
	-4.041299743705579274e+231L,
	-3.212808739737076200e+233L,
	-2.570227590770701195e+235L,
	-2.069017785679180413e+237L,
	-1.675892065632839139e+239L,
	-1.365842098493762348e+241L,
	-1.119982472873732385e+243L,
	-9.239789806518615833e+244L,
	-7.668971748218962759e+246L,
	-6.403547029139290970e+248L,
	-5.378942666188703975e+250L,
	-4.545175791002442342e+252L,
	-3.863373580709339622e+254L,
	-3.303162573962020773e+256L,
	-2.840701250554979375e+258L,
	-2.457190709357169321e+260L,
	-2.137742266081240880e+262L,
	-1.870512673954401294e+264L,
	-1.646040878763553552e+266L,
	-1.456737187157873994e+268L,
	-1.296488184434947562e+270L,
	-1.160349922432479422e+272L,
	-1.044308697493208592e+274L,
	-9.450937926592517823e+275L,
	-8.600303303298245448e+277L,
	-7.869232080094039597e+279L,
	-7.239652158863272692e+281L,
	-6.696640405278043338e+283L,
	-6.227840760744857081e+285L,
	-5.822998904673121916e+287L,
	-5.473589016694665295e+289L,
	-5.172513612645525325e+291L,
	-4.913861603046343504e+293L,
	-4.692712948797552807e+295L,
	-4.504980791675643798e+297L,
	-4.347283887459592393e+299L,
	-4.216843696343488227e+301L,
	-4.111401687051486619e+303L,
	-4.029153362975001895e+305L,
	-3.968696278528173706e+307L,
	-3.928989926523217418e+309L,
	-3.909325877634014065e+311L,
      };
    };


#ifdef EMSR_HAVE_FLOAT128
  template<>
    struct Airy_asymp_data<__float128>
    {
      static constexpr int s_max_cd = 201;

      static constexpr __float128
      s_c[s_max_cd]
      {
	1.000000000000000000000000000000000e+00Q,
	6.944444444444444444444444444444445e-02Q,
	3.713348765432098765432098765432099e-02Q,
	3.799305912780064014631915866483767e-02Q,
	5.764919041266972133313011227963217e-02Q,
	1.160990640255154110181092538964814e-01Q,
	2.915913992307505114690938436983388e-01Q,
	8.776669695100169164655066684332937e-01Q,
	3.079453030173166993362480862680011e+00Q,
	1.234157333234523870642339938330245e+01Q,
	5.562278536591708278103323749835620e+01Q,
	2.784650807776025672055514983345737e+02Q,
	1.533169432012795615968528330529591e+03Q,
	9.207206599726414698033224087507302e+03Q,
	5.989251356587906862599588327558073e+04Q,
	4.195248751165510686626470897960817e+05Q,
	3.148257417866826378983145912575630e+06Q,
	2.519891987160236767557016381168582e+07Q,
	2.142880369636803195620823884016894e+08Q,
	1.929375549182493052665328054052345e+09Q,
	1.833576693789056765675348223590718e+10Q,
	1.834183035288325633653415468373651e+11Q,
	1.926471158970446563579032395664927e+12Q,
	2.119699938864764905493571816510304e+13Q,
	2.438268268797160418199984201202274e+14Q,
	2.926599219297925046400592148165286e+15Q,
	3.659030701264312805075099317724813e+16Q,
	4.757681020363067632401403572743319e+17Q,
	6.424049357901937699484057869724498e+18Q,
	8.995207427058378952098438694307239e+19Q,
	1.304513299317609817937424037496177e+21Q,
	1.957062178658161503299043185284924e+22Q,
	3.033871086594338299189753708716216e+23Q,
	4.854832179436167359995522969659059e+24Q,
	8.011464687609593661835749240413277e+25Q,
	1.362107954526321589052986810339313e+27Q,
	2.383951672727105666951726336845792e+28Q,
	4.291560449285803546171319066670932e+29Q,
	7.940171107576632357848623628433817e+30Q,
	1.508773895252729247570260010478430e+32Q,
	2.942371035655192298256376857240314e+33Q,
	5.885240440388239473257038331156990e+34Q,
	1.206571599149304506030147504684986e+36Q,
	2.533995217967924040264412819835998e+37Q,
	5.448489654744983644276862627802042e+38Q,
	1.198751805674370861365049853157002e+40Q,
	2.697372533752490593092703511670458e+41Q,
	6.204355375581932928326145485753351e+42Q,
	1.458113275347627819362204463247607e+44Q,
	3.499678509541130409867726398367846e+45Q,
	8.574698414835427994510633526889885e+46Q,
	2.143791361584875999553509561372892e+48Q,
	5.466954268964722377387030779679476e+49Q,
	1.421479741931207334923266721673421e+51Q,
	3.767104119582453965238905174924491e+52Q,
	1.017165676733216894528206908658589e+54Q,
	2.797331747633004844984132516532928e+55Q,
	7.832869698897222655239804005906784e+56Q,
	2.232461648545130118891994043551321e+58Q,
	6.474401546982448099273469637085500e+59Q,
	1.910023391562912263916452171061022e+61Q,
	5.730287618153167906533500007282686e+62Q,
	1.747801906865772586581276345187251e+64Q,
	5.418378570224246183810530443928050e+65Q,
	1.706848042790887376960471283961076e+67Q,
	5.462096092490966837830875680393816e+68Q,
	1.775238701609358950481459584298028e+70Q,
	5.858471716005495788553004175818681e+71Q,
	1.962647854025607485159829039873276e+73Q,
	6.673200232657881231295145771503102e+74Q,
	2.302320282650229519650310362059160e+76Q,
	8.058346177096876306796396352336370e+77Q,
	2.860790616115697964648704466618304e+79Q,
	1.029911836324255107627078591853086e+81Q,
	3.759274853469072082695033095116603e+82Q,
	1.390966503884051754597927939944024e+84Q,
	5.216251488112698104939104585457117e+85Q,
	1.982222609597977811425425210588653e+87Q,
	7.631733526885405277140613683381215e+88Q,
	2.976443161750488133728977351758507e+90Q,
	1.175720886071666341365773562428566e+92Q,
	4.702984343402412250437491178174693e+93Q,
	1.904748487874923120067516222987124e+95Q,
	7.809628166793867767772136754895049e+96Q,
	3.241060252944379014728014234570755e+98Q,
	1.361271785487071738301899377328823e+100Q,
	5.785515010137358430397599315806522e+101Q,
	2.487817635034432399004799104258776e+103Q,
	1.082220303639244463772256524486913e+105Q,
	4.761853778920262972406899701541145e+106Q,
	2.119061674318428444503264864868692e+108Q,
	9.535939245488660481067843850461071e+109Q,
	4.338924336915075185707133051022668e+111Q,
	1.995937594356210235628280808705277e+113Q,
	9.281257267774873178578984209120842e+114Q,
	4.362258761302054245344246693668441e+116Q,
	2.072104467309746759407768801647412e+118Q,
	9.946249789626540297588519326430030e+119Q,
	4.824001628763866563465193696726490e+121Q,
	2.363794636489558111233771687446041e+123Q,
	1.170094760302862442547600553144731e+125Q,
	5.850554253574289032212875400525164e+126Q,
	2.954569730259901407135096119172891e+128Q,
	1.506850482670752656564243782753088e+130Q,
	7.760380603441520328785594012413897e+131Q,
	4.035449239058131850864863288532028e+133Q,
	2.118637288197074896047939919766058e+135Q,
	1.122891512986455127634484273389301e+137Q,
	6.007541796863915234743129911261451e+138Q,
	3.244110844655372941817576741848712e+140Q,
	1.768060890834934915998545344362736e+142Q,
	9.724445514012239318704962265951030e+143Q,
	5.397127555697886647769886900870839e+145Q,
	3.022424599378843181346562978935192e+147Q,
	1.707688310104939299770635973954737e+149Q,
	9.733926488872918799713142602023519e+150Q,
	5.597066004129280278937066773118756e+152Q,
	3.246331503347048856182115401407294e+154Q,
	1.899123034516305990523176999673793e+156Q,
	1.120493673015382184186863757404408e+158Q,
	6.667002197825379052751202022282931e+159Q,
	4.000239582022809074092484018707352e+161Q,
	2.420167717157850166107112283071242e+163Q,
	1.476315971466498390452292575285715e+165Q,
	9.079425903504822320053466065861946e+166Q,
	5.629294501428009309668260368979217e+168Q,
	3.518340089045622154572808609472069e+170Q,
	2.216573494616289232217672643612031e+172Q,
	1.407536194762195530346851781865534e+174Q,
	9.008307418237028687260890747812261e+175Q,
	5.810406406063194598521988532823673e+177Q,
	3.776794965501753244735849596741095e+179Q,
	2.473820571905781696815263976656646e+181Q,
	1.632734494231896463871010812229020e+183Q,
	1.085776900182112433178118702418846e+185Q,
	7.274761083941356085854772781448162e+186Q,
	4.910500878112028620595383876113036e+188Q,
	3.339165488138473217218874597110657e+190Q,
	2.287345162743856460926658260024207e+192Q,
	1.578279589877007679944175008061088e+194Q,
	1.096912143732327267162788655280075e+196Q,
	7.678439030562053888621945862209508e+197Q,
	5.413337067597845028126988176489658e+199Q,
	3.843495606538893036428145904978289e+201Q,
	2.748117894051497462467017590510021e+203Q,
	1.978658045201245086278275974638434e+205Q,
	1.434536494196155737613994555908787e+207Q,
	1.047218417674069520637626522764264e+209Q,
	7.697104507405240277996644921039320e+210Q,
	5.695893209382317241028283445171450e+212Q,
	4.243466810865795707145890742039051e+214Q,
	3.182619623725187774636952759853157e+216Q,
	2.402892356389598847645922346082956e+218Q,
	1.826209097230261743668270472503557e+220Q,
	1.397058194451033126500669860015331e+222Q,
	1.075741068948596811317372874066017e+224Q,
	8.337041171685535036499317363616066e+225Q,
	6.502928990423792370644965152066054e+227Q,
	5.104827839273246735782117989682823e+229Q,
	4.032836288744939315252568277004871e+231Q,
	3.206122353181952211091569786157574e+233Q,
	2.564911711575725597169542455779864e+235Q,
	2.064764922810364008523252808343359e+237Q,
	1.672468384191301162360585470279629e+239Q,
	1.363068815045044048700750297613391e+241Q,
	1.117722165158548262375897654608207e+243Q,
	9.221254621350072970000351429368717e+244Q,
	7.653679680924408494545443051776191e+246Q,
	6.390854170806011898205425708166369e+248Q,
	5.368343764383406034975374020836526e+250Q,
	4.536272410412819535036360064657603e+252Q,
	3.855849971009808287771220530105424e+254Q,
	3.296767293083101080425400381488402e+256Q,
	2.835233105703092985287715014497742e+258Q,
	2.452487952018686699407435459583229e+260Q,
	2.133674250351305121685656974581297e+262Q,
	1.866973387910968178228113069523844e+264Q,
	1.642943906272935388900446924725668e+266Q,
	1.454011766789009870266895920263165e+268Q,
	1.294076113394138357046345046315099e+270Q,
	1.158203114065351800769276063153581e+272Q,
	1.042387246347866902121909419376833e+274Q,
	9.433644353076337553141853273039645e+275Q,
	8.584652159889149460614524934193433e+277Q,
	7.854989126102919164182174880475904e+279Q,
	7.226619481709603434798217189493089e+281Q,
	6.684650001687572355815272936635608e+283Q,
	6.216749325730173453907148480777505e+285Q,
	5.812683583318519397833495247776379e+287Q,
	5.463943925916348704713771517239779e+289Q,
	5.163446980546234307930774002547136e+291Q,
	4.905293404959078737946931385804023e+293Q,
	4.684572943682552714555356607887959e+295Q,
	4.497206881767661668698799206542388e+297Q,
	4.339820739154940731305469165295673e+299Q,
	4.209641572182355026961502047135824e+301Q,
	4.104415447991076275844674535292913e+303Q,
	4.022341607500826066868516537658598e+305Q,
	3.962020590927621689424926459755537e+307Q,
	3.922414211165019563896148387670753e+309Q,
	3.902815759602983233505708951803189e+311Q,
      };

      static constexpr __float128
      s_d[s_max_cd]
      {
	-1.000000000000000000000000000000000e+00Q,
	-9.722222222222222222222222222222222e-02Q,
	-4.388503086419753086419753086419753e-02Q,
	-4.246283078989483310470964791952445e-02Q,
	-6.266216349203230579688055682568713e-02Q,
	-1.241058960272750945365995472686526e-01Q,
	-3.082537649010791121244706347668154e-01Q,
	-9.204799924129445709272387010397958e-01Q,
	-3.210493584648620907973650261091926e+00Q,
	-1.280729308073562507270352766191764e+01Q,
	-5.750830351391427202784792351524963e+01Q,
	-2.870332371092211077349530828987144e+02Q,
	-1.576357303337099717826796734206481e+03Q,
	-9.446354823095931962917203933936062e+03Q,
	-6.133570666385205823144156720993206e+04Q,
	-4.289524004000690702056279232746453e+05Q,
	-3.214536521400864829067001615998275e+06Q,
	-2.569790838391132545132402844162020e+07Q,
	-2.182934208321603255352054236989173e+08Q,
	-1.963523788991032752712502001911679e+09Q,
	-1.864393108810721585266530546676276e+10Q,
	-1.863529963852938843791870115867630e+11Q,
	-1.955882932389842694320697012392636e+12Q,
	-2.150644463519724977106616660546951e+13Q,
	-2.472369922906211612860123840379929e+14Q,
	-2.965882430295212630916036338073545e+15Q,
	-3.706244000635465228366390921824489e+16Q,
	-4.816782647945217540878439641969944e+17Q,
	-6.500984080751062701873088502894851e+18Q,
	-9.099198264365412234781657638750099e+19Q,
	-1.319088866907750709757953915010101e+21Q,
	-1.978219607616628114145519327828545e+22Q,
	-3.065639370223598386092264218755129e+23Q,
	-4.904119815775620835731518126711435e+24Q,
	-8.090395374187028082149401942289269e+25Q,
	-1.375142480406956245407560846801890e+27Q,
	-2.406127967357125254551277279514124e+28Q,
	-4.330398100410561949304091184921347e+29Q,
	-8.010128562268937488754778902693146e+30Q,
	-1.521724744139018769008631341040477e+32Q,
	-2.966993387417997254727141517133538e+33Q,
	-5.933283219493449591406075378758271e+34Q,
	-1.216185715477187410460666608307974e+36Q,
	-2.553715025111643293496042491585695e+37Q,
	-5.489923036149888462864519377823350e+38Q,
	-1.207664458504663581523897807455567e+40Q,
	-2.716989788543417797406104991755334e+41Q,
	-6.248514488575398643118502393125261e+42Q,
	-1.468274343468517211831627490866057e+44Q,
	-3.523567100049943586726891766274794e+45Q,
	-8.632054257075129854005687931752024e+46Q,
	-2.157849009857563711025991591283533e+48Q,
	-5.502111531144559820328426476011819e+49Q,
	-1.430448068378722838613634335059373e+51Q,
	-3.790429841685131698769796228639194e+52Q,
	-1.023349054707279003309533394425510e+54Q,
	-2.814032235678575023163142262900288e+55Q,
	-7.878810283641487890754406961953157e+56Q,
	-2.245328862657782165686760579825391e+58Q,
	-6.511083708721725425614962382904681e+59Q,
	-1.920664190401702861487017364214565e+61Q,
	-5.761686454417020881363820555267796e+62Q,
	-1.757224019571248449581714492600659e+64Q,
	-5.447123284124640062769737501986023e+65Q,
	-1.715761087400761462479847113120142e+67Q,
	-5.490178848750560497665481725023091e+68Q,
	-1.784227251997254438838327734091942e+70Q,
	-5.887691026309762600465986740286605e+71Q,
	-1.972292315224750519484938764884938e+73Q,
	-6.705515972283343125877688850299728e+74Q,
	-2.313309878271471665328832129897152e+76Q,
	-8.096267806165567489416614688112071e+77Q,
	-2.874065746584912340354730937461081e+79Q,
	-1.034625391639240256861069798223123e+81Q,
	-3.776246748970061121443091935275144e+82Q,
	-1.397162345772176706734221605600790e+84Q,
	-5.239180066082424250455320429788797e+85Q,
	-1.990822273847860578503192782001185e+87Q,
	-7.664417610512323501025584191661221e+88Q,
	-2.989028545098270324569269010751143e+90Q,
	-1.180629950314136764503000174380252e+92Q,
	-4.722378093272112919511460213960981e+93Q,
	-1.912507137520034823204247449964669e+95Q,
	-7.841055241911750535449288210649155e+96Q,
	-3.253947172439187678802479499916961e+98Q,
	-1.366620594074447265760845936768229e+100Q,
	-5.807983029594202540806910381110625e+101Q,
	-2.497367798700591448521132306194511e+103Q,
	-1.086327401565769110693593361391987e+105Q,
	-4.779721898165742383185161989351806e+106Q,
	-2.126924611885472705892887369005496e+108Q,
	-9.570933517949169326869927681104966e+109Q,
	-4.354673608555420286199717925981008e+111Q,
	-2.003104336167184060531793486653949e+113Q,
	-9.314227986310485516691165325316653e+114Q,
	-4.377591832519284664484296769920352e+116Q,
	-2.079311787196041530744839301827055e+118Q,
	-9.980488171002191038716190649412579e+119Q,
	-4.840437750156588425691650915454689e+121Q,
	-2.371766962413637565234560124840463e+123Q,
	-1.174001587549282684425889703572593e+125Q,
	-5.869894928792716433972256806807892e+126Q,
	-2.964240989606087663786929494358400e+128Q,
	-1.511734925078113281058779419001883e+130Q,
	-7.785293542778411244768854346322128e+131Q,
	-4.048280556193451824953463807732449e+133Q,
	-2.125310161545727100444941305340125e+135Q,
	-1.126395074649439387003078608095664e+137Q,
	-6.026112250640928883073093218560559e+138Q,
	-3.254046865619095370429575445499091e+140Q,
	-1.773426781247180545485642598822107e+142Q,
	-9.753691966685960339212345611111786e+143Q,
	-5.413214374045719394857129484778055e+145Q,
	-3.031353475595619675235326828208266e+147Q,
	-1.712688861525451567119891130540256e+149Q,
	-9.762181718158471539334951433959726e+150Q,
	-5.613172668889364538732569123544997e+152Q,
	-3.255593504783131734516443833365660e+154Q,
	-1.904495376905319585970201545641753e+156Q,
	-1.123636712771386061281357063876791e+158Q,
	-6.685547405607925308808924420119601e+159Q,
	-4.011274725697354754296877078069304e+161Q,
	-2.426789243059786828668280852929166e+163Q,
	-1.480322256328008562475229597199652e+165Q,
	-9.103865811724216189017270819740446e+166Q,
	-5.644325995423811737731459996132700e+168Q,
	-3.527660195241769498028630619033585e+170Q,
	-2.222398917729604052801687552005230e+172Q,
	-1.411206432558185609956621929927765e+174Q,
	-9.031614811298444026684592923097672e+175Q,
	-5.825324009159634122523328683100498e+177Q,
	-3.786417373057171724340272143484384e+179Q,
	-2.480075491177351309196592077735423e+181Q,
	-1.636831694970245012086496410252179e+183Q,
	-1.088481201303363024543444029199466e+185Q,
	-7.292745660168652392618319809338021e+186Q,
	-4.922551187015371022118317333477730e+188Q,
	-3.347299874224072421158506447529927e+190Q,
	-2.292876831819415968691898062345910e+192Q,
	-1.582068976647420663569491154539025e+194Q,
	-1.099526952179841754092854897604938e+196Q,
	-7.696612850752733306109808455966215e+197Q,
	-5.426059363878921044644325398996096e+199Q,
	-3.852465257896043311892388952597841e+201Q,
	-2.754486649310017734685944630117228e+203Q,
	-1.983211918722997088778341051680179e+205Q,
	-1.437815434754318379299969400608007e+207Q,
	-1.049595758009311449174828853122412e+209Q,
	-7.714459872698149500720425405641437e+210Q,
	-5.708649969089780437536745446168474e+212Q,
	-4.252907226462827510721298730341696e+214Q,
	-3.189653037258282112260459837775483e+216Q,
	-2.408167641474976671680271242561732e+218Q,
	-1.830192105075911169499608030786007e+220Q,
	-1.400085406139984444217897747035949e+222Q,
	-1.078056980830079258704493160124286e+224Q,
	-8.354874414833525485775251732308293e+225Q,
	-6.516750306025118178021468797447702e+227Q,
	-5.115608890676146940081552240980990e+229Q,
	-4.041299743705579271842815010010127e+231Q,
	-3.212808739737076199018768054741844e+233Q,
	-2.570227590770701194262121818382516e+235Q,
	-2.069017785679180412248326449555188e+237Q,
	-1.675892065632839138127956167250518e+239Q,
	-1.365842098493762347884271661392869e+241Q,
	-1.119982472873732384241167417307111e+243Q,
	-9.239789806518615830241558165910161e+244Q,
	-7.668971748218962757271807573358159e+246Q,
	-6.403547029139290968509706593386163e+248Q,
	-5.378942666188703973840083545063250e+250Q,
	-4.545175791002442340796980987257519e+252Q,
	-3.863373580709339621015652179920263e+254Q,
	-3.303162573962020772143005425875382e+256Q,
	-2.840701250554979374844682642298123e+258Q,
	-2.457190709357169320115791040522027e+260Q,
	-2.137742266081240879782293117526161e+262Q,
	-1.870512673954401293257929397617728e+264Q,
	-1.646040878763553551744745599418836e+266Q,
	-1.456737187157873993735062548042477e+268Q,
	-1.296488184434947561812507851620440e+270Q,
	-1.160349922432479422272092144827638e+272Q,
	-1.044308697493208592264069621071536e+274Q,
	-9.450937926592517823633405708003970e+275Q,
	-8.600303303298245448692217778193785e+277Q,
	-7.869232080094039597843429957321735e+279Q,
	-7.239652158863272692570621548716701e+281Q,
	-6.696640405278043337619425892575762e+283Q,
	-6.227840760744857081835617969592452e+285Q,
	-5.822998904673121916729384325412184e+287Q,
	-5.473589016694665295542922040659443e+289Q,
	-5.172513612645525325152777117564777e+291Q,
	-4.913861603046343504301423842373113e+293Q,
	-4.692712948797552806153193891307399e+295Q,
	-4.504980791675643797771744408282305e+297Q,
	-4.347283887459592392064377968675372e+299Q,
	-4.216843696343488226323283915479940e+301Q,
	-4.111401687051486618441856960033836e+303Q,
	-4.029153362975001894246786675740999e+305Q,
	-3.968696278528173705750831980327997e+307Q,
	-3.928989926523217417314247546744802e+309Q,
	-3.909325877634014064587453253641059e+311Q,
      };
    };
#endif // EMSR_HAVE_FLOAT128


  /**
   * A class encapsulating the asymptotic expansions of Airy functions
   * and their derivatives.
   *
   * @tparam Tp A real type
   */
  template<typename Tp>
    class Airy_asymp : public Airy_asymp_data<Tp>
    {

    public:

      using Cmplx = std::complex<Tp>;

      constexpr Airy_asymp() = default;

      AiryState<Cmplx>
      operator()(Cmplx t, bool return_fock_airy = false) const;

      AiryState<Cmplx>
      s_absarg_ge_pio3(Cmplx z) const;

      AiryState<Cmplx>
      s_absarg_lt_pio3(Cmplx z) const;

    private:
      std::pair<Cmplx, Cmplx>
      s_absarg_ge_pio3_help(Cmplx z, int sign = -1) const;
    };

  /**
   * Return the Airy functions for a given argument using asymptotic series.
   *
   *
   * @tparam Tp A real type
   */
  template<typename Tp>
    AiryState<std::complex<Tp>>
    Airy_asymp<Tp>::operator()(std::complex<Tp> t,
				 bool return_fock_airy) const
    {
      const auto s_pi = emsr::pi_v<Tp>;
      const auto s_sqrt_pi = emsr::sqrtpi_v<Tp>;
      const auto s_i = Cmplx(Tp{0}, Tp{1});
      if (std::real(t) > Tp{0})
	{
	  auto zeta0 = Tp{2} * std::pow(t, Tp{1.5L}) / Tp{3};
	  auto t0p1d4 = std::pow(t, Tp{0.25L});
	  auto ezeta0 = std::exp(-zeta0);
	  auto Ai = Cmplx{Tp{1}};
	  auto Aip = Cmplx{Tp{1}};
	  auto fact0 = -Tp{1} / zeta0;
	  auto izeta0 = Cmplx{Tp{1}};
	  auto prev_Ai0 = Tp{1};
	  auto prev_Aip0 = Tp{1};
	  for (int n = 1; n < Airy_asymp_data<Tp>::s_max_cd; ++n)
	    {
	      izeta0 *= fact0;
	      auto term = Airy_asymp_data<Tp>::s_c[n] * izeta0;
	      auto termp = Airy_asymp_data<Tp>::s_d[n] * izeta0;
	      if (std::abs(term) > prev_Ai0
	       || std::abs(termp) > prev_Aip0)
		break;
	      prev_Ai0 = std::abs(term);
	      prev_Aip0 = std::abs(termp);
	      Ai += term;
	      Aip += termp;
	    }
	  Ai *= Tp{0.5L} * ezeta0 / t0p1d4 / s_sqrt_pi;
	  Aip *= Tp{-0.5L} * t0p1d4 * ezeta0 / s_sqrt_pi;

	  auto t1 = t * std::exp(+Tp{2} * s_pi * s_i / Tp{3});
	  auto t2 = t * std::exp(-Tp{2} * s_pi * s_i / Tp{3});
	  auto zeta1 = (Tp{2} / Tp{3}) * std::pow(t1, Tp{1.5L});
	  auto zeta2 = (Tp{2} / Tp{3}) * std::pow(t2, Tp{1.5L});
	  auto t1p1d4 = std::pow(t1, Tp{+0.25L});
	  auto t2p1d4 = std::pow(t2, Tp{+0.25L});
	  auto ezeta1 = std::exp(-zeta1);
	  auto ezeta2 = std::exp(-zeta2);
	  auto Ai1 = Cmplx{Tp{1}};
	  auto Ai1p = Cmplx{Tp{1}};
	  auto Ai2 = Ai1;
	  auto Ai2p = Ai1p;
	  auto sign = Tp{1};
	  auto izeta1 = Cmplx{Tp{1}};
	  auto izeta2 = Cmplx{Tp{1}};
	  auto prev_Ai1 = Tp{1};
	  auto prev_Ai2 = Tp{1};
	  auto prev_Ai1p = Tp{1};
	  auto prev_Ai2p = Tp{1};
	  for (int n = 1; n < Airy_asymp_data<Tp>::s_max_cd; ++n)
	    {
	      sign = -sign;
	      izeta1 /= zeta1;
	      izeta2 /= zeta2;
	      const auto term1 = Airy_asymp_data<Tp>::s_c[n] * izeta1;
	      const auto term2 = Airy_asymp_data<Tp>::s_c[n] * izeta2;
	      const auto term1p = Airy_asymp_data<Tp>::s_d[n] * izeta1;
	      const auto term2p = Airy_asymp_data<Tp>::s_d[n] * izeta2;
	      if (std::abs(term1) > prev_Ai1
		  || std::abs(term2) > prev_Ai2
		  || std::abs(term1p) > prev_Ai1p
		  || std::abs(term2p) > prev_Ai2p)
		break;
	      prev_Ai1 = std::abs(term1);
	      prev_Ai2 = std::abs(term2);
	      prev_Ai1p = std::abs(term1p);
	      prev_Ai2p = std::abs(term2p);
	      Ai1 += sign * term1;
	      Ai2 += sign * term2;
	      Ai1p += sign * term1p;
	      Ai2p += sign * term2p;
	    }
	  Ai1 *= Tp{+0.5L} * ezeta1 / t1p1d4 / s_sqrt_pi;
	  Ai2 *= Tp{+0.5L} * ezeta2 / t2p1d4 / s_sqrt_pi;
	  Ai1p *= Tp{-0.5L} * t1p1d4 * ezeta1 / s_sqrt_pi;
	  Ai2p *= Tp{-0.5L} * t2p1d4 * ezeta2 / s_sqrt_pi;

	  auto Bi = std::exp(+s_i * s_pi / Tp{6}) * Ai1
		   + std::exp(-s_i * s_pi / Tp{6}) * Ai2;
	  auto Bip = std::exp(+s_i * Tp{5} * s_pi / Tp{6}) * Ai1p
		    + std::exp(-s_i * Tp{5} * s_pi / Tp{6}) * Ai2p;

	  if (return_fock_airy)
	    {
	      auto w1 = s_sqrt_pi * (Bi - s_i * Ai);
	      auto w2 = s_sqrt_pi * (Bi + s_i * Ai);
	      auto w1p = s_sqrt_pi * (Bip - s_i * Aip);
	      auto w2p = s_sqrt_pi * (Bip + s_i * Aip);
	      return AiryState<Cmplx>{t, w1, w1p, w2, w2p};
	    }
	  else
	    return AiryState<Cmplx>{t, Ai, Aip, Bi, Bip};
	}
      else // Argument t is on or left of the imaginary axis.
	{
	  auto zeta = (Tp{2} / Tp{3}) * std::pow(-t, Tp{1.5L});
	  auto tp1d4 = std::pow(-t, Tp{+0.25L});
	  auto mezeta = std::exp(-s_i * (zeta + (s_pi / Tp{4})));
	  auto pezeta = std::exp(+s_i * (zeta + (s_pi / Tp{4})));
	  auto w1 = Cmplx{Tp{1}};
	  auto w2 = Cmplx{Tp{1}};
	  auto w1p = +s_i;
	  auto w2p = -s_i;
	  auto ipn = Cmplx{Tp{1}};
	  auto imn = Cmplx{Tp{1}};
	  auto ixn = Cmplx{Tp{1}};
	  auto prev_w1 = Tp{1};
	  auto prev_w2 = Tp{1};
	  auto prev_w1p = Tp{1};
	  auto prev_w2p = Tp{1};
	  for (int n = 1; n < Airy_asymp_data<Tp>::s_max_cd; ++n)
	    {
	      ipn *= +s_i;
	      imn *= -s_i;
	      ixn /= zeta;
	      const auto term = Airy_asymp_data<Tp>::s_c[n] * ixn;
	      const auto termp = Airy_asymp_data<Tp>::s_d[n] * ixn;
	      if (std::abs(term) > prev_w1
	       || std::abs(term) > prev_w2
	       || std::abs(termp) > prev_w1p
	       || std::abs(termp) > prev_w2p)
		break;
	      prev_w1 = std::abs(term);
	      prev_w2 = std::abs(term);
	      prev_w1p = std::abs(termp);
	      prev_w2p = std::abs(termp);
	      w1 += ipn * term;
	      w2 += imn * term;
	      w1p += +s_i * ipn * termp;
	      w2p += -s_i * imn * termp;
	    }
	  w1 *= mezeta / tp1d4;
	  w2 *= pezeta / tp1d4;
	  w1p *= tp1d4 * mezeta;
	  w2p *= tp1d4 * pezeta;

	  if (return_fock_airy)
	    return AiryState<Cmplx>{t, w1, w1p, w2, w2p};
	  else
	    {
	      auto Bi = (w1 + w2) / (Tp{2} * s_sqrt_pi);
	      auto Ai = (w2 - w1) / (Tp{2} * s_i * s_sqrt_pi);
	      auto Bip = (w1p + w2p) / (Tp{2} * s_sqrt_pi);
	      auto Aip = (w2p - w1p) / (Tp{2} * s_i * s_sqrt_pi);
	      return AiryState<Cmplx>{t, Ai, Aip, Bi, Bip};
	    }
	}
    }

  /**
   * @brief This function evaluates @f$ Ai(z) @f$ and @f$ Ai'(z) @f$
   * or @f$ Bi(z) @f$ and @f$ Bi'(z) @f$ from their asymptotic expansions
   * for @f$ |arg(z)| < 2*\pi/3 @f$ i.e. roughly along the negative real axis.
   *
   * For speed, the number of terms needed to achieve about 16 decimals accuracy
   * is tabled and determined from @f$ |z| @f$.
   *
   * Note that for speed and since this function
   * is called by another, checks for valid arguments are not
   * made.
   *
   * @see Digital Library of Mathematical Functions §9.7 Asymptotic Expansions
   * 	  http://dlmf.nist.gov/9.7
   *
   * @tparam Tp A real type
   * @param[in]  z Complex arument at which @f$ Ai(z) @f$ or @f$ Bi(z) @f$
   * 		   and their derivative are evaluated. This function assumes
   * 		   @f$ |z| > 15 @f$ and @f$ |arg(z)| < 2\pi/3 @f$.
   * @param[inout] Ai  The value computed for @f$ Ai(z) @f$ or @f$ Bi(x) @f$.
   * @param[inout] Aip The value computed for @f$ Ai'(z) @f$ or @f$ Bi'(x) @f$.
   * @param[in]    sign  The sign of the series terms and exponent.
   * 			 The default (-1) gives the Airy @f$ Ai(x) @f$ and
   * 			 @f$ Ai'(x) @f$ functions for @f$ |arg(z)| < \pi @f$.
   * 			 The value +1 gives the Airy @f$ Bi(x) @f$ and
   * 			 @f$ Bi'(x) @f$ functions for @f$ |arg(z)| < \pi/3 @f$.
   * @return A pair containing the Airy function @f$ Ai(z) @f$ or @f$ Bi(x) @f$
   *         and the derivative @f$ Ai'(z) @f$ or @f$ Bi'(x) @f$
   */
  template<typename Tp>
    std::pair<std::complex<Tp>, std::complex<Tp>>
    Airy_asymp<Tp>::s_absarg_ge_pio3_help(std::complex<Tp> z,
					     int sign) const
    {
      using Real = emsr::num_traits_t<Tp>;
      const auto s_sqrt_pi = emsr::sqrtpi_v<Real>;
      const auto s_pmhd2 = Tp{1} / (Tp{2} * s_sqrt_pi);

      constexpr int s_num_nterms = 5;
      constexpr int s_max_nterms = 40;
      static_assert(Airy_asymp_data<Tp>::s_max_cd > s_max_nterms, "");
      constexpr int s_nterms[s_num_nterms]{s_max_nterms, 24, 22, 22, 18};

      auto zeta = Tp{2} * std::pow(z, Tp{1.5L}) / Tp{3};
      auto z1d4 = std::pow(z, Tp{0.25L});

      // Compute outer factors in the expansions.
      auto exp = std::exp(Tp(sign) * zeta);
      auto fact = s_pmhd2 * exp / z1d4;
      auto factp = s_pmhd2 * exp * z1d4;
      if (sign == +1)
	{
	  fact *= Tp{2};
	  factp *= -Tp{2};
	}

      // Determine number of terms to use.
      auto iterm = std::min(s_num_nterms - 1, (int(std::abs(z)) - 10) / 5);
      if (iterm < 0 || iterm >= s_num_nterms)
	iterm = 0;
      auto nterm = s_nterms[iterm];
      // Power series is in terms of +-1 / \zeta.
      auto zetam = Tp(sign) / zeta;

      emsr::Polynomial<Tp>
	cpoly(std::begin(Airy_asymp_data<Tp>::s_c),
		std::begin(Airy_asymp_data<Tp>::s_c) + nterm);
      auto Ai = fact * cpoly(zetam);

      emsr::Polynomial<Tp>
	dpoly(std::begin(Airy_asymp_data<Tp>::s_d),
		std::begin(Airy_asymp_data<Tp>::s_d) + nterm);
      auto Aip = factp * dpoly(zetam);

      return std::make_pair(Ai, Aip);
    }


  /**
   * @brief This function evaluates @f$ Ai(z), Ai'(z) @f$
   * and @f$ Bi(z), Bi'(z) @f$ from their asymptotic expansions
   * for @f$ |arg(z)| < 2*\pi/3 @f$ i.e. roughly along the negative real axis.
   *
   * @tparam Tp A real type
   * @param[in]  z Complex argument at which Ai(z) and Bi(z)
   * 		   and their derivative are evaluated. This function
   * 		   assumes @f$ |z| > 15 @f$ and @f$ |(arg(z)| < 2\pi/3 @f$.
   * @return A struct containing @f$ z, Ai(z), Ai'(z), Bi(z), Bi'(z) @f$.
   */
  template<typename Tp>
    AiryState<std::complex<Tp>>
    Airy_asymp<Tp>::s_absarg_ge_pio3(std::complex<Tp> z) const
    {
      Cmplx Ai, Aip;
      std::tie(Ai, Aip) = s_absarg_ge_pio3_help(z, -1);
      Cmplx Bi, Bip;
      std::tie(Bi, Bip) = s_absarg_ge_pio3_help(z, +1);
      return AiryState<Cmplx>{z, Ai, Aip, Bi, Bip};
    }


  /**
   * @brief This function evaluates @f$ Ai(z) @f$ and @f$ Ai'(z) @f$
   * from their asymptotic expansions for @f$ |arg(-z)| < \pi/3 @f$
   * i.e. roughly along the negative real axis.
   *
   * For speed, the number of terms needed to achieve about 16 decimals
   * accuracy is tabled and determined for @f$ |z| @f$.
   * This function assumes @f$ |z| > 15 @f$ and @f$ |arg(-z)| < \pi/3 @f$.
   *
   * Note that for speed and since this function
   * is called by another, checks for valid arguments are not
   * made.  Hence, an error return is not needed.
   *
   * @tparam Tp A real type
   * @param[in] z  The value at which the Airy function and their derivatives
   * 		   are evaluated.
   * @return A struct containing @f$ z, Ai(z), Ai'(z), Bi(z), Bi'(z) @f$.
   */
  template<typename Tp>
    AiryState<std::complex<Tp>>
    Airy_asymp<Tp>::s_absarg_lt_pio3(std::complex<Tp> z) const
    {
      const Tp s_pimh = Tp{1} / emsr::sqrtpi_v<Tp>;
      const Tp s_pid4 = emsr::pi_v<Tp> / Tp{4};

      const Cmplx s_zone{1};

      /// @todo Revisit these numbers of terms for the Airy asymptotic
      /// expansions.
      constexpr int s_num_nterms = 5;
      constexpr int s_max_nterms = 40;
      static_assert(Airy_asymp_data<Tp>::s_max_cd > s_max_nterms, "");
      constexpr int s_nterms[s_num_nterms]{s_max_nterms, 28, 24, 24, 20};

      auto zeta = Tp{2} * std::pow(-z, Tp{1.5L}) / Tp{3};
      auto z1d4 = std::pow(-z, Tp{0.25L});

      auto zetaarg = zeta - s_pid4;
      auto sinzeta = std::sin(zetaarg);
      auto coszeta = std::cos(zetaarg);

      // Determine number of terms to use.
      auto iterm = std::min(s_num_nterms - 1, (int(std::abs(z)) - 10) / 5);
      if (iterm < 0 || iterm >= s_num_nterms)
	iterm = 0;
      auto nterm = s_nterms[iterm];
      // Power series is in terms of 1 / \zeta^2.
      auto zetam2 = Tp{1} / (zeta * zeta);

      emsr::Polynomial<Tp>
	cpoly(std::begin(Airy_asymp_data<Tp>::s_c),
		std::begin(Airy_asymp_data<Tp>::s_c) + nterm);

      emsr::Polynomial<Tp>
	dpoly(std::begin(Airy_asymp_data<Tp>::s_d),
		std::begin(Airy_asymp_data<Tp>::s_d) + nterm);

      // Complete evaluation of the Airy functions.
      zeta = s_zone / zeta;
      auto Ai = coszeta * cpoly.even(zetam2)
	       + sinzeta * cpoly.odd(zetam2);
      Ai *= s_pimh / z1d4;
      auto Aip = sinzeta * dpoly.even(zetam2)
		- coszeta * dpoly.odd(zetam2);
      Aip *= s_pimh * z1d4;
      auto Bi = -sinzeta * cpoly.even(zetam2)
	       + coszeta * cpoly.odd(zetam2);
      Bi *= s_pimh / z1d4;
      auto Bip = coszeta * dpoly.even(zetam2)
		+ sinzeta * dpoly.odd(zetam2);
      Bip *= s_pimh * z1d4;

      // I think we're computing d/d(-z) above.
      return AiryState<Cmplx>{z, Ai, -Aip, Bi, -Bip};
    }


  /**
   * Class to manage the asymptotic series for Airy functions.
   *
   * @tparam Sum A sum type
   */
  template<typename Sum>
    class Airy_asymp_series
    {
    public:

      using value_type = typename Sum::value_type;
      using scalar_type = emsr::num_traits_t<value_type>;
      static constexpr scalar_type s_sqrt_pi = emsr::sqrtpi_v<scalar_type>;

      Airy_asymp_series(Sum proto)
      : m_Asum(proto),
	m_Bsum(proto),
	m_Csum(proto),
	m_Dsum(proto)
      { }
      Airy_asymp_series(const Airy_asymp_series&) = default;
      Airy_asymp_series(Airy_asymp_series&&) = default;

      AiryState<value_type>
      operator()(value_type y);

    private:

      static constexpr int s_max_iter = 10000;
      static constexpr scalar_type s_eps
	   = std::numeric_limits<scalar_type>::epsilon();

      Sum m_Asum;
      Sum m_Bsum;
      Sum m_Csum;
      Sum m_Dsum;
    };

  template<typename Sum>
    constexpr int
    Airy_asymp_series<Sum>::s_max_iter;

  template<typename Sum>
    constexpr typename Airy_asymp_series<Sum>::scalar_type
    Airy_asymp_series<Sum>::s_eps;

  template<typename Sum>
    constexpr typename Airy_asymp_series<Sum>::scalar_type
    Airy_asymp_series<Sum>::s_sqrt_pi;


  /**
   * Return an AiryState containing, not actual Airy functions, but
   * four asymptotic Airy components:
   *
   * @tparam Sum A sum type
   */
  template<typename Sum>
    AiryState<typename Airy_asymp_series<Sum>::value_type>
    Airy_asymp_series<Sum>::operator()(typename Sum::value_type y)
    {
      using Cmplx = value_type;
      using Real = scalar_type;

      m_Asum.reset(Real{1});
      m_Bsum.reset(Real{1});
      m_Csum.reset(Real{1});
      m_Dsum.reset(Real{1});

      auto zeta = Real{2} * std::pow(y, Real{1.5L}) / Real{3};
      auto sign = Real{1};
      auto numerAB = Cmplx{1};
      auto numerCD = Cmplx{1};
      for (int k = 1; k < s_max_iter; ++k)
	{
	  sign = -sign;
	  auto denom = Cmplx(2 * k) * zeta;
	  numerAB *= Real(k + Real{1} / Real{6})
		     * Real(k + Real{5} / Real{6})
		     / denom;
	  if (k > 1 && (std::abs(m_Asum.term()) < std::abs(numerAB)
		|| std::isinf(numerAB)))
	    break;
	  auto Aterm = sign * numerAB;
	  m_Asum += Aterm;
	  auto Bterm = numerAB;
	  m_Bsum += Bterm;
	  numerCD *= Real(k - Real{1} / Real{6})
		     * Real(k + Real{7} / Real{6})
		     / denom;
	  if (k > 1 && (std::abs(m_Csum.term()) < std::abs(numerCD)
		|| std::isinf(numerCD)))
	    break;
	  auto Cterm = sign * numerCD;
	  m_Csum += Cterm;
	  auto Dterm = numerCD;
	  m_Dsum += Dterm;
	  if (std::abs(Aterm) < std::abs(m_Asum()) * s_eps
	   && std::abs(Bterm) < std::abs(m_Bsum()) * s_eps
	   && std::abs(Cterm) < std::abs(m_Csum()) * s_eps
	   && std::abs(Dterm) < std::abs(m_Dsum()) * s_eps)
	    break;
	}

      auto expzeta = std::exp(zeta);
      auto y1o4 = std::pow(y, Real{0.25L});
      auto AA = Real{0.5L} * m_Asum() / s_sqrt_pi / y1o4 / expzeta;
      auto BB = Real{0.5L} * expzeta * m_Bsum() / s_sqrt_pi / y1o4;
      auto CC = Real{-0.5L} * y1o4 * m_Csum() / s_sqrt_pi / expzeta;
      auto DD = Real{0.5L} * y1o4 * expzeta * m_Dsum() / s_sqrt_pi;

      return AiryState<value_type>{y, AA, CC, BB, DD};
    }


  template<typename Tp>
    struct Airy_default_radii
    {};

  template<>
    struct Airy_default_radii<float>
    {
      static constexpr float inner_radius{2.0F};
      static constexpr float outer_radius{6.0F};
    };

  template<>
    struct Airy_default_radii<double>
    {
      static constexpr double inner_radius{4.0};
      static constexpr double outer_radius{12.0};
    };

  template<>
    struct Airy_default_radii<long double>
    {
      static constexpr long double inner_radius{5.0L};
      static constexpr long double outer_radius{15.0L};
    };

  /**
   * Class to manage the asymptotic expansions for Airy functions.
   * The parameters describing the various regions are adjustable.
   */
  template<typename Tp>
    class Airy
    {
    public:

      using value_type = Tp;
      using scalar_type = emsr::num_traits_t<value_type>;

      constexpr Airy() = default;
      Airy(const Airy&) = default;
      Airy(Airy&&) = default;

      constexpr AiryState<value_type>
      operator()(value_type y) const;

      scalar_type inner_radius{Airy_default_radii<scalar_type>::inner_radius};
      scalar_type outer_radius{Airy_default_radii<scalar_type>::outer_radius};
    };

  /**
   * Return the Airy functions for complex argument.
   */
  template<typename Tp>
    constexpr AiryState<Tp>
    Airy<Tp>::operator()(typename Airy<Tp>::value_type y) const
    {
      using Cmplx = value_type;
      using Real = scalar_type;
      const auto s_NaN = emsr::quiet_NaN(y.real());
      const auto s_cNaN = value_type(s_NaN, s_NaN);
      const auto s_pi = emsr::pi_v<Real>;
      const auto s_pi_3 = emsr::pi_v<Real> / Real{3};
      const auto s_2pi_3 = Real{2} * s_pi_3;
      const auto s_pi_6 = s_pi_3 / Real{2};
      const auto s_5pi_6 = Real{5} * s_pi_6;
      const auto s_i = Cmplx{0, 1};

      using OuterSum = emsr::KahanSum<Cmplx>;
      using InnerSum = emsr::WenigerDeltaSum<OuterSum>;

      if (std::isnan(y))
	return AiryState<Tp>{y, s_cNaN, s_cNaN, s_cNaN, s_cNaN};

      auto absargy = std::abs(std::arg(y));
      auto absy = std::abs(y);
      auto sign = std::copysign(Real{1}, std::arg(y));

      AiryState<Tp> sums;
      if (absy >= inner_radius)
	{
	  if (absy < outer_radius)
	    {
	      auto beta = Real{1};
	      Airy_asymp_series<InnerSum> asymp(InnerSum{beta});
	      sums = asymp(y);
	    }
	  else
	    {
	      Airy_asymp_series<OuterSum> asymp(OuterSum{});
	      sums = asymp(y);
	    }
	}

      Cmplx Bi, Bip;
      if (absy < inner_radius
	  || (absy < outer_radius && absargy < s_pi_3))
	std::tie(Bi, Bip) = AirySeries<Real>::s_Bi(y);
      else if (absy < outer_radius)
	{
	  Bi = Real{2} * sums.Bi_value
	      + sign * s_i * sums.Ai_value;
	  Bip = Real{2} * sums.Bi_deriv
	       + sign * s_i * sums.Ai_deriv;
	  if (absargy > s_5pi_6)
	    {
	      Bi -= sums.Bi_value;
	      Bip -= sums.Bi_deriv;
	    }
	}
      else
	{
	  Bi = Real{2} * sums.Bi_value;
	  Bip = Real{2} * sums.Bi_deriv;
	  if (absargy > s_pi_6)
	    {
	      Bi += sign * s_i * sums.Ai_value;
	      Bip += sign * s_i * sums.Ai_deriv;
	    }
	  if (absargy > s_5pi_6)
	    {
	      Bi -= sums.Bi_value;
	      Bip -= sums.Bi_deriv;
	    }
	}

      Cmplx Ai, Aip;
      if ((absy < inner_radius
	          + outer_radius * absargy / s_pi && absargy < s_2pi_3)
	  || (absy < outer_radius && absargy >= s_2pi_3))
	std::tie(Ai, Aip) = AirySeries<Real>::s_Ai(y);
      else if (absy < outer_radius)
	{
	  Ai = sums.Ai_value;
	  Aip = sums.Ai_deriv;
	}
      else
	{
	  Ai = sums.Ai_value;
	  Aip = sums.Ai_deriv;
	  if (absargy >= s_5pi_6)
	    {
	      Ai += sign * s_i * sums.Bi_value;
	      Aip += sign * s_i * sums.Bi_deriv;
	    }
	}

      return AiryState<Tp>{y, Ai, Aip, Bi, Bip};
    }


  /**
   * @brief  Return the complex Airy Ai function.
   */
  template<typename Tp>
    std::complex<Tp>
    airy_ai(const std::complex<Tp>& z)
    {
      auto airy = Airy<std::complex<Tp>>()(z);
      return airy.Ai_value;
    }


  /**
   * @brief  Return the complex Airy Bi function.
   */
  template<typename Tp>
    std::complex<Tp>
    airy_bi(const std::complex<Tp>& z)
    {
      auto airy = Airy<std::complex<Tp>>()(z);
      return airy.Bi_value;
    }

} // namespace detail
} // namespace emsr

#endif // SF_AIRY_TCC
