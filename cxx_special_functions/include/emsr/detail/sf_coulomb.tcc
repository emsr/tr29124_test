
// Copyright (C) 2016-2019 Free Software Foundation, Inc.
//
// This file is part of the GNU ISO C++ Library.  This library is free
// software; you can redistribute it and/or modify it under the
// terms of the GNU General Public License as published by the
// Free Software Foundation; either version 3, or (at your option)
// any later version.
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

#include <emsr/sf_gamma.h>

namespace emsr
{
namespace detail
{

/**
 * 
 */
template<typename Tp>
  Tp
  coulomb_norm(unsigned int l, Tp eta)
  {
    const auto s_2pi = emsr::tau_v<Tp>;
    auto _Ck = std::sqrt(s_2pi * eta / (std::exp(s_2pi * eta) - Tp{1}));
    if (l == 0)
      return _Ck;
    else
      {
	for (int k = 0; k < l; ++k)
	  _Ck *= std::hypot(Tp(k + 1), eta) / Tp(k) / Tp(k + 1);
	return _Ck;
      }
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
   *    S_l = l / x + \eta / l
   * @f]
   */
  template<typename Tp>
    std::pair<Tp, Tp>
    coulomb_f_recur(unsigned int l_min, unsigned int k_max,
		      Tp eta, Tp x,
		      Tp _F_l_max, Tp _Fp_l_max)
    {
      const auto x_inv = Tp{1}/x;
      auto fcl = _F_l_max;
      auto fpl = _Fp_l_max;
      const auto l_max = l_min + k_max;
      auto l = l_max;

      for (int k = k_max - 1; k >= 0; --k)
	{
	  const auto el = eta / l;
	  const auto rl = std::hypot(Tp{1}, el);
	  const auto sl = el  + l * x_inv;
	  const auto fc_lm1 = (fcl * sl + fpl) / rl;
	  fpl = fc_lm1 * sl - fcl * rl;
	  fcl = fc_lm1;
	  l -= 1;
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
   *    S_l = l / x + \eta / l
   * @f]
   */
  template<typename Tp>
    std::pair<Tp, Tp>
    coulomb_g_recur(unsigned int l_min, unsigned int k_max,
		      Tp eta, Tp x,
		      Tp _G_l_min, Tp _Gp_l_min)
    {
      const auto x_inv = Tp{1} / x;
      auto gcl = _G_l_min;
      auto gpl = _Gp_l_min;
      auto l = l_min + 1;

      for (int k = 1; k <= k_max; ++k)
	{
	  const auto el = eta / l;
	  const auto rl = std::hypot(Tp{1}, el);
	  const auto sl = el + l * x_inv;
	  const auto gc_lm1 = (sl * gcl - gpl) / rl;
	  gpl = rl * gcl - sl * gc_lm1;
	  gcl = gc_lm1;
	  l += 1;
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
    coulomb_CF1(unsigned int l,
        	  Tp eta, Tp x)
    {
      const auto _CF1_small = 1.e-30;
      const auto _CF1_abort = 1.0e+05;
      const auto _CF1_acc = Tp{2} * std::numeric_limits<Tp>::epsilon();
      const auto x_inv = Tp{1} / x;
      const auto px = l + Tp{1} + _CF1_abort;

      auto pk = l + 1;
      auto F = eta / pk + pk * x_inv;

      auto fcl_sign = Tp{1};

      if (std::abs(F) < _CF1_small)
	F = _CF1_small;
      auto D = Tp{0};
      auto C = F;

      Tp df;
      do
	{
	  const auto pk1 = pk + 1;
	  const auto ek  = eta / pk;
	  const auto rk2 = Tp{1} + ek * ek;
	  const auto tk  = (pk + pk1) * (x_inv + ek / pk1);
	  D = tk - rk2 * D;
	  C = tk - rk2 / C;
	  if (std::abs(C) < _CF1_small)
	    C = _CF1_small;
	  if (std::abs(D) < _CF1_small)
	    D = _CF1_small;
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
      while (std::abs(df - Tp{1}) > _CF1_acc);

      return std::make_pair(F, fcl_sign);
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
    coulomb_CF2(unsigned int l, Tp eta, Tp x)
    {
      const auto s_i = std::complex<Tp>{0, 1};
      const auto _CF2_acc = Tp{4} * std::numeric_limits<Tp>::epsilon();
      const auto _CF2_abort = 2.0e+05;

      const auto wi = Tp{2} * eta;
      const auto x_inv = Tp{1} / x;
      const auto e2mm1 = eta * eta + Tp(l * (l + 1));

      auto a = std::complex<Tp>(-e2mm1, eta);
      auto b = Tp{2} * std::complex<Tp>(x - eta, Tp{2});
      auto d = std::conj(b) / std::norm(b);

      auto dpq = x_inv * std::conj(a * d);

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
	  if (pk > _CF2_abort)
	    throw std::runtime_error("coulomb_CF2: Too many iterations.");
        }
      while (std::abs(dpq) > (std::abs(PQ)) * _CF2_acc);

      //if (Q < CF2_abort * std::numeric_limits<Tp>::epsilon() * std::abs(P))
	//status = GSL_ELOSS;

      return PQ;
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
	  const auto lag = assoc_laguerre(n - l - 1, 2 * l + 1,
                                        	rho);
	  const auto sphh = sph_legendre(l, m, theta)
 			    * std::polar(Tp(1), Tp(m) * phi);

	  const auto psi = norm * ea * pp * lag * sphh;

	  return psi;
	}
    }

} // namespace detail
} // namespace emsr

#endif // SF_COULOMB_TCC
