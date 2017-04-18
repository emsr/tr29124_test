// Special functions -*- C++ -*-

// Copyright (C) 2016-2017 Free Software Foundation, Inc.
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

/** @file bits/sf_hydrogen.tcc
 *  This is an internal header file, included by other library headers.
 *  You should not attempt to use it directly.
 */

#ifndef _GLIBCXX_BITS_SF_HYDROGEN_TCC
#define _GLIBCXX_BITS_SF_HYDROGEN_TCC 1

#pragma GCC system_header

#include <complex>

namespace std _GLIBCXX_VISIBILITY(default)
{
// Implementation-space details.
namespace __detail
{
_GLIBCXX_BEGIN_NAMESPACE_VERSION

template<typename _Tp>
  _Tp
  __coulomb_norm(unsigned int __l, _Tp __eta)
  {
    const auto _S_2pi = __gnu_cxx::__const_2_pi(__eta);
    auto _Ck = std::sqrt(_S_2pi * __eta / (std::exp(_S_2pi * __eta) - _Tp{1}));
    if (__l == 0)
      return _Ck;
    else
      {
	for (int __k = 0; __k < __l; ++__k)
	  _Ck *= std::hypot(_Tp(__k + 1), __eta) / _Tp(__k) / _Tp(__k + 1);
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
  template<typename _Tp>
    std::pair<_Tp, _Tp>
    __coulomb_f_recur(unsigned int __l_min, unsigned int __k_max,
		      _Tp __eta, _Tp __x,
		      _Tp _F_l_max, _Tp _Fp_l_max)
    {
      const auto __x_inv = _Tp{1}/__x;
      auto __fcl = _F_l_max;
      auto __fpl = _Fp_l_max;
      const auto __l_max = __l_min + __k_max;
      auto __l = __l_max;

      for (int __k = __k_max - 1; __k >= 0; --__k)
	{
	  const auto __el = __eta / __l;
	  const auto __rl = std::hypot(_Tp{1}, __el);
	  const auto __sl = __el  + __l * __x_inv;
	  const auto __fc_lm1 = (__fcl * __sl + __fpl) / __rl;
	  __fpl = __fc_lm1 * __sl - __fcl * __rl;
	  __fcl = __fc_lm1;
	  __l -= 1;
	}

      return std::make_pair(__fcl, __fpl);
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
  template<typename _Tp>
    std::pair<_Tp, _Tp>
    __coulomb_g_recur(unsigned int __l_min, unsigned int __k_max,
		      _Tp __eta, _Tp __x,
		      _Tp _G_l_min, _Tp _Gp_l_min)
    {
      const auto __x_inv = _Tp{1} / __x;
      auto __gcl = _G_l_min;
      auto __gpl = _Gp_l_min;
      auto __l = __l_min + 1;

      for (int k = 1; k <= __k_max; ++k)
	{
	  const auto __el = __eta / __l;
	  const auto __rl = std::hypot(_Tp{1}, __el);
	  const auto __sl = __el + __l * __x_inv;
	  const auto __gc_lm1 = (__sl * __gcl - __gpl) / __rl;
	  __gpl = __rl * __gcl - __sl * __gc_lm1;
	  __gcl = __gc_lm1;
	  __l += 1;
	}

      return std::make_pair(__gcl, __gpl);
    }

  /**
   * Evaluate the first continued fraction, giving the ratio F'/F
   * at the upper l value.
   * We also determine the sign of F at that point, since it is the sign
   * of the last denominator in the continued fraction.
   */
  template<typename _Tp>
    std::pair<_Tp, _Tp>
    __coulomb_CF1(unsigned int __l,
        	  _Tp __eta, _Tp __x)
    {
      const auto _CF1_small = 1.e-30;
      const auto _CF1_abort = 1.0e+05;
      const auto _CF1_acc = _Tp{2} * std::numeric_limits<_Tp>::epsilon();
      const auto __x_inv = _Tp{1} / __x;
      const auto __px = __l + _Tp{1} + _CF1_abort;

      auto __pk = __l + 1;
      auto __F = __eta / __pk + __pk * __x_inv;

      auto __fcl_sign = _Tp{1};

      if (std::abs(__F) < _CF1_small)
	__F = _CF1_small;
      auto __D = _Tp{0};
      auto __C = __F;

      _Tp __df;
      do
	{
	  const auto __pk1 = __pk + 1;
	  const auto __ek  = __eta / __pk;
	  const auto __rk2 = _Tp{1} + __ek * __ek;
	  const auto __tk  = (__pk + __pk1) * (__x_inv + __ek / __pk1);
	  __D = __tk - __rk2 * __D;
	  __C = __tk - __rk2 / __C;
	  if (std::abs(__C) < _CF1_small)
	    __C = _CF1_small;
	  if (std::abs(__D) < _CF1_small)
	    __D = _CF1_small;
	  __D = _Tp{1} / __D;
	  __df = __D * __C;
	  __F *= __df;
	  if (__D < _Tp{0})
	    {
	      // sign of result depends on sign of denominator.
	      __fcl_sign = -__fcl_sign;
	    }
	  __pk = __pk1;
	  if (__pk > __px)
	    std::__throw_runtime_error("__coulomb_CF1: Too many iterations.");
	}
      while (std::abs(__df - _Tp{1}) > _CF1_acc);

      return std::make_pair(__F, __fcl_sign);
    }

  /**
   * Evaluate the second continued fraction to obtain the ratio
   * @f[
   *    (G' + i F') / (G + i F) := P + i Q
   * @f]
   * at the specified l value.
   */
  template<typename _Tp>
    std::complex<_Tp>
    __coulomb_CF2(unsigned int __l, _Tp __eta, _Tp __x)
    {
      const auto _S_i = std::complex{0, 1};
      const auto _CF2_acc = _Tp{4} * std::numeric_limits<_Tp>::epsilon();
      const auto _CF2_abort = 2.0e+05;

      const auto __wi = _Tp{2} * __eta;
      const auto __x_inv = _Tp{1} / __x;
      const auto __e2mm1 = __eta * __eta + _Tp(__l * (__l + 1));

      auto __a = std::complex<_Tp>(-__e2mm1, __eta);
      auto __b = _Tp{2} * std::complex<_Tp>(__x - __eta, _Tp{2});
      auto __d = std::conj(__b) / std::norm(__b);

      auto __dpq = __x_inv * std::conj(__a * __d);

      auto __pk = _Tp{0};
      auto __PQ = std::complex<_Tp>(_Tp{0}, _Tp{1} - __eta * __x_inv);

      do
        {
	  __PQ += __dpq;
	  __pk += _Tp{2};
	  __a += std::complex(__pk, __wi);
	  __b += _Tp{2} * _S_i;
	  __d = __b + __a * __d;
	  __d = std::conj(__d) / std::norm(__d);
	  __dpq *= __b * __d - _Tp{1};
	  if (__pk > _CF2_abort)
	    std::__throw_runtime_error("__coulomb_CF2: Too many iterations.");
        }
      while (std::abs(__dpq) > (std::abs(__PQ)) * _CF2_acc);

      //if (Q < CF2_abort * std::numeric_limits<_Tp>::epsilon() * std::abs(P))
	//status = GSL_ELOSS;

      return __PQ;
    }

  /**
   * Return the bound-state Coulomb wave-function.
   */
  template <typename _Tp>
    std::complex<_Tp>
    __hydrogen(unsigned int __n,
               unsigned int __l, unsigned int __m,
               _Tp __Z, _Tp __r, _Tp __theta, _Tp __phi)
    {
      const auto _S_NaN = __gnu_cxx::__quiet_NaN(__r);

      if (__isnan(__Z) || __isnan(__r) || __isnan(__theta) || __isnan(__phi))
	return std::complex<_Tp>{_S_NaN, _S_NaN};
      else if(__n < 1)
	std::__throw_domain_error(__N("__hydrogen: "
				      "level number less than one"));
      else if(__l > __n - 1)
	std::__throw_domain_error(__N("__hydrogen: "
				      "angular momentum number too large"));
      else if(__Z <= _Tp(0))
	std::__throw_domain_error(__N("__hydrogen: non-positive charge"));
      else if(__r < _Tp(0))
	std::__throw_domain_error(__N("__hydrogen: negative radius"));
      else
	{
	  const auto __A = _Tp(2) * __Z / __n;

	  const auto __pre = std::sqrt(__A * __A * __A / (_Tp(2) * __n));
	  const auto __ln_a = __log_gamma(__n + __l + 1);
	  const auto __ln_b = __log_gamma(__n - __l);
	  const auto __ex = std::exp((__ln_b - __ln_a) / _Tp(2));
	  const auto __norm = __pre * __ex;

	  const auto __rho = __A * __r;
	  const auto __ea = std::exp(-__rho / _Tp(2));
	  const auto __pp = std::pow(__rho, __l);
	  const auto __lag = __assoc_laguerre(__n - __l - 1, 2 * __l + 1,
                                        	__rho);
	  const auto __sphh = __sph_legendre(__l, __m, __theta)
 			    * std::polar(_Tp(1), _Tp(__m) * __phi);

	  const auto __psi = __norm * __ea * __pp * __lag * __sphh;

	  return __psi;
	}
    }

_GLIBCXX_END_NAMESPACE_VERSION
} // namespace __detail
}

#endif // _GLIBCXX_BITS_SF_HYDROGEN_TCC
