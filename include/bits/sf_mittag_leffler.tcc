// Special functions -*- C++ -*-

// Copyright (C) 2019 Free Software Foundation, Inc.
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

/** @file bits/sf_mittag_leffler.tcc
 *  This is an internal header file, included by other library headers.
 *  Do not attempt to use it directly. @headername{cmath}
 */

#ifndef _GLIBCXX_BITS_SF_MITTAG_LEFFLER_TCC
#define _GLIBCXX_BITS_SF_MITTAG_LEFFLER_TCC 1

#pragma GCC system_header

#include <ext/integration.h>

namespace std _GLIBCXX_VISIBILITY(default)
{
_GLIBCXX_BEGIN_NAMESPACE_VERSION

// Implementation-space details.
namespace __detail
{

  /* Monotone integrand for the Mittag-Leffler function. */
  template<typename _Tp>
    std::complex<_Tp>
    __mittag_leffler_K(_Tp __alpha, _Tp __beta, _Tp __chi,
		       const std::complex<_Tp>& __z)
    {
      const auto _S_pi = __gnu_cxx::__const_pi(__alpha);
      const auto __chip1 = std::pow(__chi, _Tp{1} / __alpha);
      const auto __chip2 = std::pow(__chi, (_Tp{1} - __beta) / __alpha);
      return __chip2
	   * std::exp(-__chip1)
	   * (__chi * std::sin(_S_pi * (_Tp{1} - __beta))
	      - __z * std::sin(_S_pi * (_Tp{1} - __beta + __alpha)))
	   / (__chi * __chi - _Tp{2} * __chi * __z * std::cos(__alpha * _S_pi)
		 + __z * __z)
	   / _S_pi / __alpha;
    }

  /* Monotone integral for the Mittag-Leffler function. */
  template<typename _Tp>
    std::complex<_Tp>
    __mittag_leffler_K_integral(_Tp __alpha, _Tp __beta,
				_Tp __chi_min, _Tp __chi_max,
				const std::complex<_Tp>& __z)
    {
      const auto _S_eps = __gnu_cxx::__epsilon(__chi_min);
      auto __func = [__alpha, __beta, __z](_Tp __chi)
		    -> std::complex<_Tp>
		    { return __mittag_leffler_K(__alpha, __beta, __chi, __z); };

      const auto __epsabs = _Tp{100} * _S_eps;
      const auto __epsrel = _Tp{0};
      auto __ws = __gnu_cxx::cquad_workspace<_Tp, std::complex<_Tp>>();

      auto __quad
	= __gnu_cxx::cquad_integrate(__ws, __func, __chi_min, __chi_max,
				     __epsabs, __epsrel);

      return __quad.__result;
    }

  /* Oscillatory integrand for the Mittag-Leffler function. */
  template<typename _Tp>
    std::complex<_Tp>
    __mittag_leffler_P(_Tp __alpha, _Tp __beta, _Tp __epsilon, _Tp __phi,
		       const std::complex<_Tp>& __z)
    {
      const auto _S_i = std::complex<_Tp>{0, 1};
      const auto _S_pi = __gnu_cxx::__const_pi(__alpha);
      const auto __epsp1 = std::pow(__epsilon, _Tp{1} / __alpha);
      const auto __rat = _Tp{1} + (_Tp{1} - __beta) / __alpha;
      const auto __epsp2 = std::pow(__epsilon, __rat);
      const auto __omega = __phi * __rat + __epsp1 * std::sin(__phi / __alpha);
      return __epsp2
	   * std::exp(__epsp1 * std::cos(__phi / __alpha))
	   * std::polar(_Tp{1}, __omega)
	   / (__epsilon * _S_i - __z)
	   / _Tp{2} / _S_pi / __alpha;
    }

  /* Oscillatory integral for the Mittag-Leffler function. */
  template<typename _Tp>
    std::complex<_Tp>
    __mittag_leffler_P_integral(_Tp __alpha, _Tp __beta, _Tp __epsilon,
				_Tp __phi_min, _Tp __phi_max,
				const std::complex<_Tp>& __z)
    {
      const auto _S_eps = __gnu_cxx::__epsilon(__phi_min);
      auto __func = [__alpha, __beta, __epsilon, __z](_Tp __phi)
		    -> std::complex<_Tp>
		    {
		      return __mittag_leffler_P(__alpha, __beta,
						__epsilon, __phi, __z);
		    };

      const auto __epsabs = _Tp{100} * _S_eps;
      const auto __epsrel = _Tp{0};
      auto __ws = __gnu_cxx::cquad_workspace<_Tp, std::complex<_Tp>>();

      auto __quad
	= __gnu_cxx::cquad_integrate(__ws, __func, __phi_min, __phi_max,
				     __epsabs, __epsrel);

      return __quad.__result;
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
  template<typename _Tp>
    std::complex<_Tp>
    __mittag_leffler(_Tp __alpha, _Tp __beta, const std::complex<_Tp>& __z)
    {
      using _Cmplx = std::complex<_Tp>;
      const auto _S_eps = __gnu_cxx::__epsilon(__alpha);
      const auto _S_2pi = __gnu_cxx::__const_2_pi(__alpha);
      const auto _S_pi = __gnu_cxx::__const_pi(__alpha);

      const auto __az = std::abs(__z);
      if (__alpha == _Tp{0})
	{
	  if (__az < _Tp{1})
	    return __gamma_reciprocal(__beta) / (_Tp{1} - __z);
	  else
	    {
	      const auto _S_inf = __gnu_cxx::__infinity(__alpha);
	      return std::complex<_Tp>{_S_inf, _S_inf};
	    }
	}
      else if (__alpha < _Tp{0})
	return __gamma_reciprocal(__beta)
	     - __mittag_leffler(-__alpha, __beta, _Tp{1} / __z);
      else if (__alpha > _Tp{1})
	{
          unsigned int __k0 = _Tp{1} + std::floor(__alpha);
	  const auto __alpha0 = __alpha / __k0;
	  const auto __rho0 = std::pow(__z, _Tp{1} / _Tp(__k0));
	  const auto __lamb = _S_2pi / _Tp(__k0);

	  auto __E = _Cmplx{0};
	  for (auto __k = 0u; __k < __k0; ++__k)
	    {
	      auto __zk = __rho0 * std::polar(_Tp{1}, __lamb * _Tp(__k));
	      __E += __mittag_leffler(__alpha0, __beta, __zk);
	    }
	  return __E / _Tp(__k0);
	}
      else if (__az < _S_eps)
	return __gamma_reciprocal(__beta);
      else if (__az < _Tp{1})
	{
	  unsigned int __k0 = std::max(std::ceil((_Tp{1} - __beta) / __alpha),
				std::ceil(std::log(_S_eps * (_Tp{1} - __az))
					    / std::log(__az)));
	  auto __E = _Cmplx{0};
	  auto __zk = _Cmplx{1};
	  for (auto __k = 0u; __k <= __k0; ++__k)
	    {
	      const auto __arg = __beta + __alpha * __k;
	      const auto __term = __zk * __gamma_reciprocal(__arg);
	      __E += __term;
	      if (std::abs(__term) < _S_eps)
		break;
	      __zk *= __z;
	    }
	  return __E;
	}
      else if (__az > std::floor(_Tp{10} + _Tp{5} * __alpha))
	{
	  unsigned int __k0 = std::floor(-std::log(_S_eps) / std::log(__az));
	  auto __E = _Cmplx{0};
	  auto __zk = _Cmplx{1};
	  for (auto __k = 1u; __k <= __k0; ++__k)
	    {
	      __zk /= __z;
	      __E += __zk * __gamma_reciprocal(__beta - __alpha * __k);
	    }
	  if (std::arg(__z)
	      < _S_pi * (__alpha / _Tp{4} + std::min(_Tp{1}, __alpha) / _Tp{2}))
	    {
	      const auto __zp1 = std::pow(__z, _Tp{1} / __alpha);
	      const auto __zp2 = std::pow(__z, (_Tp{1} - __beta) / __alpha);
	      const auto __extra = __zp2 * std::exp(__zp1) / __alpha;
	      return __extra - __E;
	    }
	  else
	    return -__E;
	}
      else
	{
	  auto __chi0 = _Tp{0};
	  if (__beta >= _Tp{0})
	    __chi0 = std::max({_Tp{1}, _Tp{2} * __az,
			std::pow(-std::log(_S_pi * _S_eps / _Tp{6}), __alpha)});
	  else
	    {
	      const auto __abeta = std::abs(__beta);
	      __chi0 = std::max({std::pow(_Tp{1} + __abeta, __alpha),
			 _Tp{2} * __az, 
		std::pow(-_Tp{2}
			  * std::log(_S_pi * _S_eps
			  / (_Tp{6} * (__abeta + _Tp{2})
			   * std::pow(_Tp{2} * __abeta, __abeta))),
			 __alpha)});
	    }

	  const auto __absarz = std::abs(std::arg(__z));
	  if (__absarz > __alpha * _S_pi + _S_eps)
	    {
	      if (__beta <= _Tp{1})
		return __mittag_leffler_K_integral(__alpha, __beta,
						   _Tp{0}, __chi0, __z);
	      else
		{
		  const auto __api = _S_pi * __alpha;
		  return __mittag_leffler_K_integral(__alpha, __beta,
						     _Tp{1}, __chi0, __z)
		       + __mittag_leffler_P_integral(__alpha, __beta, _Tp{1},
						     -__api, __api, __z);
		}
	    }
	  else if (__absarz < __alpha * _S_pi - _S_eps)
	    {
	      const auto __zp1 = std::pow(__z, _Tp{1} / __alpha);
	      const auto __zp2 = std::pow(__z, (_Tp{1} - __beta) / __alpha);
	      const auto __extra = __zp2 * std::exp(__zp1) / __alpha;
	      if (__beta <= _Tp{1})
		return __mittag_leffler_K_integral(__alpha, __beta,
						   _Tp{0}, __chi0, __z)
			+ __extra;
	      else
		{
		  const auto __lo = __az / _Tp{2};
		  const auto __api = _S_pi * __alpha;
		  return __mittag_leffler_K_integral(__alpha, __beta,
						     __lo, __chi0, __z)
		       + __mittag_leffler_P_integral(__alpha, __beta, __lo,
							-__api, __api, __z)
		       + __extra;
		}
	    }
	  else
	    {
	      const auto __lo = (__az + _Tp{1}) / _Tp{2};
	      const auto __api = _S_pi * __alpha;
	      return __mittag_leffler_K_integral(__alpha, __beta,
						 __lo, __chi0, __z)
		   + __mittag_leffler_P_integral(__alpha, __beta, __lo,
						 -__api, __api, __z);
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
  template<typename _Tp>
    std::complex<_Tp>
    __mittag_leffler_deriv(_Tp __alpha, _Tp __beta,
			   const std::complex<_Tp>& __z)
    {
      using _Cmplx = std::complex<_Tp>;
      const auto _S_eps = __gnu_cxx::__epsilon(__alpha);

      const auto __az = std::abs(__z);
      if (__az < _Tp{1})
	{
	  auto __k1 = _Tp{0};
	  if (__alpha > _Tp{1})
	    __k1 = _Tp{1} + (_Tp{2} - __alpha - __beta) / (__alpha - _Tp{1});
	  else
	    {
	      const auto __D = _Tp{1}
			     + __alpha * (__alpha - _Tp{4} * __beta + _Tp{6});
	      const auto __omega = __alpha + __beta - _Tp{3} / _Tp{2};
	      const auto __rat = _Tp{1} + (_Tp{3} - __alpha - __beta) / __alpha;
	      if (__D <= _Tp{0})
		__k1 = __rat;
	      else
		__k1 = std::max(__rat,
			_Tp{1}
			+ (_Tp{1} - _Tp{2} * __omega * __alpha + std::sqrt(__D))
				  / (2 * __alpha * __alpha));
	    }
	  __k1 = std::ceil(__k1);
	  unsigned int __k0 = std::max(__k1,
				 std::ceil(std::log(_S_eps * (_Tp{1} - __az))
					 / std::log(__az)));
	  auto __Ep = _Cmplx{0};
	  auto __zk = _Cmplx{1};
	  for (auto __k = 0u; __k <= __k0; ++__k)
	    {
	      __Ep += _Tp(__k + 1) * __zk
		    * __gamma_reciprocal(__beta + __alpha * _Tp(__k + 1));
	      __zk *= __z;
	    }
	  return __Ep;
	}
      else
	return (__mittag_leffler(__alpha, __beta - _Tp{1}, __z)
	      - (__beta - _Tp{1}) * __mittag_leffler(__alpha, __beta, __z))
	     / __alpha / __z;
    }

} // namespace __detail

_GLIBCXX_END_NAMESPACE_VERSION
} // namespace std

#endif // _GLIBCXX_BITS_SF_MITTAG_LEFFLER_TCC
