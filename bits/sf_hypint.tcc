// Special functions -*- C++ -*-

// Copyright (C) 2016 Free Software Foundation, Inc.
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

/** @file bits/sf_hypint.tcc
 *  This is an internal header file, included by other library headers.
 *  Do not attempt to use it directly. @headername{cmath}
 */

#ifndef _GLIBCXX_BITS_SF_HYPINT_TCC
#define _GLIBCXX_BITS_SF_HYPINT_TCC 1

#pragma GCC system_header

#include <complex>
#include <ext/math_const.h>

namespace std _GLIBCXX_VISIBILITY(default)
{
// Implementation-space details.
namespace __detail
{
_GLIBCXX_BEGIN_NAMESPACE_VERSION

  /**
   *  @brief This function computes the hyperbolic cosine @f$ Chi(x) @f$
   *    and hyperbolic sine @f$ Shi(x) @f$ integrals
   *    by continued fraction for positive argument.
   */
  ////FIXME!!!!
  template<typename _Tp>
    void
    __chshint_cont_frac(_Tp __t, _Tp& _Chi, _Tp& _Shi)
    {
      constexpr unsigned int _S_max_iter = 100;
      constexpr auto _S_eps = _Tp{5} * __gnu_cxx::__epsilon<_Tp>();
      constexpr auto _S_fp_min = __gnu_cxx::__min<_Tp>();
      constexpr auto _S_pi_2 = __gnu_cxx::__math_constants<_Tp>::__pi_half;

      // Evaluate Chi and Shi by Lentz's modified method of continued fracions.
      std::complex<_Tp> __b(_Tp{1}, __t);
      std::complex<_Tp> __c(_Tp{1} / _S_fp_min);
      std::complex<_Tp> __d(_Tp{1} / __b);
      std::complex<_Tp> __h(__d);
      unsigned int __i = 2;
      while (true)
	{
	  _Tp __a = -(__i - 1) * (__i - 1);
	  __b += _Tp{2};
	  __d = _Tp{1} / (__a * __d + __b);
	  __c = __b + __a / __c;
	  auto __del = __c * __d;
	  __h *= __del;
	  if (std::abs(__del.real() - _Tp{1}) + std::abs(__del.imag()) < _S_eps)
	    break;
	  if (__i > _S_max_iter)
	    std::__throw_runtime_error(__N("__chshint_cont_frac: "
				       "Continued fraction evaluation failed"));
	  ++__i;
	}
      __h *= std::polar(_Tp{1}, -__t);
      _Chi = -__h.real();
      _Shi = _S_pi_2 + __h.imag();

      return;
    }


  /**
   *  @brief This function computes the hyperbolic cosine @f$ Chi(x) @f$
   *    and hyperbolic sine @f$ Shi(x) @f$ integrals
   *    by series summation for positive argument.
   */
  template<typename _Tp>
    void
    __chshint_series(_Tp __t, _Tp& _Chi, _Tp& _Shi)
    {
      constexpr auto _S_max_iter = 100;
      constexpr auto _S_eps = _Tp{5} * __gnu_cxx::__epsilon<_Tp>();
      constexpr auto _S_fp_min = __gnu_cxx::__min<_Tp>();
      constexpr auto _S_gamma_e = __gnu_cxx::__math_constants<_Tp>::__gamma_e;

      // Evaluate Chi and Shi by series simultaneously.
      _Tp _Csum(0), _Ssum(0);
      if (__t * __t < _S_fp_min)
	{
	  // Avoid underflow.
	  _Csum = _Tp{0};
	  _Ssum = __t;
	}
      else
	{
	  // Evaluate Shi and Chi by series expansion.
	  _Tp __sum(0);
	  _Tp __fact(1);
	  auto __odd = true;
	  auto __k = 1;
	  while (true)
	    {
	      __fact *= __t / __k;
	      _Tp __term = __fact / __k;
	      __sum += __term;
	      _Tp __err = __term / std::abs(__sum);
	      if (__odd)
		{
		  _Ssum = __sum;
		  __sum = _Csum;
		}
	      else
		{
		  _Csum = __sum;
		  __sum = _Ssum;
		}
	      if (__err < _S_eps)
		break;
	      __odd = !__odd;
	      ++__k;
	      if (__k > _S_max_iter)
		std::__throw_runtime_error(__N("__chshint_series: "
					       "Series evaluation failed"));
	    }
	}
      _Chi = _S_gamma_e + std::log(__t) + _Csum;
      _Shi = _Ssum;

      return;
    }


  /**
   *  @brief This function returns the hyperbolic cosine @f$ Ci(x) @f$
   *    and hyperbolic sine @f$ Si(x) @f$ integrals as a pair.
   *
   *  The hyperbolic cosine integral is defined by:
   *  @f[
   *      Chi(x) = \gamma_E + \log(x) + \int_0^x dt \frac{\cosh(t) - 1}{t}
   *  @f]
   *
   *  The hyperbolic sine integral is defined by:
   *  @f[
   *      Shi(x) = \int_0^x dt \frac{\sinh(t)}{t}
   *  @f]
   */
  template<typename _Tp>
    std::pair<_Tp, _Tp>
    __chshint(_Tp __x, _Tp& _Chi, _Tp& _Shi)
    {
      constexpr auto _S_NaN = __gnu_cxx::__quiet_NaN<_Tp>();
      if (__isnan(__x))
	return std::make_pair(_S_NaN, _S_NaN);

      auto __t = std::abs(__x);
      if (__t == _Tp{0})
	{
	  _Chi = -__gnu_cxx::__infinity<_Tp>();
	  _Shi = _Tp{0};
	}
      else if (__t > _Tp{2})
	__chshint_cont_frac(__t, _Chi, _Shi);
      else
	__chshint_series(__t, _Chi, _Shi);

      if (__x < _Tp{0})
	_Shi = -_Shi;

      return std::make_pair(_Chi, _Shi);
    }

_GLIBCXX_END_NAMESPACE_VERSION
} // namespace __detail
}

#endif // _GLIBCXX_BITS_SF_HYPINT_TCC
