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

/** @file bits/sf_trigint.tcc
 *  This is an internal header file, included by other library headers.
 *  Do not attempt to use it directly. @headername{cmath}
 */

#ifndef _GLIBCXX_BITS_SF_TRIGINT_TCC
#define _GLIBCXX_BITS_SF_TRIGINT_TCC 1

#pragma GCC system_header

#include <complex>
#include <ext/math_const.h>

namespace std _GLIBCXX_VISIBILITY(default)
{
// Implementation-space details.
namespace __detail
{
_GLIBCXX_BEGIN_NAMESPACE_VERSION

  enum
  {
    SININT,
    COSINT
  };

  /**
   *  @brief This function computes the sine @f$ Si(x) @f$
   *         and cosine @f$ Ci(x) @f$ integrals by continued fraction
   *         for positive argument.
   */
  template<typename _Tp>
    void
    __sincosint_cont_frac(_Tp __t, _Tp& _Si, _Tp& _Ci)
    {
      constexpr auto _S_max_iter = 100;
      constexpr auto _S_eps = _Tp{5} * __gnu_cxx::__epsilon<_Tp>();
      constexpr auto _S_fp_min = __gnu_cxx::__min<_Tp>();
      constexpr auto _S_pi_2 = __gnu_cxx::__math_constants<_Tp>::__pi_half;

      // Evaluate Ci and Si by Lentz's modified method of continued fractions.
      std::complex<_Tp> __b(_Tp{1}, __t);
      std::complex<_Tp> __c(_Tp{1} / _S_fp_min);
      std::complex<_Tp> __d(_Tp{1} / __b);
      std::complex<_Tp> __h(__d);
      int __i = 2;
      while (true)
	{
	  auto __a = -_Tp(__i - 1) * _Tp(__i - 1);
	  __b += _Tp{2};
	  __d = _Tp{1} / (__a * __d + __b);
	  __c = __b + __a / __c;
	  std::complex<_Tp> __del = __c * __d;
	  __h *= __del;
	  if (std::abs(__del - _Tp{1}) < _S_eps)
	    break;
	  if (__i > _S_max_iter)
	    std::__throw_runtime_error(__N("__sincosint_cont_frac: "
				      "continued fraction evaluation failed"));
	  ++__i;
	}
      __h *= std::polar(_Tp{1}, -__t);
      _Ci = -__h.real();
      _Si = _S_pi_2 + __h.imag();

      return;
    }


  /**
   *  @brief This function computes the sine @f$ Si(x) @f$
   *         and cosine @f$ Ci(x) @f$ integrals by series summation
   *         for positive argument.
   */
  template<typename _Tp>
    void
    __sincosint_series(_Tp __t, _Tp& _Si, _Tp& _Ci)
    {
      constexpr auto _S_max_iter = 100;
      constexpr auto _S_eps = _Tp{5} * __gnu_cxx::__epsilon<_Tp>();
      constexpr auto _S_fp_min = __gnu_cxx::__min<_Tp>();
      constexpr auto _S_gamma_e
	= __gnu_cxx::__math_constants<_Tp>::__gamma_e;

      // Evaluate Ci and Si by series simultaneously.
      _Tp __sumc(0), __sums(0);
      if (__t * __t < _S_fp_min)
	{
	  // Avoid underflow.
	  __sumc = _Tp{0};
	  __sums = __t;
	}
      else
	{
	  // Evaluate Si and Ci by series expansion.
	  _Tp __sum(0);
	  _Tp __sign(1), __fact(1);
	  bool __odd = true;
	  unsigned int __k = 1;
	  while (true)
	    {
	      __fact *= __t / __k;
	      _Tp __term = __fact / __k;
	      __sum += __sign * __term;
	      _Tp __err = __term / std::abs(__sum);
	      if (__odd)
		{
		  __sign = -__sign;
		  __sums = __sum;
		  __sum = __sumc;
		}
	      else
		{
		  __sumc = __sum;
		  __sum = __sums;
		}
	      if (__err < _S_eps)
		break;
	      __odd = !__odd;
	      ++__k;
	      if (__k > _S_max_iter)
		std::__throw_runtime_error(__N("__sincosint_series: "
					       "series evaluation failed"));
	    }
	}
      _Si = __sums;
      _Ci = _S_gamma_e + std::log(__t) + __sumc;

      return;
    }


  /**
   *  @brief This function computes the sine @f$ Si(x) @f$
   *         and cosine @f$ Ci(x) @f$ integrals by asymptotic series summation
   *         for positive argument.
   *
   *   The asymptotic series is very good for x > 50.
   */
  template<typename _Tp>
    void
    __sincosint_asymp(_Tp __t, _Tp& _Si, _Tp& _Ci)
    {
      constexpr auto _S_max_iter = 100;
      constexpr auto _S_eps = _Tp{5} * __gnu_cxx::__epsilon<_Tp>();
      constexpr auto _S_pi_2 = __gnu_cxx::__math_constants<_Tp>::__pi_half;

      auto __invt = _Tp{1} / __t;
      auto __term = _Tp{1}; // 0!
      auto __sume = _Tp{__term};
      __term *= __invt; // 1! / t
      auto __sumo = _Tp{__term};
      auto __sign = _Tp{1};
      auto __even = true;
      auto __k = 2;
      while (true)
	{
	  __term *= __k * __invt;

	  if (__even)
	    {
	      __sign = -__sign;
	      __sume += __sign * __term;

	    }
	  else
	    {
	      __sumo += __sign * __term;
	      if (__term / std::abs(__sumo) < _S_eps)
		break;
	    }

	  __even = !__even;

	  if (__k > _S_max_iter)
	    std::__throw_runtime_error(__N("__sincosint_asymp: "
	    				   "Series evaluation failed"));
	  ++__k;
	}

      _Si = _S_pi_2
	   - std::cos(__t) * __invt * __sume
	   - std::sin(__t) * __invt * __sumo;
      _Ci = std::sin(__t) * __invt * __sume
	   - std::cos(__t) * __invt * __sumo;

      return;
    }


  /**
   *  @brief This function returns the sine @f$ Si(x) @f$
   *         and cosine @f$ Ci(x) @f$ integrals as a @c pair.
   *
   *  The sine integral is defined by:
   *  @f[
   *      Si(x) = \int_0^x dt \frac{\sin(t)}{t}
   *  @f]
   *
   *  The cosine integral is defined by:
   *  @f[
   *      Ci(x) = \gamma_E + \log(x) + \int_0^x dt \frac{\cos(t) - 1}{t}
   *  @f]
   */
  template<typename _Tp>
    std::pair<_Tp, _Tp>
    __sincosint(_Tp __x)
    {
      constexpr auto _S_NaN = __gnu_cxx::__quiet_NaN<_Tp>();
      if (__isnan(__x))
	return std::make_pair(_S_NaN, _S_NaN);

      auto __t = std::abs(__x);
      _Tp _Ci, _Si;
      if (__t == _Tp{0})
	{
	  _Si = _Tp{0};
	  _Ci = -__gnu_cxx::__infinity<_Tp>();
	}
      else if (__t > _Tp{1000}) // Check this!
	__sincosint_asymp(__t, _Si, _Ci);
      else if (__t > _Tp{2})
	__sincosint_cont_frac(__t, _Si, _Ci);
      else
	__sincosint_series(__t, _Si, _Ci);

      if (__x < _Tp{0})
	_Si = -_Si;

      return std::make_pair(_Si, _Ci);
    }

_GLIBCXX_END_NAMESPACE_VERSION
} // namespace __detail
} // namespace std

#endif // _GLIBCXX_BITS_SF_TRIGINT_TCC
