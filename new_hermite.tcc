// Special functions -*- C++ -*-

// Copyright (C) 2016
// Free Software Foundation, Inc.
//
// This file is part of the GNU ISO C++ Library.  This library is free
// software; you can redistribute it and/or modify it under the
// terms of the GNU General Public License as published by the
// Free Software Foundation; either version 2, or (at your option)
// any later version.
//
// This library is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License along
// with this library; see the file COPYING.  If not, write to the Free
// Software Foundation, 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301,
// USA.
//
// As a special exception, you may use this file as part of a free software
// library without restriction.  Specifically, if other files instantiate
// templates or use macros or inline functions from this file, or you compile
// this file and link it with other files to produce an executable, this
// file does not by itself cause the resulting executable to be covered by
// the GNU General Public License.  This exception does not however
// invalidate any other reasons why the executable file might be covered by
// the GNU General Public License.

/** @file tr1/sf_hermite.tcc
 *  This is an internal header file, included by other library headers.
 *  You should not attempt to use it directly.
 */

//
// ISO C++ TR29124: Special functions
//

#ifndef _GLIBCXX_BITS_NEW_HERMITE_TCC
#define _GLIBCXX_BITS_NEW_HERMITE_TCC 1



  /**
   *   @brief This routine returns the Hermite polynomial
   *          of order n: \f$ H_n(x) \f$ by recursion on n.
   *
   *   The Hermite polynomial is defined by:
   *   @f[
   *     H_n(x) = (-1)^n e^{x^2} \frac{d^n}{dx^n} e^{-x^2}
   *   @f]
   *
   *   @param __n The order of the Hermite polynomial.
   *   @param __x The argument of the Hermite polynomial.
   *   @return The value of the Hermite polynomial of order n
   *           and argument x.
   */
  template<typename _Tp>
    _Tp
    __poly_hermite_recursion(unsigned int __n, _Tp __x)
    {
      // Compute H_0.
      auto __H_nm2 = _Tp{1};
      if (__n == 0)
	return __H_nm2;

      // Compute H_1.
      auto __H_nm1 = _Tp{2} * __x;
      if (__n == 1)
	return __H_nm1;

      // Compute H_n.
      _Tp __H_n;
      for (unsigned int __i = 2; __i <= __n; ++__i)
	{
	  __H_n = _Tp{2} * (__x * __H_nm1 - _Tp(__i - 1) * __H_nm2);
	  __H_nm2 = __H_nm1;
	  __H_nm1 = __H_n;
	}

      return __H_n;
    }


  /**
   *   @brief This routine returns the Hermite polynomial
   *          of large order n: \f$ H_n(x) \f$.  We assume here
   *          that x >= 0.
   *
   *   The Hermite polynomial is defined by:
   *   @f[
   *     H_n(x) = (-1)^n e^{x^2} \frac{d^n}{dx^n} e^{-x^2}
   *   @f]
   *
   *  @see "Asymptotic analysis of the Hermite polynomials
   *        from their differential-difference equation", 
   *        Diego Dominici, arXiv:math/0601078v1 [math.CA] 4 Jan 2006
   *   @param __n The order of the Hermite polynomial.
   *   @param __x The argument of the Hermite polynomial.
   *   @return The value of the Hermite polynomial of order n
   *           and argument x.
   */
  template<typename _Tp>
    _Tp
    __poly_hermite_asymp(unsigned int __n, _Tp __x)
    {
      constexpr auto _S_pi = __gnu_cxx::__math_constants<_Tp>::__pi;
      constexpr auto _S_sqrt_2 = __gnu_cxx::__math_constants<_Tp>::__root_2;
      constexpr auto _S_sqrt_2pi = _S_sqrt_2
				 * __gnu_cxx::__math_constants<_Tp>::__root_pi;
      // __x >= 0 in this routine.
      const auto __xturn = std::sqrt(_Tp(2 * __n));
      if (std::abs(__x - __xturn) < _Tp{0.05L} * __xturn)
	{
	  // Transition region x ~ sqrt(2n).
	  const auto __n_2 = _Tp(__n) / _Tp{2};
	  const auto __n6th = std::pow(_Tp(__n), _Tp{1} / _Tp{6});
	  const auto __exparg = _Tp(__n) * std::log(__xturn) - _Tp{3} * __n_2
			      + __xturn * __x;
	  const auto __airyarg = _S_sqrt_2 * (__x - __xturn) * __n6th;
	  _Tp _Ai, _Bi, _Aip, _Bip;
	  std::__detail::__airy(__airyarg, _Ai, _Bi, _Aip, _Bip);
	  return _S_sqrt_2pi * __n6th * std::exp(__exparg) * _Ai;
	}
      else if (__x < __xturn)
	{
	  // Oscillatory region |x| < sqrt(2n).
	  const auto __theta = std::asin(__x / __xturn);
	  const auto __2theta = _Tp{2} * __theta;
	  const auto __n_2 = _Tp(__n) / _Tp{2};
	  const auto __exparg = __n_2 * (_Tp{2} * std::log(__xturn)
					- std::cos(__2theta));
	  const auto __arg = __theta / _Tp{2}
			   + __n_2 * (std::sin(__2theta) + __2theta - _S_pi);
	  return std::sqrt(_Tp{2} / std::cos(__theta))
	       * std::exp(__exparg) * std::cos(__arg);
	}
      else
	{
	  // Exponential region |x| > sqrt(2n).
	  const auto __sigma = std::sqrt((__x - __xturn) * (__x + __xturn));
	  const auto __exparg = _Tp{0.5L} * (__x * (__x - __sigma) - __n)
			      + __n * std::log(__sigma + __x);
	  return std::exp(__exparg)
	       * std::sqrt(_Tp{0.5L} * (_Tp{1} + __x / __sigma));
	}
    }


  /**
   *  @brief  Compute the normalized Hermite polynomial by recursion.
   *  @todo  Tabulate sqrt(int) or even sqrt(i)/sqrt(i+1) and add a helper.
   */
  template<typename _Tp>
    _Tp
    __poly_hermite_norm_recursion(unsigned int __n, _Tp __x)
    {
      const auto _S_inv_root4_pi = _Tp{7.511255444649424828587030047762276930510e-1L};
      const auto _S_root_2 = __gnu_cxx::__math_constants<_Tp>::__root_2;

      auto __H_nm2 = _S_inv_root4_pi;
      if (__n == 0)
	return __H_nm2;
      auto __H_nm1 = _S_root_2 * __x * __H_nm2;
      if (__n == 1)
	return __H_nm1;

      _Tp __H_n;
      for (int __i = 2; __i < __n; ++__i)
	{
	  __H_n = __x * std::sqrt(_Tp{2} / _Tp(__i)) * __H_nm1
		- std::sqrt(_Tp(__i - 1) / _Tp(__i)) * __H_nm2;
	  __H_nm2 = __H_nm1;
	  __H_nm1 = __H_n;
	}

      return __H_n;
    }

  /**
   *  @brief  Compute the normalized probabalists Hermite polynomial by recursion.
   */
  template<typename _Tp>
    _Tp
    __poly_hermite_prob_recursion(unsigned int __n, _Tp __x)
    {
      // Compute He_0.
      auto __He_nm2 = _Tp{1};
      if (__n == 0)
	return __He_nm2;

      // Compute He_1.
      auto __He_nm1 = _Tp{2} * __x;
      if (__n == 1)
	return __He_nm1;

      // Compute H_n.
      _Tp __He_n;
      for (unsigned int __i = 2; __i <= __n; ++__i)
	{
	  __He_n = __x * __He_nm1 - _Tp(__i - 1) * __He_nm2;
	  __He_nm2 = __He_nm1;
	  __He_nm1 = __He_n;
	}

      return __He_n;
    }

  /**
   *  @brief  Compute the normalized probabalists Hermite polynomial by recursion.
   */
  template<typename _Tp>
    _Tp
    __poly_hermite_prob_norm_recursion(unsigned int __n, _Tp __x)
    {
      const auto _S_inv_root4_pi = _Tp{7.511255444649424828587030047762276930510e-1L};
      const auto _S_root_2 = __gnu_cxx::__math_constants<_Tp>::__root_2;

      // Compute He_0.
      auto __He_nm2 = _S_inv_root4_pi;
      if (__n == 0)
	return __He_nm2;

      // Compute He_1.
      auto __He_nm1 = _Tp{2} * __x;
      if (__n == 1)
	return __He_nm1;

      // Compute H_n.
      _Tp __He_n;
      for (unsigned int __i = 2; __i <= __n; ++__i)
	{
	  __He_n = (__x * __He_nm1 - _Tp(__i - 1) * __He_nm2) / std::sqrt(_Tp(__i));
	  __He_nm2 = __He_nm1;
	  __He_nm1 = __He_n;
	}

      return __He_n;
    }

#endif // _GLIBCXX_BITS_NEW_HERMITE_TCC
