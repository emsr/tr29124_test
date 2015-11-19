// Special functions -*- C++ -*-

// Copyright (C) 2006-2015 Free Software Foundation, Inc.
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

/** @file bits/sf_beta.tcc
 *  This is an internal header file, included by other library headers.
 *  Do not attempt to use it directly. @headername{cmath}
 */

//
// ISO C++ 14882 TR29124: Mathematical Special Functions
//

// Written by Edward Smith-Rowland.
//
// References:
//   (1) Handbook of Mathematical Functions,
//       ed. Milton Abramowitz and Irene A. Stegun,
//       Dover Publications,
//       Section 6, pp. 253-266
//   (2) The Gnu Scientific Library, http://www.gnu.org/software/gsl
//   (3) Numerical Recipes in C, by W. H. Press, S. A. Teukolsky,
//       W. T. Vetterling, B. P. Flannery, Cambridge University Press (1992),
//       2nd ed, pp. 213-216
//   (4) Gamma, Exploring Euler's Constant, Julian Havil,
//       Princeton, 2003.

#ifndef _GLIBCXX_BITS_SF_BETA_TCC
#define _GLIBCXX_BITS_SF_BETA_TCC 1

namespace std _GLIBCXX_VISIBILITY(default)
{
// Implementation-space details.
namespace __detail
{
_GLIBCXX_BEGIN_NAMESPACE_VERSION

  /**
   *   @brief  Return the beta function: \f$B(x,y)\f$.
   *
   *   The beta function is defined by
   *   @f[
   *     B(x,y) = \frac{\Gamma(x)\Gamma(y)}{\Gamma(x+y)}
   *   @f]
   *
   *   @param __x The first argument of the beta function.
   *   @param __y The second argument of the beta function.
   *   @return  The beta function.
   */
  template<typename _Tp>
    _Tp
    __beta_gamma(_Tp __x, _Tp __y)
    {

      _Tp __bet;
      if (__x > __y)
	{
	  __bet = __gamma(__x) / __gamma(__x + __y);
	  __bet *= __gamma(__y);
	}
      else
	{
	  __bet = __gamma(__y) / __gamma(__x + __y);
	  __bet *= __gamma(__x);
	}

      return __bet;
    }

  /**
   *   @brief  Return the beta function \f$B(x,y)\f$ using
   *           the log gamma functions.
   *
   *   The beta function is defined by
   *   @f[
   *     B(x,y) = \frac{\Gamma(x)\Gamma(y)}{\Gamma(x+y)}
   *   @f]
   *
   *   @param __x The first argument of the beta function.
   *   @param __y The second argument of the beta function.
   *   @return  The beta function.
   */
  template<typename _Tp>
    _Tp
    __beta_lgamma(_Tp __x, _Tp __y)
    {
      _Tp __bet = __log_gamma(__x)
		+ __log_gamma(__y)
		- __log_gamma(__x + __y);
      __bet = std::exp(__bet);
      return __bet;
    }


  /**
   *   @brief  Return the beta function \f$B(x,y)\f$ using
   *           the product form.
   *
   *   The beta function is defined by
   *   @f[
   *     B(x,y) = \frac{\Gamma(x)\Gamma(y)}{\Gamma(x+y)}
   *   @f]
   *
   *   @param __x The first argument of the beta function.
   *   @param __y The second argument of the beta function.
   *   @return  The beta function.
   */
  template<typename _Tp>
    _Tp
    __beta_product(_Tp __x, _Tp __y)
    {

      _Tp __bet = (__x + __y) / (__x * __y);

      const unsigned int __max_iter = 1000000;

      for (unsigned int __k = 1; __k < __max_iter; ++__k)
	{
	  _Tp __term = (_Tp{1} + (__x + __y) / __k)
		     / ((_Tp{1} + __x / __k) * (_Tp{1} + __y / __k));
	  __bet *= __term;
	}

      return __bet;
    }


  /**
   *   @brief  Return the beta function \f$ B(x,y) \f$.
   *
   *   The beta function is defined by
   *   @f[
   *     B(x,y) = \frac{\Gamma(x)\Gamma(y)}{\Gamma(x+y)}
   *   @f]
   *
   *   @param __x The first argument of the beta function.
   *   @param __y The second argument of the beta function.
   *   @return  The beta function.
   */
  template<typename _Tp>
    inline _Tp
    __beta(_Tp __x, _Tp __y)
    {
      if (__isnan(__x) || __isnan(__y))
	return __gnu_cxx::__math_constants<_Tp>::__NaN;
      else
	return __beta_lgamma(__x, __y);
    }


  template<typename _Tp>
    _Tp
    __beta_inc_cont_frac(_Tp __a, _Tp __b, _Tp __x)
    {
      constexpr unsigned int _S_itmax = 100;
      constexpr _Tp _S_eps = __gnu_cxx::__math_constants<_Tp>::__eps;

      _Tp __apb = __a + __b;
      _Tp __ap1 = __a + _Tp{1};
      _Tp __am1 = __a - _Tp{1};
      _Tp __c = _Tp{1};
      _Tp __d = _Tp{1} - __apb * __x / __ap1;
      if (std::abs(__d) < _S_eps)
	__d = _S_eps;
      __d = _Tp{1} / __d;
      _Tp __h = __d;
      for (unsigned int __m = 1; __m <= _S_itmax; ++__m)
	{
          _Tp __m2(2 * __m);

          _Tp __aa = _Tp(__m) * (__b - _Tp(__m)) * __x
                   / ((__am1 + __m2) * (__a + __m2));
          __d = _Tp{1} + __aa * __d;
          if (std::abs(__d) < _S_eps)
            __d = _S_eps;
          __c = _Tp{1} + __aa / __c;
          if (std::abs(__c) < _S_eps)
            __c = _S_eps;
          __d = _Tp{1} / __d;
          __h *= __d * __c;

          __aa = -(__a + _Tp(__m)) * (__apb + __m) * __x
               / ((__a + __m2) * (__ap1 + __m2));
          __d = _Tp{1} + __aa * __d;
          if (std::abs(__d) < _S_eps)
            __d = _S_eps;
          __c = _Tp{1} + __aa / __c;
          if (std::abs(__c) < _S_eps)
            __c = _S_eps;
          __d = _Tp{1} / __d;
          _Tp __del = __d * __c;
          __h *= __del;

          if (std::abs(__del - _Tp{1}) < _S_eps)
            return __h;
	}
      throw std::runtime_error("__beta_inc_cont_frac: "
			       "continued fractions failed to converge");
    }


  template<typename _Tp>
    _Tp
    __beta_inc(_Tp __a, _Tp __b, _Tp __x)
    {
      if (__isnan(__x) || __isnan(__a) || __isnan(__b))
	return __gnu_cxx::__math_constants<_Tp>::__NaN;

      if (__x < _Tp{0} || __x > _Tp{1})
	throw std::domain_error("__beta_inc: x out of range");

      _Tp __fact;
      if (__x == _Tp{0} || __x == _Tp{1})
	__fact = _Tp{0};
      else
	__fact = std::exp(std::lgamma(__a + __b)
                	- std::lgamma(__a) - std::lgamma(__b)
                	+ __a * std::log(__x) + __b * std::log(_Tp{1} - __x));

      if (__x < (__a + _Tp{1}) / (__a + __b + _Tp{2}))
	return __fact * __beta_inc_cont_frac(__a, __b, __x) / __a;
      else
	return _Tp{1} - __fact * __beta_inc_cont_frac(__b, __a, _Tp{1} - __x) / __b;
    }

_GLIBCXX_END_NAMESPACE_VERSION
} // namespace __detail
} // namespace std

#endif // __GLIBCXX_BITS_SF_BETA_TCC
