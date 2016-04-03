// Special functions -*- C++ -*-

// Copyright (C) 2006-2016 Free Software Foundation, Inc.
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
// (1) Handbook of Mathematical Functions,
//     ed. Milton Abramowitz and Irene A. Stegun,
//     Dover Publications,
//     Section 6, pp. 253-266
// (2) The Gnu Scientific Library, http://www.gnu.org/software/gsl
// (3) Numerical Recipes in C, by W. H. Press, S. A. Teukolsky,
//     W. T. Vetterling, B. P. Flannery, Cambridge University Press (1992),
//     2nd ed, pp. 213-216
// (4) Gamma, Exploring Euler's Constant, Julian Havil,
//     Princeton, 2003.

#ifndef _GLIBCXX_BITS_SF_BETA_TCC
#define _GLIBCXX_BITS_SF_BETA_TCC 1

#include <ext/math_const.h>

namespace std _GLIBCXX_VISIBILITY(default)
{
// Implementation-space details.
namespace __detail
{
_GLIBCXX_BEGIN_NAMESPACE_VERSION

  /**
   * @brief  Return the beta function: @f$ B(a,b) @f$.
   *
   * The beta function is defined by
   * @f[
   *   B(a,b) = \int_0^1 t^{a - 1} (1 - t)^{b - 1} dt
   *          = \frac{\Gamma(a)\Gamma(b)}{\Gamma(a+b)}
   * @f]
   *
   * @param __a The first argument of the beta function.
   * @param __b The second argument of the beta function.
   * @return  The beta function.
   */
  template<typename _Tp>
    _Tp
    __beta_gamma(_Tp __a, _Tp __b)
    {

      _Tp __bet;
      if (__a > __b)
	{
	  __bet = __gamma(__a) / __gamma(__a + __b);
	  __bet *= __gamma(__b);
	}
      else
	{
	  __bet = __gamma(__b) / __gamma(__a + __b);
	  __bet *= __gamma(__a);
	}

      return __bet;
    }

  /**
   * @brief  Return the beta function @f$B(a,b)@f$ using
   *	     the log gamma functions.
   *
   * The beta function is defined by
   * @f[
   *   B(a,b) = \int_0^1 t^{a - 1} (1 - t)^{b - 1} dt
   *          = \frac{\Gamma(a)\Gamma(b)}{\Gamma(a+b)}
   * @f]
   *
   * @param __a The first argument of the beta function.
   * @param __b The second argument of the beta function.
   * @return  The beta function.
   */
  template<typename _Tp>
    _Tp
    __beta_lgamma(_Tp __a, _Tp __b)
    {
      _Tp __bet = __log_gamma(__a)
		+ __log_gamma(__b)
		- __log_gamma(__a + __b);
      __bet = std::exp(__bet);
      return __bet;
    }


  /**
   * @brief  Return the beta function @f$B(x,y)@f$ using
   *	     the product form.
   *
   * The beta function is defined by
   * @f[
   *   B(a,b) = \int_0^1 t^{a - 1} (1 - t)^{b - 1} dt
   *          = \frac{\Gamma(a)\Gamma(b)}{\Gamma(a+b)}
   * @f]
   * Here, we employ the product form:
   * @f[
   *   B(a,b) = \frac{a + b}{a b} \prod_{k=1}^{\infty}
   *            \frac{1 + (a + b) / k}{(1 + a / k) (1 + b / k)}
   * @f]
   *
   * @param __a The first argument of the beta function.
   * @param __b The second argument of the beta function.
   * @return  The beta function.
   */
  template<typename _Tp>
    _Tp
    __beta_product(_Tp __a, _Tp __b)
    {

      _Tp __bet = (__a + __b) / (__a * __b);

      const unsigned int _S_max_iter = 1000000;

      for (unsigned int __k = 1; __k < _S_max_iter; ++__k)
	{
	  _Tp __term = (_Tp{1} + (__a + __b) / __k)
		     / ((_Tp{1} + __a / __k) * (_Tp{1} + __b / __k));
	  __bet *= __term;
	}

      return __bet;
    }


  /**
   * @brief  Return the beta function @f$ B(a,b) @f$.
   *
   * The beta function is defined by
   * @f[
   *   B(a,b) = \int_0^1 t^{a - 1} (1 - t)^{b - 1} dt
   *          = \frac{\Gamma(a)\Gamma(b)}{\Gamma(a+b)}
   * @f]
   *
   * @param __a The first argument of the beta function.
   * @param __b The second argument of the beta function.
   * @return  The beta function.
   */
  template<typename _Tp>
    _Tp
    __beta(_Tp __a, _Tp __b)
    {
      if (__isnan(__a) || __isnan(__b))
	return __gnu_cxx::__quiet_NaN<_Tp>();
      else
	return __beta_lgamma(__a, __b);
    }


  /**
   * Return the regularized incomplete beta function, @f$ I_x(a,b) @f$,
   * of arguments @c a, @c b, and @c x.
   *
   *
   * @param __a The first parameter
   * @param __b The second parameter
   * @param __x The argument
   */
  template<typename _Tp>
    _Tp
    __beta_inc_cont_frac(_Tp __a, _Tp __b, _Tp __x)
    {
      constexpr unsigned int _S_itmax = 100;
      constexpr _Tp _S_eps = __gnu_cxx::__epsilon<_Tp>();

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
      std::__throw_runtime_error("__beta_inc_cont_frac: "
				 "continued fractions failed to converge");
    }

  /**
   * Return the regularized incomplete beta function, @f$ I_x(a,b) @f$,
   * of arguments @c a, @c b, and @c x.
   *
   * The regularized incomplete beta function is defined by:
   * @f[
   *   I_x(a,b) = \frac{B_x(a,b)}{B(a,b)}
   * @f]
   * where
   * @f[
   *   B_x(a,b) = \int_0^x t^{a - 1} (1 - t)^{b - 1} dt
   * @f]
   * is the non-regularized beta function and @f$ B(a,b) @f$
   * is the usual beta function.
   *
   * @param __a The first parameter
   * @param __b The second parameter
   * @param __x The argument
   */
  template<typename _Tp>
    _Tp
    __beta_inc(_Tp __a, _Tp __b, _Tp __x)
    {
      constexpr auto _S_NaN = __gnu_cxx::__quiet_NaN<_Tp>();

      if (__isnan(__x) || __isnan(__a) || __isnan(__b))
	return _S_NaN;

      if (__x < _Tp{0} || __x > _Tp{1})
	std::__throw_domain_error("__beta_inc: argument out of range");
      else if (__x == _Tp{0} || __x == _Tp{1})
	return _Tp{0};
      else
	{
	  auto __fact = std::exp(std::lgamma(__a + __b)
		      - std::lgamma(__a) - std::lgamma(__b)
		      + __a * std::log(__x) + __b * std::log(_Tp{1} - __x));

	  if (__x < (__a + _Tp{1}) / (__a + __b + _Tp{2}))
	    return __fact * __beta_inc_cont_frac(__a, __b, __x) / __a;
	  else
	    return _Tp{1}
		 - __fact * __beta_inc_cont_frac(__b, __a, _Tp{1} - __x) / __b;
	}
    }

  /**
   * @brief  Return the Students T probability function.
   *
   * The students T propability function is related to the incomplete beta function:
   * @f[
   *   A(t|\nu) = 1 - I_{\frac{\nu}{\nu + t^2}}(\frac{\nu}{2}, \frac{1}{2})
   *   A(t|\nu) = 
   * @f]
   *
   * @param __t 
   * @param __nu 
   */
  template<typename _Tp>
    _GLIBCXX14_CONSTEXPR _Tp
    __students_t_cdf(_Tp __t, unsigned int __nu)
    {
      if (__isnan(__t))
	return std::numeric_limits<_Tp>::quiet_NaN();
      else
	return __beta_inc(_Tp{0.5L}, _Tp(__nu) / _Tp{2},
			  __t * __t / (_Tp(__nu) + __t * __t));
    }

  /**
   * @brief  Return the complement of the Students T probability function.
   *
   * The complement of the students T propability function is:
   * @f[
   *   A_c(t|\nu) = I_{\frac{\nu}{\nu + t^2}}(\frac{\nu}{2}, \frac{1}{2})
   * 		  = 1 - A(t|\nu)
   * @f]
   *
   * @param __t 
   * @param __nu 
   */
  template<typename _Tp>
    _GLIBCXX14_CONSTEXPR _Tp
    __students_t_cdfc(_Tp __t, unsigned int __nu)
    {
      if (__isnan(__t))
	return std::numeric_limits<_Tp>::quiet_NaN();
      else
	return __beta_inc(_Tp(__nu) / _Tp{2}, _Tp{0.5L},
			  _Tp(__nu) / (_Tp(__nu) + __t * __t));
    }

  /**
   * @brief  Return the F-distribution propability function.
   * This returns the probability that the observed chi-square for a correct model
   * exceeds the value @f$ \chi^2 @f$.
   *
   * The f-distribution propability function is related to the incomplete beta function:
   * @f[
   *   Q(F|\nu_1, \nu_2) = I_{\frac{\nu_2}{\nu_2 + \nu_1 F}}
   * 			     (\frac{\nu_2}{2}, \frac{\nu_1}{2})
   * @f]
   *
   * @param __nu1 
   * @param __nu2 
   */
  template<typename _Tp>
    _GLIBCXX14_CONSTEXPR _Tp
    __f_cdf(_Tp __F, unsigned int __nu1, unsigned int __nu2)
    {
      if (__isnan(__F))
	return std::numeric_limits<_Tp>::quiet_NaN();
      else if (__F < _Tp{0})
	std::__throw_domain_error(__N("__f_cdf: F is negative"));
      else
	return __beta_inc(_Tp(__nu2) / _Tp{2}, _Tp(__nu1) / _Tp{2},
			  _Tp(__nu2) / (_Tp(__nu2) + __nu1 * __F));
    }

  /**
   * @brief  Return the F-distribution propability function.
   * This returns the probability that the observed chi-square for a correct model
   * exceeds the value @f$ \chi^2 @f$.
   *
   * The f-distribution propability function is related to the incomplete beta function:
   * @f[
   *   P(F|\nu_1, \nu_2) = 1 - I_{\frac{\nu_2}{\nu_2 + \nu_1 F}}
   * 			     (\frac{\nu_2}{2}, \frac{\nu_1}{2})
   * 			 = 1 - Q(F|\nu_1, \nu_2)
   * @f]
   *
   * @param __F 
   * @param __nu1 
   * @param __nu2 
   */
  template<typename _Tp>
    _GLIBCXX14_CONSTEXPR _Tp
    __f_cdfc(_Tp __F, unsigned int __nu1, unsigned int __nu2)
    {
      if (__isnan(__F))
	return std::numeric_limits<_Tp>::quiet_NaN();
      else if (__F < _Tp{0})
	std::__throw_domain_error(__N("__f_cdfc: F is negative"));
      else
	return __beta_inc(_Tp(__nu1) / _Tp{2}, _Tp(__nu2) / _Tp{2},
			  __nu1 * __F / (_Tp(__nu2) + __nu1 * __F));
    }

  /**
   * @brief  Return the binomial cumulative distribution function.
   *
   * The binomial cumulative distribution function is related
   * to the incomplete beta function:
   * @f[
   *   P(p|n, k) = I_p(k, n-k+1)
   * @f]
   *
   * @param __p 
   * @param __n 
   * @param __k 
   */
  template<typename _Tp>
    _GLIBCXX14_CONSTEXPR _Tp
    __binomial_cdf(_Tp __p, unsigned int __n, unsigned int __k)
    {
      if (__isnan(__p))
	return std::numeric_limits<_Tp>::quiet_NaN();
      else if (__p < _Tp{0} || __p > _Tp{1})
	std::__throw_domain_error(__N("__binomial_cdf: "
				      "probability is out of range"));
      else if (__k == 0)
	return _Tp{1};
      else if (__k > __n)
	return _Tp{0};
      else
	return __beta_inc(_Tp(__k), _Tp(__n - __k - 1), __p);
    }

  /**
   * @brief  Return the complementary binomial cumulative distribution function.
   *
   * The binomial cumulative distribution function is related
   * to the incomplete beta function:
   * @f[
   *   Q(p|n, k) = I_{1-p}(n-k+1, k)
   * @f]
   *
   * @param __p 
   * @param __n 
   * @param __k 
   */
  template<typename _Tp>
    _GLIBCXX14_CONSTEXPR _Tp
    __binomial_cdfc(_Tp __p, unsigned int __n, unsigned int __k)
    {
      if (__isnan(__p))
	return std::numeric_limits<_Tp>::quiet_NaN();
      else if (__p < _Tp{0} || __p > _Tp{1})
	std::__throw_domain_error(__N("__binomial_cdfc: "
				      "probability is out of range"));
      else if (__k == 0)
	return _Tp{1};
      else if (__k > __n)
	return _Tp{0};
      else
	return __beta_inc(_Tp(__n - __k - 1), _Tp(__k), _Tp{1} - __p);
    }

_GLIBCXX_END_NAMESPACE_VERSION
} // namespace __detail
} // namespace std

#endif // __GLIBCXX_BITS_SF_BETA_TCC
