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

/** @file bits/sf_gamma.tcc
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

#ifndef _GLIBCXX_BITS_SF_GAMMA_TCC
#define _GLIBCXX_BITS_SF_GAMMA_TCC 1

#include <ext/math_const.h>

namespace std _GLIBCXX_VISIBILITY(default)
{
// Implementation-space details.
namespace __detail
{
_GLIBCXX_BEGIN_NAMESPACE_VERSION

  /**
   *  @brief This returns Bernoulli numbers from a table or by summation
   *         for larger values.
   *
   *  Recursion is unstable.
   *
   *  @param __n the order n of the Bernoulli number.
   *  @return  The Bernoulli number of order n.
   */
  template<typename _Tp>
    _GLIBCXX14_CONSTEXPR _Tp
    __bernoulli_series(unsigned int __n)
    {
      constexpr _Tp
      __num[28]
      {
	 _Tp{1UL},	                 -_Tp{1UL} / _Tp{2UL},
	 _Tp{1UL} / _Tp{6UL},             _Tp{0UL},
	-_Tp{1UL} / _Tp{30UL},            _Tp{0UL},
	 _Tp{1UL} / _Tp{42UL},            _Tp{0UL},
	-_Tp{1UL} / _Tp{30UL},            _Tp{0UL},
	 _Tp{5UL} / _Tp{66UL},            _Tp{0UL},
	-_Tp{691UL} / _Tp{2730UL},        _Tp{0UL},
	 _Tp{7UL} / _Tp{6UL},             _Tp{0UL},
	-_Tp{3617UL} / _Tp{510UL},        _Tp{0UL},
	 _Tp{43867UL} / _Tp{798UL},       _Tp{0UL},
	-_Tp{174611UL} / _Tp{330UL},      _Tp{0UL},
	 _Tp{854513UL} / _Tp{138UL},      _Tp{0UL},
	-_Tp{236364091UL} / _Tp{2730UL},  _Tp{0UL},
	 _Tp{8553103UL} / _Tp{6UL},       _Tp{0UL}
      };

      if (__n == 0)
	return _Tp{1};

      if (__n == 1)
	return -_Tp{1} / _Tp{2};

      // Take care of the rest of the odd ones.
      if (__n % 2 == 1)
	return _Tp{0};

      // Take care of some small evens that are painful for the series.
      if (__n < 28)
	return __num[__n];

      auto __fact = _Tp{1};
      if ((__n / 2) % 2 == 0)
	__fact *= -_Tp{1};
      for (unsigned int __k = 1; __k <= __n; ++__k)
	__fact *= __k / (_Tp{2} * __gnu_cxx::__math_constants<_Tp>::__pi);
      __fact *= _Tp{2};

      auto __sum = _Tp{0};
      for (unsigned int __i = 1; __i < 1000; ++__i)
	{
	  auto __term = std::pow(_Tp(__i), -_Tp(__n));
	  if (__term < __gnu_cxx::__epsilon<_Tp>())
	    break;
	  __sum += __term;
	}

      return __fact * __sum;
    }


  /**
   *   @brief This returns Bernoulli number \f$B_n\f$.
   *
   *   @param __n the order n of the Bernoulli number.
   *   @return  The Bernoulli number of order n.
   */
  template<typename _Tp>
    inline _GLIBCXX14_CONSTEXPR _Tp
    __bernoulli(int __n)
    { return __bernoulli_series<_Tp>(__n); }


  /**
   *  @brief Return \f$log(\Gamma(x))\f$ by asymptotic expansion
   *         with Bernoulli number coefficients.  This is like
   *         Sterling's approximation.
   *
   *  @param __x The argument of the log of the gamma function.
   *  @return  The logarithm of the gamma function.
   */
  template<typename _Tp>
    _GLIBCXX14_CONSTEXPR _Tp
    __log_gamma_bernoulli(_Tp __x)
    {
      auto __lg = (__x - _Tp{0.5L}) * std::log(__x) - __x
		+ _Tp{0.5L} * std::log(_Tp{2}
		* __gnu_cxx::__math_constants<_Tp>::__pi);

      const auto __xx = __x * __x;
      auto __help = _Tp{1} / __x;
      for ( unsigned int __i = 1; __i < 20; ++__i )
	{
	  const auto __2i = _Tp(2 * __i);
	  __help /= __2i * (__2i - _Tp{1}) * __xx;
	  __lg += __bernoulli<_Tp>(2 * __i) * __help;
	}

      return __lg;
    }


  /**
   *  @brief Return \f$log(\Gamma(x))\f$ by the Lanczos method.
   *         This method dominates all others on the positive axis I think.
   *
   *  @param __x The argument of the log of the gamma function.
   *  @return  The logarithm of the gamma function.
   */
  template<typename _Tp>
    _GLIBCXX14_CONSTEXPR _Tp
    __log_gamma_lanczos(_Tp __x)
    {
      constexpr int _S_num_lanczos_cheb_7 = 9;
      constexpr _Tp
      _S_lanczos_cheb_7[_S_num_lanczos_cheb_7]
      {
       _Tp{ 0.99999999999980993227684700473478L},
       _Tp{ 676.520368121885098567009190444019L},
       _Tp{-1259.13921672240287047156078755283L},
       _Tp{ 771.3234287776530788486528258894L},
       _Tp{-176.61502916214059906584551354L},
       _Tp{ 12.507343278686904814458936853L},
       _Tp{-0.13857109526572011689554707L},
       _Tp{ 9.984369578019570859563e-6L},
       _Tp{ 1.50563273514931155834e-7L}
      };
      constexpr auto _S_log_root_2pi
	  = _Tp{0.9189385332046727417803297364056176L};

      const auto __xm1 = __x - _Tp{1};
      auto __sum = _S_lanczos_cheb_7[0];
      for(unsigned int __k = 1; __k < _S_num_lanczos_cheb_7; ++__k)
	__sum += _S_lanczos_cheb_7[__k] / (__xm1 + __k);

      auto __term1 = (__xm1 + _Tp{0.5L})
		   * std::log((__xm1 + _Tp{7.5L})
		   / __gnu_cxx::__math_constants<_Tp>::__e);
      auto __term2 = _S_log_root_2pi + std::log(__sum);
      auto __result = __term1 + (__term2 - _Tp{7});

      return __result;
    }


  /**
   *  @brief Return \f$ log(|\Gamma(x)|) \f$.
   *         This will return values even for \f$ x < 0 \f$.
   *         To recover the sign of \f$ \Gamma(x) \f$ for
   *         any argument use @a __log_gamma_sign.
   *
   *  @param __x The argument of the log of the gamma function.
   *  @return  The logarithm of the gamma function.
   */
  template<typename _Tp>
    _Tp
    __log_gamma(_Tp __x)
    {
#if _GLIBCXX_USE_C99_MATH_TR1
	return std::lgamma(__x);
#else
      if (__x > _Tp{0.5L})
	return __log_gamma_lanczos(__x);
      else
	{
	  const auto __sin_fact = std::abs(
			std::sin(__gnu_cxx::__math_constants<_Tp>::__pi * __x));
	  if (__sin_fact == _Tp{0})
	    std::__throw_domain_error(__N("__log_gamma: "
					  "argument is nonpositive integer"));
	  return __gnu_cxx::__math_constants<_Tp>::__ln_pi
		     - std::log(__sin_fact)
		     - __log_gamma_lanczos(_Tp{1} - __x);
	}
#endif
    }


  /**
   *   @brief Return the sign of \f$ \Gamma(x) \f$.
   *          At nonpositive integers zero is returned.
   *
   *   @param __x The argument of the gamma function.
   *   @return  The sign of the gamma function.
   */
  template<typename _Tp>
    _Tp
    __log_gamma_sign(_Tp __x)
    {
      if (__x > _Tp{0})
	return _Tp{1};
      else
	{
	  const auto __sin_fact
		  = std::sin(__gnu_cxx::__math_constants<_Tp>::__pi * __x);
	  if (__sin_fact > _Tp{0})
	    return _Tp{1};
	  else if (__sin_fact < _Tp{0})
	    return -_Tp{1};
	  else
	    return _Tp{0};
	}
    }


  /**
   *  @brief Return the logarithm of the binomial coefficient.
   *  The binomial coefficient is given by:
   *  @f[
   *    \left( __n \over __k \right) = \frac{n!}{(n-k)! k!}
   *  @f]
   *
   *  @param __n The first argument of the binomial coefficient.
   *  @param __k The second argument of the binomial coefficient.
   *  @return  The logarithm of the binomial coefficient.
   */
  template<typename _Tp>
    _Tp
    __log_bincoef(unsigned int __n, unsigned int __k)
    {
      if (__k > __n)
	return __gnu_cxx::__quiet_NaN<_Tp>();
      else if (__k == 0 || __k == __n)
	return _Tp{0};
      else
	_Tp __coeff = std::lgamma(_Tp(1 + __n))
		    - std::lgamma(_Tp(1 + __k))
		    - std::lgamma(_Tp(1 + __n - __k));
    }


  /**
   *  @brief Return the binomial coefficient.
   *  The binomial coefficient is given by:
   *  @f[
   *    \left( __n \over __k \right) = \frac{n!}{(n-k)! k!}
   *  @f]
   *
   *  @param __n The first argument of the binomial coefficient.
   *  @param __k The second argument of the binomial coefficient.
   *  @return  The binomial coefficient.
   */
  template<typename _Tp>
    _Tp
    __bincoef(unsigned int __n, unsigned int __k)
    {
      // Max e exponent before overflow.
      constexpr auto __max_bincoeff
                      = std::numeric_limits<_Tp>::max_exponent10
                      * std::log(_Tp(10)) - _Tp(1);

      if (__k > __n)
	return _Tp{0};
      else if (__k == 0 || __k == __n)
	return _Tp{1};
      else
        {
	  const auto __log_coeff = __log_bincoef<_Tp>(__n, __k);
	  if (__log_coeff > __max_bincoeff)
	    return __gnu_cxx::__quiet_NaN<_Tp>();
	  else
	    return std::exp(__log_coeff);
	}
    }


  /**
   *  @brief Return \f$ \Gamma(x) \f$.
   *
   *  @param __x The argument of the gamma function.
   *  @return  The gamma function.
   */
  template<typename _Tp>
    inline _Tp
    __gamma(_Tp __x)
    {
#if _GLIBCXX_USE_C99_MATH_TR1
      return std::lgamma(__x);
#else
      return std::exp(__log_gamma(__x));
#endif
    }


  template<typename _Tp>
    std::pair<_Tp, _Tp>
    __gamma_series(_Tp __a, _Tp __x)
    {
      constexpr auto _S_eps = 3.0 * __gnu_cxx::__epsilon<_Tp>();
      const auto _S_itmax = 10 * (10 + std::sqrt(std::abs(__a)));

      _Tp __lngam = std::lgamma(__a);

      if (__x < _Tp{0})
	throw std::domain_error("gamma_series: argument less than 0");
      else if (__x == _Tp{0})
	return std::make_pair(_Tp{0}, __lngam);
      else
	{
	  _Tp __aa = __a;
	  _Tp __term, __sum;
	  __term = __sum = _Tp{1} / __a;
	  for (unsigned int __n = 1; __n <= _S_itmax; ++__n)
	    {
	      __aa += _Tp{1};
	      __term *= __x / __aa;
	      __sum += __term;
	      if (std::abs(__term) < _S_eps * std::abs(__sum))
		{
		  _Tp __gamser = std::exp(-__x + __a * std::log(__x) - __lngam) * __sum;
		  return std::make_pair(__gamser, __lngam);
		}
	    }
	  throw std::logic_error("__gamma_series: a too large, itmax too small in routine.");
	}
    }


  template<typename _Tp>
    std::pair<_Tp, _Tp>
    __gamma_cont_frac(_Tp __a, _Tp __x)
    {
      constexpr auto _S_fpmin = 3 * __gnu_cxx::__min<_Tp>();
      constexpr auto _S_eps = 3 * __gnu_cxx::__epsilon<_Tp>();
      const auto _S_itmax = 10 * (10 + std::sqrt(std::abs(__a)));

      auto __lngam = std::lgamma(__a);

      auto __b = __x + _Tp{1} - __a;
      auto __c = _Tp{1} / _S_fpmin;
      auto __d = _Tp{1} / __b;
      auto __h = __d;
      for (unsigned int __n = 1; __n <= _S_itmax; ++__n)
	{
	  auto __an = -_Tp{__n} * (_Tp{__n} - __a);
	  __b += _Tp{2};
	  __d = __an * __d + __b;
	  if (std::abs(__d) < _S_fpmin)
	    __d = _S_fpmin;
	  __c = __b + __an / __c;
	  if (std::abs(__c) < _S_fpmin)
	    __c = _S_fpmin;
	  __d = _Tp{1} / __d;
	  _Tp __del = __d * __c;
	  __h *= __del;
	  if (std::abs(__del - _Tp{1}) < _S_eps)
	    {
	      _Tp __gamcf = std::exp(-__x + __a * std::log(__x) - __lngam) * __h;
	      return std::make_pair(__gamcf, __lngam);
	    }
	}
      throw std::logic_error("__gamma_cont_fraction: a too large, itmax too small in routine.");
    }


  /**
   *  @brief  Return the regularized lower incomplete gamma function.
   *  The regularized lower incomplete gamma function is defined by
   *  @f[
   *    P(a,x) = \frac{\gamma(a,x)}{\Gamma(a)}
   *  @f]
   *  where @f$ \Gamma(a) @f$ is the gamma function and
   *  @f[
   *    \gamma(a,x) = \int_0^x e^{-t}t^{a-1}dt  (a > 0)
   *  @f]
   *  is the lower incomplete gamma function.
   */
  template<typename _Tp>
    _Tp
    __gamma_p(_Tp __a, _Tp __x)
    {
      if (__x < _Tp{0} || __a <= _Tp{0})
	throw std::domain_error("gamma_p: invalid arguments");

      if (__x < __a + _Tp{1})
	return __gamma_series(__a, __x).first;
      else
	return _Tp{1} - __gamma_cont_frac(__a, __x).first;
    }


  /**
   *  @brief  Return the regularized upper incomplete gamma function.
   *  The regularized upper incomplete gamma function is defined by
   *  @f[
   *    Q(a,x) = \frac{\Gamma(a,x)}{\Gamma(a)}
   *  @f]
   *  where @f$ \Gamma(a) @f$ is the gamma function and
   *  @f[
   *    \Gamma(a,x) = \int_x^\infty e^{-t}t^{a-1}dt  (a > 0)
   *  @f]
   *  is the upper incomplete gamma function.
   */
  template<typename _Tp>
    _Tp
    __gamma_q(_Tp __a, _Tp __x)
    {
      if (__x < _Tp{0} || __a <= _Tp{0})
	throw std::domain_error("__gamma_q: invalid arguments");

      if (__x < __a + _Tp{1})
	return _Tp{1} - __gamma_series(__a, __x).first;
      else
	return __gamma_cont_frac(__a, __x).first;
    }


  /**
   *  @brief  Return the lower incomplete gamma function.
   *  The lower incomplete gamma function is defined by
   *  @f[
   *    \gamma(a,x) = \int_0^x e^{-t}t^{a-1}dt  (a > 0)
   *  @f]
   */
  template<typename _Tp>
    _Tp
    __gamma_l(_Tp __a, _Tp __x)
    {
      if (__x < _Tp{0} || __a <= _Tp{0})
	throw std::domain_error("__gamma_l: invalid arguments");

      if (__x < __a + _Tp{1})
      {
	std::pair<_Tp, _Tp> __gp = __gamma_series(__a, __x);
	return std::exp(__gp.second) * __gp.first;
      }
      else
      {
	std::pair<_Tp, _Tp> __gp = __gamma_cont_frac(__a, __x);
	return std::exp(__gp.second) * (_Tp{1} - __gp.first);
      }
    }


  /**
   *  @brief  Return the upper incomplete gamma function.
   *  The lower incomplete gamma function is defined by
   *  @f[
   *    \Gamma(a,x) = \int_x^\infty e^{-t}t^{a-1}dt  (a > 0)
   *  @f]
   */
  template<typename _Tp>
    _Tp
    __gamma_u(_Tp __a, _Tp __x)
    {
      if (__x < _Tp{0} || __a <= _Tp{0})
	throw std::domain_error("__gamma_u: invalid arguments");

      if (__x < __a + _Tp{1})
      {
	std::pair<_Tp, _Tp> __gp = __gamma_series(__a, __x);
	return std::exp(__gp.second) * (_Tp{1} - __gp.first);
      }
      else
      {
	std::pair<_Tp, _Tp> __gp = __gamma_cont_frac(__a, __x);
	return std::exp(__gp.second) * __gp.first;
      }
    }


  /**
   *  @brief  Return the logarithm of the (upper) Pochhammer symbol
   *  or the rising factorial function.
   *  The Pochammer symbol is defined by
   *  @f[
   *    (a)_n = \prod_{k=0}^{n-1} (a + k), (a)_0 = 1
   *          = \Gamma(a + n) / \Gamma(n)
   *  @f]
   *  Thus this function returns
   *  @f[
   *    ln[(a)_n] = \Gamma(a + n) - \Gamma(n), ln[(a)_0] = 0
   *  @f]
   *  Many notations exist: @f[ a^{\overline{n}} @f],
   *   @f[ \left[ a \over n \right] @f], and others.
   */
  template<typename _Tp>
    _Tp
    __log_pochhammer_u(_Tp __n, _Tp __a)
    {
      if (__isnan(__n) || __isnan(__a))
	return __gnu_cxx::__quiet_NaN<_Tp>();
      else if (__n == _Tp{0})
	return _Tp{0};
      else
	return std::lgamma(__a + __n) - std::lgamma(__a);
    }


  /**
   *  @brief  Return the (upper) Pochhammer function
   *  or the rising factorial function.
   *  The Pochammer symbol is defined by
   *  @f[
   *    (a)_n = \prod_{k=0}^{n-1} (a + k), (a)_0 = 1
   *          = \Gamma(a + n) / \Gamma(n)
   *  @f]
   *  Many notations exist: @f[ a^{\overline{n}} @f],
   *   @f[ \left[ a \over n \right] @f], and others.
   */
  template<typename _Tp>
    _Tp
    __pochhammer_u(_Tp __n, _Tp __a)
    {
      constexpr auto __log10{2.3025850929940456840179914546843642L};
      if (__isnan(__n) || __isnan(__a))
	return __gnu_cxx::__quiet_NaN<_Tp>();
      else if (__n == _Tp{0})
	return _Tp{1};
      else
	{
          _Tp __logpoch = std::lgamma(__a + __n) - std::lgamma(__a);
          if (std::abs(__logpoch)
              > std::numeric_limits<_Tp>::max_digits10 * __log10)
            return __gnu_cxx::__infinity<_Tp>();
          else
            return std::exp(__logpoch);
	}
    }


  /**
   *  @brief  Return the logarithm of the lower Pochhammer symbol
   *  or the falling factorial function.
   *  The lower Pochammer symbol is defined by
   *  @f[
   *    (a)_n = \prod_{k=0}^{n-1} (a - k), (a)_0 = 1
   *          = \Gamma(a + 1) / \Gamma(a - n + 1)
   *  @f]
   *  In particular, $f[ (n)_n = n! $f].
   *  Thus this function returns
   *  @f[
   *    ln[(a)_n] = \Gamma(a + 1) - \Gamma(a - n + 1), ln[(a)_0] = 0
   *  @f]
   *  Many notations exist: @f[ a^{\underline{n}} @f],
   *   @f[ \left{ a \over n \right} @f], and others.
   */
  template<typename _Tp>
    _Tp
    __log_pochhammer_l(_Tp __n, _Tp __a)
    {
      if (__isnan(__n) || __isnan(__a))
	return __gnu_cxx::__quiet_NaN<_Tp>();
      else if (__n == _Tp{0})
	return _Tp{0};
      else
	return std::lgamma(__a + 1) - std::lgamma(__a - __n + 1);
    }


  /**
   *  @brief  Return the logarithm of the lower Pochhammer symbol
   *  or the falling factorial function.
   *  The lower Pochammer symbol is defined by
   *  @f[
   *    (a)_n = \prod_{k=0}^{n-1} (a - k), (a)_0 = 1
   *          = \Gamma(a + 1) / \Gamma(a - n + 1)
   *  @f]
   *  In particular, $f[ (n)_n = n! $f].
   */
  template<typename _Tp>
    _Tp
    __pochhammer_l(_Tp __n, _Tp __a)
    {
      constexpr auto __log10{2.3025850929940456840179914546843642L};
      if (__isnan(__n) || __isnan(__a))
	return __gnu_cxx::__quiet_NaN<_Tp>();
      else if (__n == _Tp{0})
	return _Tp{1};
      else
	{
          auto __logpoch = std::lgamma(__a + 1) - std::lgamma(__a - __n + 1);
          if (std::abs(__logpoch)
              > std::numeric_limits<_Tp>::max_digits10 * __log10)
            return __gnu_cxx::__infinity<_Tp>();
          else
            return std::exp(__logpoch);
	}
    }


  /**
   *  @brief  Return the digamma function by series expansion.
   *  The digamma or @f$ \psi(x) @f$ function is defined by
   *  @f[
   *    \psi(x) = \frac{\Gamma'(x)}{\Gamma(x)}
   *  @f]
   *
   *  The series is given by:
   *  @f[
   *    \psi(x) = -\gamma_E - \frac{1}{x}
   *    	 \sum_{k=1}^{\infty} \frac{x - 1}{(k + 1)(x + k)}
   *  @f]
   */
  template<typename _Tp>
    _Tp
    __psi_series(_Tp __x)
    {
      _Tp __sum = -__gnu_cxx::__math_constants<_Tp>::__gamma_e;
      const unsigned int _S_max_iter = 100000;
      for (unsigned int __k = 0; __k < _S_max_iter; ++__k)
	{
	  const auto __term = (__x - _Tp{1}) / ((__k + 1) * (__k + __x));
	  __sum += __term;
	  if (std::abs(__term) < __gnu_cxx::__epsilon<_Tp>())
	    break;
	}
      return __sum;
    }


  /**
   *  @brief  Return the digamma function for large argument.
   *  The digamma or @f$ \psi(x) @f$ function is defined by
   *  @f[
   *    \psi(x) = \frac{\Gamma'(x)}{\Gamma(x)}
   *  @f]
   *
   *  The asymptotic series is given by:
   *  @f[
   *    \psi(x) = \ln(x) - \frac{1}{2x}
   *    	- \sum_{n=1}^{\infty} \frac{B_{2n}}{2 n x^{2n}}
   *  @f]
   */
  template<typename _Tp>
    _Tp
    __psi_asymp(_Tp __x)
    {
      auto __sum = std::log(__x) - _Tp{0.5L} / __x;
      const auto __xx = __x * __x;
      auto __xp = __xx;
      const unsigned int __max_iter = 100;
      for (unsigned int __k = 1; __k < __max_iter; ++__k)
	{
	  const _Tp __term = __bernoulli<_Tp>(2 * __k) / (2 * __k * __xp);
	  __sum -= __term;
	  if (std::abs(__term / __sum) < __gnu_cxx::__epsilon<_Tp>())
	    break;
	  __xp *= __xx;
	}
      return __sum;
    }


  /**
   *  @brief  Return the digamma function.
   *  The digamma or @f$ \psi(x) @f$ function is defined by
   *  @f[
   *    \psi(x) = \frac{\Gamma'(x)}{\Gamma(x)}
   *  @f]
   *  For negative argument the reflection formula is used:
   *  @f[
   *    \psi(x) = \psi(1-x) - \pi \cot(\pi x)
   *  @f]
   */
  template<typename _Tp>
    _Tp
    __psi(_Tp __x)
    {
      constexpr auto _S_eps = _Tp{4} * __gnu_cxx::__epsilon<_Tp>();
      constexpr auto _S_x_asymp = _Tp{20};
      constexpr auto __gamma_E = __gnu_cxx::__math_constants<_Tp>::__gamma_e;
      constexpr auto __2_ln_2 = 2 * __gnu_cxx::__math_constants<_Tp>::__ln_2;

      const auto __n = std::nearbyint(__x);
      const bool __integral = (std::abs(__x - _Tp{__n}) < _S_eps);
      const auto __m = std::nearbyint(2 * __x);
      const bool __half_integral = !__integral
				&& (std::abs(2 * __x - _Tp{__m}) < _S_eps);
      if (__integral)
	{
	  if (__n <= 0)
	    return __gnu_cxx::__quiet_NaN<_Tp>();
	  else
	    {
	      _Tp __sum = -__gamma_E;
	      for (int __k = 1; __k < __n; ++__k)
		__sum += _Tp{1} / __k;
	      return __sum;
	    }
	}
      if (__half_integral)
	{
	  _Tp __sum = -__gamma_E - __2_ln_2;
	  for (int __k = 1; __k < __m / 2; ++__k)
	    __sum += _Tp{2} / (2 * __k - 1);
	  return __sum;
	}
      else if (__x < _Tp{0})
	{
	  constexpr auto __pi = __gnu_cxx::__math_constants<_Tp>::__pi;
	  return __psi(_Tp{1} - __x) - __pi / std::tan(__pi * __x);
	}
      else if (__x > _S_x_asymp)
	return __psi_asymp(__x);
      else
	{
	  //return __psi_series(__x);
	  // The series does not converge quickly enough.
	  // Reflect to larger argument and use asymptotic expansion.
	  auto __w = _Tp{0};
	  auto __y = __x;
	  while (__y <= _S_x_asymp)
	    {
	      __w += 1 / __y;
	      __y += 1;
	    }
	  return __psi_asymp(__y) - __w;
	}
    }


  /**
   *  @brief  Return the polygamma function @f$ \psi^{(n)}(x) @f$.
   *
   *  The polygamma function is related to the Hurwitz zeta function:
   *  @f[
   *    \psi^{(n)}(x) = (-1)^{n+1} m! \zeta(m+1,x)
   *  @f]
   */
  template<typename _Tp>
    _Tp
    __psi(unsigned int __n, _Tp __x)
    {
      if (__x <= _Tp{0})
	std::__throw_domain_error(__N("__psi: argument out of range"));
      else if (__n == 0)
	return __psi(__x);
      else
	{
	  const auto __hzeta = __hurwitz_zeta(_Tp{__n + 1}, __x);
	  const auto __ln_nfact = __log_gamma(_Tp{__n + 1});
	  auto __result = std::exp(__ln_nfact) * __hzeta;
	  if (__n % 2 == 1)
	    __result = -__result;
	  return __result;
	}
    }

_GLIBCXX_END_NAMESPACE_VERSION
} // namespace __detail
} // namespace std

#endif // _GLIBCXX_BITS_SF_GAMMA_TCC

