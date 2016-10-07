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

/** @file bits/sf_distributions.tcc
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

#ifndef _GLIBCXX_BITS_SF_DISTRIBUTIONS_TCC
#define _GLIBCXX_BITS_SF_DISTRIBUTIONS_TCC 1

#pragma GCC system_header

#include <ext/math_const.h>

namespace std _GLIBCXX_VISIBILITY(default)
{
// Implementation-space details.
namespace __detail
{
_GLIBCXX_BEGIN_NAMESPACE_VERSION

  /**
   *  @brief  Return the chi-squared propability function.
   *  This returns the probability that the observed chi-squared for a correct model
   *  is less than the value @f$ \chi^2 @f$.
   *
   *  The chi-squared propability function is related
   *  to the normalized lower incomplete gamma function:
   *  @f[
   *    P(\chi^2|\nu) = \Gamma_P(\frac{\nu}{2}, \frac{\chi^2}{2})
   *  @f]
   */
  template<typename _Tp>
    _Tp
    __chi_squared_pdf(_Tp __chi2, unsigned int __nu)
    {
      if (__isnan(__chi2))
	return std::numeric_limits<_Tp>::quiet_NaN();
      else if (__chi2 < _Tp{0})
	std::__throw_domain_error(__N("__chi_squared_cdf: "
				      "chi-squared is negative"));
      else
	return __pgamma(_Tp(__nu) / _Tp{2}, __chi2 / _Tp{2});
    }

  /**
   *  @brief  Return the complementary chi-squared propability function.
   *  This returns the probability that the observed chi-squared for a correct model
   *  is greater than the value @f$ \chi^2 @f$.
   *
   *  The complementary chi-squared propability function is related
   *  to the normalized upper incomplete gamma function:
   *  @f[
   *    Q(\chi^2|\nu) = \Gamma_Q(\frac{\nu}{2}, \frac{\chi^2}{2})
   *  @f]
   */
  template<typename _Tp>
    _Tp
    __chi_squared_pdfc(_Tp __chi2, unsigned int __nu)
    {
      if (__isnan(__chi2) || __isnan(__nu))
	return std::numeric_limits<_Tp>::quiet_NaN();
      else if (__chi2 < _Tp{0})
	std::__throw_domain_error(__N("__chi_square_pdfc: "
				      "chi-squared is negative"));
      else
	return __qgamma(_Tp(__nu) / _Tp{2}, __chi2 / _Tp{2});
    }

  /**
   * @brief Return the gamma propability distribution function.
   *
   * The formula for the gamma probability density function is:
   * @f[
   *    \Gamma(x|\alpha,\beta) = \frac{1}{\beta\Gamma(\alpha)}
   *                             (x/\beta)^{\alpha - 1} e^{-x/\beta}
   * @f]
   */
  template<typename _Tp>
    _Tp
    __gamma_pdf(_Tp __alpha, _Tp __beta, _Tp __x)
    {
      if (__isnan(__alpha) || __isnan(__beta) || __isnan(__x))
	return std::numeric_limits<_Tp>::quiet_NaN();
      else
	return std::pow(__beta, __alpha) * std::pow(__x, __alpha - _Tp{1})
	     * std::exp(__beta * __x) / __gamma(__alpha);
    }

  /**
   * @brief Return the gamma cumulative propability distribution function.
   *
   * The formula for the gamma probability density function is:
   * @f[
   *    \Gamma(x|\alpha,\beta) = \frac{1}{\beta\Gamma(\alpha)}
   *                             (x/\beta)^{\alpha - 1} e^{-x/\beta}
   * @f]
   */
  template<typename _Tp>
    _Tp
    __gamma_cdf(_Tp __alpha, _Tp __beta, _Tp __x)
    {
      if (__isnan(__alpha) || __isnan(__beta) || __isnan(__x))
	return std::numeric_limits<_Tp>::quiet_NaN();
      else
	return __tgamma_lower(__alpha, __beta * __x)
	     / __gamma(__alpha);
    }

  /**
   * @brief Return the gamma complementary cumulative propability
   *        distribution function.
   *
   * The formula for the gamma probability density function is:
   * @f[
   *    \Gamma(x|\alpha,\beta) = \frac{1}{\beta\Gamma(\alpha)}
   *                            (x/\beta)^{\alpha - 1} e^{-x/\beta} 
   * @f]
   */
  template<typename _Tp>
    _Tp
    __gamma_cdfc(_Tp __alpha, _Tp __beta, _Tp __x)
    {
      if (__isnan(__alpha) || __isnan(__beta) || __isnan(__x))
	return std::numeric_limits<_Tp>::quiet_NaN();
      else
	return __tgamma(__alpha, __beta * __x)
	     / __gamma(__alpha);
    }


  /**
   * @brief Return the Rice probability density function.
   *
   * The formula for the Rice probability density function is
   * @f[
   *   p(x|\nu,\sigma) = \frac{x}{\sigma^2}
   *                     \exp\left(-\frac{x^2+\nu^2}{2\sigma^2}\right)
   *                     I_0\left(\frac{x \nu}{\sigma^2}\right)
   * @f]
   * where @f$I_0(x)@f$ is the modified Bessel function of the first kind
   * of order 0 and @f$\nu >= 0@f$ and @f$\sigma > 0@f$.
   */
  template<typename _Tp>
    _Tp
    __rice_pdf(_Tp __nu, _Tp __sigma, _Tp __x)
    {
      if (__isnan(__nu) || __isnan(__sigma))
	return std::numeric_limits<_Tp>::quiet_NaN();
      else
	{
	  auto __sigma2 = __sigma * __sigma;
	  return (__x / __sigma2)
               * std::exp(-(__x * __x + __nu * __nu) / (_Tp{2} * __sigma2))
               * __cyl_bessel_i(_Tp{0}, (__x * __nu) / (__sigma2));
	}
    }


  /**
   * @brief Return the normal probability density function.
   *
   * The formula for the normal probability density function is
   * @f[
   *   f(x|\mu,\sigma) = \frac{e^{(x-\mu)^2/2\sigma^2}}{\sigma\sqrt{2\pi}}
   * @f]
   */
  template<typename _Tp>
    _Tp
    __normal_pdf(_Tp __nu, _Tp __sigma, _Tp __x)
    {
      constexpr auto _S_sqrt_2 = __gnu_cxx::__math_constants<_Tp>::__root_2;
      constexpr auto _S_sqrt_pi = __gnu_cxx::__math_constants<_Tp>::__root_pi;
      constexpr auto _S_sqrt_2pi = _S_sqrt_2 * _S_sqrt_pi;
      if (__isnan(__nu) || __isnan(__sigma))
	return std::numeric_limits<_Tp>::quiet_NaN();
      else
	{
	  __x -= __nu;
	  __x /= __sigma;
	  __x *= __x;
	  __x /= _Tp{2};
	  return std::exp(-(__x)) / (__sigma * _S_sqrt_2pi);
	}
    }

  /**
   * @brief Return the normal cumulative probability density function.
   *
   * The formula for the normal cumulative probability density function is
   * @f[
   *     F(x|\mu,\sigma)
   *        = \frac{1}{2}\left[ 1-erf(\frac{x-\mu}{\sqrt{2}\sigma}) \right]
   * @f]
   */
  template<typename _Tp>
    _Tp
    __normal_cdf(_Tp __mu, _Tp __sigma, _Tp __x)
    {
      constexpr auto _S_sqrt_2 = __gnu_cxx::__math_constants<_Tp>::__root_2;
      if (__isnan(__mu) || __isnan(__sigma) || __isnan(__x))
	return std::numeric_limits<_Tp>::quiet_NaN();
      else
	return _Tp{0.5L}
	     * (_Tp{1} + std::erf((__x - __mu) / (__sigma * _S_sqrt_2)));
    }


  /**
   * @brief Return the lognormal probability density function.
   *
   * The formula for the lognormal probability density function is
   * @f[
   *   f(x|\mu,\sigma) = \frac{e^{(\ln{x}-\mu)^2/2\sigma^2}}{\sigma\sqrt{2\pi}}
   * @f]
   */
  template<typename _Tp>
    _Tp
    __lognormal_pdf(_Tp __nu, _Tp __sigma, _Tp __x)
    {
      constexpr auto _S_sqrt_2 = __gnu_cxx::__math_constants<_Tp>::__root_2;
      constexpr auto _S_sqrt_pi = __gnu_cxx::__math_constants<_Tp>::__root_pi;
      constexpr auto _S_sqrt_2pi = _S_sqrt_2 * _S_sqrt_pi;
      if (__isnan(__nu) || __isnan(__sigma))
	return std::numeric_limits<_Tp>::quiet_NaN();
      else
	{
	  __x -= __nu;
	  __x /= __sigma;
	  __x *= __x;
	  __x /= _Tp{2};
	  return std::exp(-(std::log(__x))) / (__sigma * _S_sqrt_2pi);
	}
    }

  /**
   * @brief Return the lognormal cumulative probability density function.
   *
   * The formula for the lognormal cumulative probability density function is
   * @f[
   *   F(x|\mu,\sigma)
   *     = \frac{1}{2}\left[ 1-erf(\frac{\ln{x}-\mu}{\sqrt{2}\sigma}) \right]
   * @f]
   */
  template<typename _Tp>
    _Tp
    __lognormal_cdf(_Tp __mu, _Tp __sigma, _Tp __x)
    {
      constexpr auto _S_sqrt_2 = __gnu_cxx::__math_constants<_Tp>::__root_2;
      if (__isnan(__mu) || __isnan(__sigma) || __isnan(__x))
	return std::numeric_limits<_Tp>::quiet_NaN();
      else
	return _Tp{0.5L} * (_Tp{1} + std::erf((std::log(__x) - __mu)
					    / (__sigma * _S_sqrt_2)));
    }


  /**
   * @brief Return the exponential probability density function.
   *
   * The formula for the exponential probability density function is
   * @f[
   *   f(x|\lambda) = \lambda e^{-\lambda x} \mbox{ for } x >= 0
   * @f]
   */
  template<typename _Tp>
    _Tp
    __exponential_pdf(_Tp __lambda, _Tp __x)
    {
      if (__isnan(__lambda) || __isnan(__x))
	return std::numeric_limits<_Tp>::quiet_NaN();
      else if (__x < _Tp{0})
	return _Tp{0};
      else
	return __lambda * std::exp(-__lambda * __x);
    }

  /**
   * @brief Return the exponential cumulative probability density function.
   *
   * The formula for the exponential cumulative probability density function is
   * @f[
   *   F(x|\lambda) = 1 - e^{-\lambda x} \mbox{ for } x >= 0
   * @f]
   */
  template<typename _Tp>
    _Tp
    __exponential_cdf(_Tp __lambda, _Tp __x)
    {
      if (__isnan(__lambda) || __isnan(__x))
	return std::numeric_limits<_Tp>::quiet_NaN();
      else if (__x < _Tp{0})
	return _Tp{0};
      else
	return _Tp{1} - std::exp(-__lambda * __x);
    }


  /**
   * @brief Return the Weibull probability density function.
   *
   * The formula for the Weibull probability density function is
   * @f[
   *   f(x | a, b) = \frac{a}{b}
   *                 \left(\frac{x}{b} \right)^{a-1}
   *                 \exp{-\left(\frac{x}{b}\right)^a}
   *                 \mbox{ for } x >= 0
   * @f]
   */
  template<typename _Tp>
    _Tp
    __weibull_pdf(_Tp __a, _Tp __b, _Tp __x)
    {
      if (__isnan(__a) || __isnan(__b) || __isnan(__x))
	return std::numeric_limits<_Tp>::quiet_NaN();
      else if (__x < _Tp{0})
	return _Tp{0};
      else
	return (__a / __b) * std::pow(__x / __b, __a - _Tp{1})
 			   * std::exp(-std::pow(__x / __b, __a));
    }

  /**
   * @brief Return the Weibull cumulative probability density function.
   *
   * The formula for the Weibull cumulative probability density function is
   * @f[
   *   F(x|\lambda) = 1 - e^{-(x / b)^a} \mbox{ for } x >= 0
   * @f]
   */
  template<typename _Tp>
    _Tp
    __weibull_cdf(_Tp __a, _Tp __b, _Tp __x)
    {
      if (__isnan(__a) || __isnan(__b) || __isnan(__x))
	return std::numeric_limits<_Tp>::quiet_NaN();
      else if (__x < _Tp{0})
	return _Tp{0};
      else
	return _Tp{1} - std::exp(-std::pow(__x / __b, __a));
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
    _Tp
    __student_t_cdf(_Tp __t, unsigned int __nu)
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
    _Tp
    __student_t_cdfc(_Tp __t, unsigned int __nu)
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
   * @param __nu1 The number of degrees of freedom of sample 1
   * @param __nu2 The number of degrees of freedom of sample 2
   * @param __F The F statistic
   */
  template<typename _Tp>
    _Tp
    __fisher_f_cdf(_Tp __F, unsigned int __nu1, unsigned int __nu2)
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
    _Tp
    __fisher_f_cdfc(_Tp __F, unsigned int __nu1, unsigned int __nu2)
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
   * @brief  Return the binomial probability mass function.
   *
   * The binomial cumulative distribution function is related
   * to the incomplete beta function:
   * @f[
   *   f(k|n,p) = \binom{n}{k}p^k(1-p)^{n-k}
   * @f]
   *
   * @param __p 
   * @param __n 
   * @param __k 
   */
  template<typename _Tp>
    _Tp
    __binomial_pdf(_Tp __p, unsigned int __n, unsigned int __k)
    {
      if (__isnan(__p))
	return std::numeric_limits<_Tp>::quiet_NaN();
      else if (__p < _Tp{0} || __p > _Tp{1})
	std::__throw_domain_error(__N("__binomial_cdf: "
				      "probability is out of range"));
      else if (__k > __n)
	return _Tp{0};
      else if (__n == 0)
	return _Tp{1};
      else if (__k == 0)
	return std::pow(_Tp{1} - __p, __n);
      else if (__k == __n)
	return std::pow(__p, __n);
      else
	return __bincoef<_Tp>(__n, __k)
	     * std::pow(__p, __k)
	     * std::pow(_Tp{1} - __p, __n - __k);
    }


  /**
   * @brief  Return the binomial cumulative distribution function.
   *
   * The binomial cumulative distribution function is related
   * to the incomplete beta function:
   * @f[
   *   P(k|n,p) = I_p(k, n-k+1)
   * @f]
   *
   * @param __p 
   * @param __n 
   * @param __k 
   */
  template<typename _Tp>
    _Tp
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
   *   Q(k|n,p) = I_{1-p}(n-k+1, k)
   * @f]
   *
   * @param __p 
   * @param __n 
   * @param __k 
   */
  template<typename _Tp>
    _Tp
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

#endif // _GLIBCXX_BITS_SF_DISTRIBUTIONS_TCC

