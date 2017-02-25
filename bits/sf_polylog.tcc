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

/** @file bits/sf_polylog.tcc
 *  This is an internal header file, included by other library headers.
 *  Do not attempt to use it directly. @headername{cmath}
 */

//
// ISO C++ 14882 TR29124: Mathematical Special Functions
//

// Written by Florian Goth and Edward Smith-Rowland.
//
// References:
// (1) David C. Wood, "The Computation of Polylogarithms."
//

#ifndef _GLIBCXX_BITS_SF_POLYLOG_TCC
#define _GLIBCXX_BITS_SF_POLYLOG_TCC 1

#pragma GCC system_header

#include <complex>
#include <ext/math_const.h>
#include <bits/complex_util.h>

namespace std _GLIBCXX_VISIBILITY(default)
{
// Implementation-space details.
namespace __detail
{
_GLIBCXX_BEGIN_NAMESPACE_VERSION


  template<typename _Tp>
    class _Terminator
    {
    private:

      using _Real = __num_traits_t<_Tp>;
      std::size_t _M_max_iter;
      _Real _M_toler;

    public:

      _Terminator(std::size_t __max_iter, _Real __mul = _Real{1})
      : _M_max_iter(__max_iter),
	_M_toler(__mul * std::numeric_limits<_Real>::epsilon())
      { }

      bool
      operator()(_Tp __term, _Tp __sum)
      {
	return (--_M_max_iter == 0
		|| std::abs(__term) < _M_toler * std::abs(__sum));
      }
    };


  template<typename _Tp>
    std::complex<_Tp>
    __clamp_pi(std::complex<_Tp> __w)
    {
      const auto _S_pi = __gnu_cxx::__const_pi(std::real(__w));
      const auto _S_i2pi = std::complex<_Tp>{0, _Tp{2} * _S_pi};
      while (__w.imag() > _S_pi)
	__w -= _S_i2pi;
      while (__w.imag() <= -_S_pi)
	__w += _S_i2pi;
      return __w;
    }

  template<typename _Tp>
    std::complex<_Tp>
    __clamp_0_m2pi(std::complex<_Tp> __w)
    {
      const auto _S_2pi = __gnu_cxx::__const_2_pi(std::real(__w));
      const auto _S_i2pi = std::complex<_Tp>{0, _S_2pi};
      while (__w.imag() > _Tp{0})
	__w = std::complex<_Tp>(__w.real(), __w.imag() - _S_2pi);
      while (__w.imag() <= -_S_2pi)
	__w = std::complex<_Tp>(__w.real(), __w.imag() + _S_2pi);
      return __w;
    }

  /**
   * This function treats the cases of positive integer index s.
   *
   * @f[
   *   Li_s(e^w) = \sum_{k=0, k != s-1} \zeta(s-k) \frac{w^k}{k!}
   *             + (H_{s-1} - \log(-w)) \frac{w^{s-1}}{(s-1)!}
   * @f]
   * The radius of convergence is @f$ |w| < 2 \pi @f$.
   * Note that this series involves a @f$ \log(-x) @f$.
   * gcc and Mathematica differ in their implementation
   * of @f$ \log(e^{i \pi}) @f$:
   * gcc: @f$ \log(e^{+- i * \pi}) = +- i \pi @f$
   * whereas Mathematica doesn't preserve the sign in this case:
   * @f$ \log(e^{+- i\pi}) = +i \pi @f$
   *
   * @param __s the positive integer index s.
   * @param __w the argument w.
   * @return the value of the polylogarithm.
   */
  template<typename _Tp>
    std::complex<_Tp>
    __polylog_exp_pos(unsigned int __s, std::complex<_Tp> __w)
    {
      const auto __proto = std::real(__w);
      const auto _S_2pi = __gnu_cxx::__const_2_pi(__proto);
      const auto _S_pi = __gnu_cxx::__const_pi(__proto);
      const auto _S_pipio6 = __gnu_cxx::__const_pi_sqr_div_6(__proto);
      std::complex<_Tp> __res = __riemann_zeta<_Tp>(__s);
      auto __wk = __w;
      auto __fact = _Tp{1};
      auto __harmonicN = _Tp{1}; // HarmonicNumber_1
      for (unsigned int __k = 1; __k <= __s - 2; ++__k)
	{
	  __res += __wk * __fact * __riemann_zeta<_Tp>(__s - __k);
	  __wk *= __w;
	  const auto __temp = _Tp{1} / _Tp(1 + __k);
	  __fact *= __temp;
	  __harmonicN += __temp;
	}
      // harmonicN now contains H_{s-1}.
      // fac should be 1/(n-1)!
      __res += (__harmonicN - std::log(-__w)) * __wk * __fact;
      __wk *= __w;
      __fact /= __s;
      __res -= __wk * __fact / _Tp{2};
      __wk *= __w;
      // Now comes the remainder of the series.
      const auto __pref = __wk / _S_pi / _S_2pi;
      // Remainder of series.
      __fact *= _Tp{2} * _Tp{3} / _Tp(__s + 2) / _Tp(__s + 3);
      auto __w2 = -(__w / _S_2pi) * (__w / _S_2pi);
      auto __w2k = __w2;
      auto __rzarg = _Tp{2};
      __fact /= _Tp(__s + 1); // (1/(s+1)!)
      auto __sum = std::complex<_Tp>(_S_pipio6 * __fact); // Sum the zeroth order term.
      const unsigned int __maxit = 200;
      _Terminator<std::complex<_Tp>> __done(__maxit);
      while (true)
	{
	  __rzarg += _Tp{2};
	  const auto __rz = __riemann_zeta(__rzarg);
	  const auto __term = (__rz * __fact) * __w2k;
	  __sum += __term;
	  if (__done(__term, __sum))
	    break;
	  __w2k *= __w2;
	  __fact *= _Tp(__rzarg) / _Tp(__rzarg + __s)
		  * _Tp(__rzarg + 1) / _Tp(__rzarg + __s + 1);
	}
      __res -= __pref * __sum;
      return __res;
    }

  /**
   * This function treats the cases of positive integer index s for real w.
   *
   * This specialization is worthwhile to catch the differing behaviour
   * of log(x).
   * @f[
   *   Li_s(e^w) = \sum_{k=0, k != s-1} \zeta(s-k) \frac{w^k}{k!}
   *             + \left(H_{s-1} - \log(-w)\right) \frac{w^{s-1}}{(s-1)!}
   * @f]
   * The radius of convergence is @f$ |w| < 2 \pi @f$.
   * Note that this series involves a @f$ \log(-x) @f$.
   * gcc and Mathematica differ in their implementation
   * of @f$ \log(e^{i\pi}) @f$:
   * gcc: @f$ \log(e^{+- i\pi}) = +- i\pi @f$
   * whereas Mathematica doesn't preserve the sign in this case:
   * @f$ \log(e^{+- i\pi}) = +i\pi @f$
   *
   * @param __s the positive integer index.
   * @param __w the argument
   * @return the value of the Polylogarithm
   */
  template<typename _Tp>
    std::complex<_Tp>
    __polylog_exp_pos(unsigned int __s, _Tp __w)
    {
      const auto __proto = std::real(__w);
      const auto _S_2pi = __gnu_cxx::__const_2_pi(__proto);
      const auto _S_pi = __gnu_cxx::__const_pi(__proto);
      const auto _S_pipio6 = __gnu_cxx::__const_pi_sqr_div_6(__proto);
      auto __res = __riemann_zeta<_Tp>(__s);
      auto __wk = __w;
      auto __fact = _Tp{1};
      auto __harmonicN = _Tp{1}; // HarmonicNumber_1
      for (unsigned int __k = 1; __k <= __s - 2; ++__k)
	{
	  __res += __wk * __fact * __riemann_zeta<_Tp>(__s - __k);
	  __wk *= __w;
	  auto __temp = _Tp{1} / (_Tp{1} + __k);
	  __fact *= __temp;
	  __harmonicN += __temp;
	}
      // harmonicN now contains H_{s-1}
      // fact should be 1/(n-1)!
      auto __imagtemp = __fact * __wk
		    * (__harmonicN - std::log(std::complex<_Tp>(-__w, _Tp{0})));
      __res += std::real(__imagtemp);
      __wk *= __w;
      __fact /= __s;
      __res -= __wk * __fact / _Tp{2};
      __wk *= __w;
      // Now comes the remainder of the series.
      const auto __pref = __wk / _S_pi / _S_2pi;
      bool __terminate = false;
      __fact /= _Tp(__s + 1); // (1/(n+1)!)
      // Subtract the zeroth order term.
      auto __sum = _S_pipio6 * __fact;
      // Remainder of series.
      __fact *= _Tp{2} * _Tp{3} / _Tp(__s + 2) / _Tp(__s + 3);
      auto __w2 = -(__w / _S_2pi) * (__w / _S_2pi);
      auto __w2k = __w2;
      auto __rzarg = _Tp{2};
      const unsigned int __maxit = 200;
      _Terminator<_Tp> __done(__maxit);
      while (true)
	{
	  __rzarg += _Tp{2};
	  auto __rz = __riemann_zeta(__rzarg);
	  auto __term = __rz * __fact * __w2k;
	  __sum += __term;
	  if (__done(__term, __sum))
	    break;
	  __w2k *= __w2;
	  __fact *= _Tp(__rzarg) / _Tp(__rzarg + __s)
		  * _Tp(__rzarg + 1) / _Tp(__rzarg + __s + 1);
	}
      __res -= __pref * __sum;
      return std::complex<_Tp>(__res, std::imag(__imagtemp));
    }

  /**
   * This function treats the cases of negative real index s.
   * Theoretical convergence is present for @f$ |w| < 2\pi @f$.
   * We use an optimized version of
   * @f[
   *   Li_s(e^w) = \Gamma(1-s)(-w)^{s-1} + \frac{(2\pi)^{-s}}{\pi} A_p(w)
   * @f]
   * @f[
   *   A_p(w) = \sum_k \frac{\Gamma(1+k-s)}{k!}
   *         \sin\left(\frac{\pi}{2} (s-k)\right)
   *         \left(\frac{w}{2\pi}\right)^k \zeta(1+k-s)
   * @f]
   * @param __s  The negative real index
   * @param __w  The complex argument
   * @return  The value of the polylogarithm.
   */
  template<typename _Tp>
    std::complex<_Tp>
    __polylog_exp_neg(_Tp __s, std::complex<_Tp> __w)
    {
      const auto _S_i = std::complex<_Tp>{0, 1};
      const auto _S_2pi = __gnu_cxx::__const_2_pi(__s);
      const auto _S_pi_2 = __gnu_cxx::__const_pi_half(__s);
      // Basic general loop, but s is a negative quantity here
      // FIXME Large s makes problems.
      // The series should be rearrangeable so that we only need
      // the ratio Gamma(1-s)/(2 pi)^s
      auto __ls = __log_gamma(_Tp{1} - __s);
      auto __res = std::exp(__ls - (_Tp{1} - __s) * std::log(-__w));
      const auto __wup = __w / _S_2pi;
      auto __w2k = __wup;
      auto __pref = _Tp{2} * std::pow(_S_2pi, -(_Tp{1} - __s));
      // here we factor up the ratio of Gamma(1 - s + k)/k! .
      // This ratio should be well behaved even for large k in the series
      // afterwards
      // Note that we have a problem for large s
      // Since s is negative we evaluate the Gamma Function
      // on the positive real axis where it is real.
      auto __gam = std::exp(__ls);

      auto __phase = __polar_pi(_Tp{1}, __s / _Tp{2});
      auto __cp = std::real(__phase);
      auto __sp = std::imag(__phase);
      // Here we add the expression that would result from ignoring
      // the zeta function in the series.
      std::complex<_Tp> __expis(__cp, __sp);
      auto __p = _S_2pi - _S_i * __w;
      auto __q = _S_2pi + _S_i * __w;
      // This can be optimized for real values of w
      __res += _S_i * __gam * (std::conj(__expis) * std::pow(__p, __s - _Tp{1})
	     - __expis * std::pow(__q, __s - _Tp{1}));
      // The above expression is the result of
      // sum_k Gamma(1+k-s) /k! * sin(pi /2* (s-k)) * (w/2/pi)^k
      // Therefore we only need to sample values of zeta(n) on the real axis
      // that really differ from one
      std::complex<_Tp> __sum = __sp * __gam * __riemann_zeta_m_1(_Tp{1} - __s);
      unsigned int __j = 1;
      __gam *= (_Tp{1} - __s);
      constexpr unsigned int __maxit = 200;
      _Terminator<std::complex<_Tp>> __done(__maxit);
      while (true)
	{
	  auto __rzarg = _Tp(1 + __j) - __s;
	  auto __rz = __riemann_zeta_m_1(__rzarg);
	  _Tp __sine;
	  // Save repeated recalculation of the sines.
	  if (__j & 1)
	    { // odd
	      __sine = __cp;
	      if (!((__j - 1) / 2 & 1))
		__sine = -__sine;
	    }
	  else
	    { // even
	      __sine = __sp;
	      if ((__j / 2) & 1)
		__sine = -__sine;
	    }
	  auto __term =  __w2k * (__gam * __sine * __rz);
	  __w2k *= __wup;
	  ++__j;
	  __gam  *= __rzarg / _Tp(__j); // == 1/(j+1) we incremented j above.
	  __sum += __term;
	  if (__done(__term, __sum))
	    break;
	}
      __res += __pref * __sum;
      return __res;
    }

  /**
   * This function treats the cases of negative integer index @f$ s = -p @f$
   * which are even.
   *
   * In that case the sine occuring in the expansion occasionally
   * takes on the value zero.
   * We use that to provide an optimized series for p = 2n:
   *
   * @f[
   *   Li_{-p}(e^w) = \Gamma(1+p) (-w)^{-p-1} - A_p(w) - \sigma B_p(w)
   * @f]
   * with
   * @f[
   *   A_p(w) = 2 (2\pi)^{-p-1} \frac{p!}{(2\pi)^{-p/2}}
   *           \left(1 + \frac{w^2}{4\pi^2}\right)^{-(p+1)/2}
   *          \cos\left[(1 + p)ArcTan\left(\frac{2\pi}{w}\right)\right]
   * @f]
   * and 
   * @f[
   *   B_p(w) = - 2 (2\pi)^{-p-1} \sum_{k = 0}^{\infty} 
   *           \frac{(2k + 1 + p)!}{(2k + 1)!}
   *           (-1)^k \left(\frac{w}{2\pi}\right)^{2k+1} [\zeta(2 + 2k + p) - 1]
   * @f]
   * This is suitable for @f$ |w| < 2 \pi @f$
   * The original series is (This might be worthwhile if we use
   * the already present table of the Bernoullis)
   * @f[
   *   Li_{-p}(e^w) = p! (-w)^{-p-1}
   *      - \sigma (2\pi)^p /\pi \sum_{k = 0}^{\infty}
   *       \frac{(2k + 1 + p)!}{(2k + 1)!}
   *       (-1)^k \left(\frac{w}{2\pi}\right)^{2k+1} \zeta(2 + 2k + p)
   * @f]
   * Where @f$ p = 4k (\sigma = 1) @f$ or @f$ p = 4k + 2 (\sigma = -1) @f$.
   *
   * @param __p the integral index @f$ p = 4k @f$ or @f$ p = 4k + 2 @f$.
   * @param __w The complex argument w
   * @return the value of the Polylogarithm.
   */
  template<typename _Tp>
    std::complex<_Tp>
    __polylog_exp_neg_even(unsigned int __p, std::complex<_Tp> __w)
    {
      const auto _S_2pi = __gnu_cxx::__const_2_pi(std::real(__w));
      const auto _S_pi = __gnu_cxx::__const_pi(std::real(__w));
      const int __pp = 1 + __p;
      const int __sigma = __p % 4 == 0 ? +1 : -1;
      const auto __lnp = __log_factorial<_Tp>(__p);
      auto __res = std::exp(__lnp - _Tp(__pp) * std::log(-__w));
      auto __wup = __w / _S_2pi;
      auto __wq = -__wup * __wup;
      auto __pref = _Tp{2} * std::pow(_S_2pi, _Tp(-__pp));
      // Subtract the expression A_p(w)
      __pref *= __sigma;
      __res -= __pref
	     * std::exp(__lnp - _Tp{0.5L} * __pp * std::log(_Tp{1} - __wq))
	     * std::cos(_Tp(__pp) * std::atan(_Tp{1} / __wup));
      unsigned int __2k = 0;
      auto __gam = __factorial<_Tp>(1 + __p);
      std::complex<_Tp> __sum;
      constexpr unsigned int __maxit = 400;
      _Terminator<std::complex<_Tp>> __done(__maxit);
      while (true)
	{
	  auto __term = __gam * __wup
		      * __riemann_zeta_m_1<_Tp>(__2k + 2 + __p);
	  __sum += __term;
	  if (__done(__term, __sum))
	    break;
	  __gam *= _Tp(__2k + 2 + __p) / _Tp(__2k + 2)
		 * _Tp(__2k + 3 + __p) / _Tp(__2k + 3);
	  __2k += 2;
	  __wup *= __wq;
	}
      __res -= __pref * __sum;
      return __res;
    }

  /**
   * This function treats the cases of negative integer index @f$ s = -p @f$
   * which are odd.
   *
   * In that case the sine occuring in the expansion occasionally vanishes.
   * We use that to provide an optimized series for @f$ p = 1 + 2k @f$:
   * In the template parameter sigma we transport whether
   * @f$ p = 4k + 1 (\sigma = 1) @f$ or @f$ p = 4k + 3  (\sigma = -1) @f$.
   *
   * @f[
   *   Li_{-p}(e^w) = \Gamma(1+p) (-w)^{-p-1} + \sigma A_p(w) - \sigma B_p(w)
   * @f]
   * with
   * @f[
   *   A_p(w) = 2 (2\pi)^{-p-1} p!
   *          \left(1 + \frac{w^2}{4\pi^2}\right)^{-(p + 1)/2}
   *           \cos((1 + p) ArcTan(2 \pi / w))
   * @f]
   * and 
   * @f[
   *   B_p(w) = 2 (2\pi)^{-p-1}
   *            \sum_{k=0}^{\infty}\frac{(2k + p)!}{(2k)!}
   *      \left(\frac{-w^2}{4\pi^2}\right)^k [\zeta(2k + p + 1) - 1]
   * @f]
   * This is suitable for @f$ |w| < 2 \pi @f$.
   * The original series is (This might be worthwhile if we use
   * the already present table of the Bernoullis)
   * @f[
   *   Li_{-p}(e^w) = p! (-w)^{-p-1}
   *      - 2\sigma (2\pi)^{-p-1} \sum_{k = 0}^{\infty}
   *       \frac{(2k + p)!}{(2k)!}
   *       (-1)^k \left(\frac{w}{2\pi}\right)^{2k} \zeta(1 + 2k + p)
   * @f]
   *
   * @param __p the integral index @f$ p = 4k + 1 @f$ or @f$ p = 4k + 3 @f$.
   * @param __w The complex argument w.
   * @return The value of the Polylogarithm.
   */
  template<typename _Tp>
    std::complex<_Tp>
    __polylog_exp_neg_odd(unsigned int __p, std::complex<_Tp> __w)
    {
      const auto _S_2pi = __gnu_cxx::__const_2_pi(std::real(__w));
      const auto _S_pi = __gnu_cxx::__const_pi(std::real(__w));
      const int __pp = 1 + __p;
      const int __sigma = __p % 4 == 1 ? +1 : -1;
      const auto __lnp = __log_factorial<_Tp>(__p);
      auto __res = std::exp(__lnp - _Tp(__pp) * std::log(-__w));
      auto __wup = __w / _S_2pi;
      auto __wq = -__wup * __wup;
      auto __pref = _Tp{2} * std::pow(_S_2pi, _Tp(-__pp));
      __pref *= __sigma;
      // Add the expression A_p(w)
      __res += __pref
	     * std::exp(__lnp - _Tp{0.5L} * __pp * std::log(_Tp{1} - __wq))
	     * std::cos(_Tp(__pp) * std::atan(_Tp{1} / __wup));
      auto __gam = __factorial<_Tp>(__p);
      unsigned int __2k = 0;
      std::complex<_Tp> __sum;
      constexpr unsigned int __maxit = 400;
      _Terminator<std::complex<_Tp>> __done(__maxit);
      while (true)
	{
	  auto __term = __gam * __wup
		      * __riemann_zeta_m_1<_Tp>(__2k + 1 + __p);
	  __sum += __term;
	  if (__done(__term, __sum))
	    break;
	  __gam *= _Tp(__2k + 1 + __p) / _Tp(__2k + 1)
		 * _Tp(__2k + 2 + __p) / _Tp(__2k + 2);
	  __2k += 2;
	  __wup *= __wq;
	}
      __res -= __pref * __sum;
      return __res;
  }

  /**
   * This function treats the cases of negative integer index s
   * and branches accordingly
   *
   * @param __s the integer index s.
   * @param __w The Argument w
   * @return The value of the Polylogarithm evaluated by a suitable function.
   */
  template<typename _Tp>
    std::complex<_Tp>
    __polylog_exp_neg(int __s, std::complex<_Tp> __w)
    { // negative integer __s
      const auto __p = -__s;
      switch (__p % 2)
      {
      case 0:
	return __polylog_exp_neg_even<_Tp>(__p, __w);
      case 1:
	return __polylog_exp_neg_odd<_Tp>(__p, __w);
      }
    }

  /**
   * This function treats the cases of positive real index s.
   *
   * The defining series is
   * @f[
   *   Li_s(e^w) = A_s(w) + B_s(w) + \Gamma(1-s)(-w)^{s-1}
   * @f]
   * with
   * @f[
   *   A_s(w) = \sum_{k=0}^{m} \zeta(s-k)w^k/k!
   * @f]
   * @f[
   *   B_s(w) = \sum_{k=m+1}^{\infty} \sin(\pi/2(s-k))
   *             \Gamma(1-s+k)\zeta(1-s+k) (w/2/\pi)^k/k!
   * @f]
   *
   * @param __s the positive real index s.
   * @param __w The complex argument w.
   * @return the value of the polylogarithm.
   */
  template<typename _Tp>
    std::complex<_Tp>
    __polylog_exp_pos(_Tp __s, std::complex<_Tp> __w)
    { // positive s
      const auto _S_2pi = __gnu_cxx::__const_2_pi(__s);
      const auto _S_pi = __gnu_cxx::__const_pi(__s);
      const auto _S_pi_2 = __gnu_cxx::__const_pi_half(__s);
      std::complex<_Tp> __res = __riemann_zeta(__s);
      auto __wk = __w;
      auto __phase = __polar_pi(_Tp{1}, __s / _Tp{2});
      auto __cp = std::real(__phase);
      auto __sp = std::imag(__phase);
      // This is \Gamma(1-s)(-w)^{s-1}
      __res += _S_pi / (_Tp{2} * __sp * __cp)
	  * std::exp(-__log_gamma(__s) + (__s - _Tp{1}) * std::log(-__w));
      auto __fact = _Tp{1};
      const auto __m = static_cast<unsigned int>(std::floor(__s));
      for (unsigned int __k = 1; __k <= __m; ++__k)
	{
	  __res += __wk * __fact
		 * __riemann_zeta(__s - _Tp(__k));
	  __wk *= __w;
	  __fact /= _Tp(1 + __k);
	}
      // fac should now be 1/(m+1)!
      const auto __pref = _Tp{2} * std::pow(_S_2pi, __s - _Tp{1});
      // Now comes the remainder of the series
      constexpr unsigned int __maxit = 100;
      unsigned int __j = 0;
      bool __terminate = false;
      auto __wup = __w / _S_2pi;
      auto __wbark = std::pow(__wup, _Tp(__m + 1));
      // It is 1 < 2 - s + m < 2 => Gamma(2-s+m) will not overflow
      // Here we factor up the ratio of Gamma(1 - s + k) / k!.
      // This ratio should be well behaved even for large k
      auto __gam = __gamma(_Tp(2 + __m) - __s) * __fact;
      while (!__terminate)
	{
	  auto __idx = __m + 1 + __j;
	  auto __zetaarg = _Tp(1 + __idx) - __s;
	  auto __rz = __riemann_zeta(__zetaarg);
	  auto __sine = __cp;
	  if (__idx & 1) // Save the repeated calculation of the sines.
	    { // odd
	      __sine = __cp;
	      if (!((__idx - 1) / 2 & 1))
		__sine = -__sine;
	    }
	  else
	    { // even
	      __sine = __sp;
	      if ((__idx / 2) & 1)
		__sine = -__sine;
	    }
	  auto __term = __wbark * __sine * __gam * __rz;
	  __wbark *= __wup;
	  __gam *= __zetaarg / _Tp(1 + __idx);
	  ++__j;
	  __terminate = (__j > __maxit
			 || __gnu_cxx::__fp_is_zero(__pref * __term));
	  __res += __pref * __term;
	}
      return __res;
    }

  /**
   * This function implements the asymptotic series for the polylog.
   * It is given by
   * @f[
   *    2 \sum_{k=0}^{\infty} \zeta(2k) w^{s-2k}/\Gamma(s-2k+1)
   *       -i \pi w^{s-1}/\Gamma(s)
   * @f]
   * for @f$ Re(w) >> 1 @f$
   *
   * Don't check this against Mathematica 8.
   * For real u the imaginary part of the polylog is given by
   * @f$ Im(Li_s(e^u)) = - \pi u^{s-1}/\Gamma(s) @f$.
   * Check this relation for any benchmark that you use.
   *
   * @param __s the real index s.
   * @param __w the large complex argument w.
   * @return the value of the polylogarithm.
   */
  template<typename _Tp>
    std::complex<_Tp>
    __polylog_exp_asymp(_Tp __s, std::complex<_Tp> __w)
    { // asymptotic expansion
      const auto _S_pi = __gnu_cxx::__const_pi(__s);
      // wgamma = w^{s-1} / \Gamma(s)
      auto __wgamma = std::exp((__s - _Tp{1}) * std::log(__w)
		    - __log_gamma(__s));
      auto __res = std::complex<_Tp>(_Tp{0}, -_S_pi) * __wgamma;
      // wgamma = w^s / Gamma(s+1)
      __wgamma *= __w / __s;
      constexpr unsigned int __maxiter = 100;
      bool __terminate = false;
      // zeta(0) * w^s / Gamma(s + 1)
      std::complex<_Tp> __oldterm = -_Tp{0.5L} * __wgamma;
      __res += _Tp{2} * __oldterm;
      std::complex<_Tp> __term;
      auto __wq = _Tp{1} / (__w * __w);
      unsigned int __k = 1;
      while (!__terminate)
	{
	  __wgamma *= __wq * (__s + _Tp(1 - 2 * __k))
		    * (__s + _Tp(2 - 2 * __k));
	  __term = __riemann_zeta<_Tp>(2 * __k) * __wgamma;
	  if (std::abs(__term) > std::abs(__oldterm))
	    __terminate = true; // Failure of asymptotic expansion.
	  if (__k > __maxiter || __gnu_cxx::__fp_is_zero(_Tp{2} * __term))
	    __terminate = true; // Precision goal reached.
	  __res += _Tp{2} * __term;
	  __oldterm = __term;
	  ++__k;
	}
      return __res;
  }

  /**
   * Theoretical convergence for Re(w) < 0.
   *
   * Seems to beat the other expansions for @f$ Re(w) < -\pi/2 - \pi/5 @f$.
   * Note that this is an implementation of the basic series:
   * @f[
   *   Li_s(e^z) = \sum_{k=1} e^{kz} * k^{-s}
   * @f]
   *
   * @param __s is an arbitrary type, integral or float.
   * @param __w something with a negative real part.
   * @return the value of the polylogarithm.
   */
  template<typename _PowTp, typename _Tp>
    _Tp
    __polylog_exp_negative_real_part(_PowTp __s, _Tp __w)
    {
      auto __ew = std::exp(__w);
      const auto __up = __ew;
      auto __res = __ew;
      unsigned int __maxiter = 500;
      bool __terminate = false;
      unsigned int __k = 2;
      while (!__terminate)
	{
	  __ew *= __up;
	  _Tp __temp = std::pow(__k, __s); // This saves us a type conversion
	  auto __term = __ew / __temp;
	  __terminate = (__k > __maxiter || __gnu_cxx::__fp_is_zero(__term));
	  __res += __term;
	  ++__k;
	}
      return __res;
  }

  /**
   * Here s is a positive integer and the function descends
   * into the different kernels depending on w.
   *
   * @param __s a positive integer.
   * @param __w an arbitrary complex number.
   * @return The value of the polylogarithm.
   */
  template<typename _Tp>
    std::complex<_Tp>
    __polylog_exp_int_pos(unsigned int __s, std::complex<_Tp> __w)
    {
      const auto _S_2pi = __gnu_cxx::__const_2_pi(std::real(__w));
      const auto _S_pi = __gnu_cxx::__const_pi(std::real(__w));
      const auto _S_pi_2 = __gnu_cxx::__const_pi_half(std::real(__w));
      auto __rw = __w.real();
      auto __iw = __w.imag();
      if (__gnu_cxx::__fp_is_real(__w)
	  && __gnu_cxx::__fp_is_equal(std::remainder(__iw, _S_2pi), _Tp{0}))
	{
	  if (__s > 1)
	    return __riemann_zeta<_Tp>(__s);
	  else
	    return std::numeric_limits<_Tp>::infinity();
	}
      else if (0 == __s)
	{
	  auto __t = std::exp(__w);
	  return __t / (_Tp{1} - __t);
	}
      else if (1 == __s)
	return -std::log(_Tp{1} - std::exp(__w));
      else
	{
	  if (__rw < -(_S_pi_2 + _S_pi / _Tp{5})  )
	    // Choose the exponentially converging series
	    return __polylog_exp_negative_real_part(__s, __w);
	  else if (__rw < _Tp{6})
	    // The transition point chosen here, is quite arbitrary
	    // and needs more testing.
	    // The reductions of the imaginary part yield the same results
	    // as Mathematica.
	    // Necessary to improve the speed of convergence
	    return __polylog_exp_pos(__s, __clamp_pi(__w));
	  else
	    // Wikipedia says that this is required for Wood's formula.
	    // FIXME: The series should terminate after a finite number
	    // of terms.
	    return __polylog_exp_asymp(static_cast<_Tp>(__s),
				       __clamp_0_m2pi(__w));
	}
    }

  /**
   * Here s is a positive integer and the function descends
   * into the different kernels depending on w.
   *
   * @param __s a positive integer
   * @param __w an arbitrary real argument w
   * @return the value of the polylogarithm.
   */
  template<typename _Tp>
    std::complex<_Tp>
    __polylog_exp_int_pos(unsigned int __s, _Tp __w)
    {
      const auto _S_pi = __gnu_cxx::__const_pi(std::real(__w));
      const auto _S_pi_2 = __gnu_cxx::__const_pi_half(std::real(__w));
      if (__gnu_cxx::__fp_is_zero(__w))
	{
	  if (__s > 1)
	    return __riemann_zeta<_Tp>(__s);
	  else
	    return std::numeric_limits<_Tp>::infinity();
	}
      else if (0 == __s)
	{
	  auto __t = std::exp(__w);
	  return __t / (_Tp{1} - __t);
	}
      else if (1 == __s)
	return -std::log(_Tp{1} - std::exp(__w));
      else
	{
	  if (__w < -(_S_pi_2 + _S_pi / _Tp{5}))
	    // Choose the exponentially converging series
	    return __polylog_exp_negative_real_part(__s,
						    std::complex<_Tp>(__w));
	  else if (__w < _Tp{6})
	    // The transition point chosen here, is quite arbitrary
	    // and needs more testing.
	    return __polylog_exp_pos(__s, __w);
	  else
	    // FIXME: The series should terminate
	    // after a finite number of terms.
	    return __polylog_exp_asymp(static_cast<_Tp>(__s),
				       std::complex<_Tp>(__w));
	}
    }

  /**
   * This treats the case where s is a negative integer.
   *
   * @param __s a negative integer.
   * @param __w an arbitrary complex number
   * @return the value of the polylogarith,.
   */
  template<typename _Tp>
    std::complex<_Tp>
    __polylog_exp_int_neg(int __s, std::complex<_Tp> __w)
    {
      const auto _S_2pi = __gnu_cxx::__const_2_pi(std::real(__w));
      const auto _S_pi = __gnu_cxx::__const_pi(std::real(__w));
      const auto _S_pi_2 = __gnu_cxx::__const_pi_half(std::real(__w));
      if ((((-__s) & 1) == 0) && __gnu_cxx::__fp_is_imag(__w))
	{
	  // Now s is odd and w on the unit-circle.
	  auto __iw = imag(__w);  //get imaginary part
	  auto __rem = std::remainder(__iw, _S_2pi);
	  if (__gnu_cxx::__fp_is_equal(std::abs(__rem), _Tp{0.5L}))
	    // Due to: Li_{-n}(-1) + (-1)^n Li_{-n}(1/-1) = 0.
	    return _Tp{0};
	  else
	    // No asymptotic expansion available... check the reduction.
	    return __polylog_exp_neg(__s, std::complex<_Tp>(__w.real(), __rem));
	}
      else
	{
	  if (std::real(__w) < -(_S_pi_2 + _S_pi / _Tp{5})  )
	    // Choose the exponentially converging series
	    return __polylog_exp_negative_real_part(__s, __w);
	  else if (std::real(__w) < _Tp{6}) // Arbitrary transition point...
	    // The reductions of the imaginary part yield the same results
	    // as Mathematica.
	    // Necessary to improve the speed of convergence.
	    return __polylog_exp_neg(__s, __clamp_pi(__w));
	  else
	    // Wikipedia says that this clamping is required for Wood's formula.
	    // FIXME: The series should terminate
	    // after a finite number of terms.
	    return __polylog_exp_asymp(static_cast<_Tp>(__s),
				       __clamp_0_m2pi(__w));
	}
    }

  /**
   * This treats the case where s is a negative integer and w is a real.
   *
   * @param __s a negative integer.
   * @param __w the argument.
   * @return the value of the polylogarithm.
   */
  template<typename _Tp>
    std::complex<_Tp>
    __polylog_exp_int_neg(const int __s, _Tp __w)
    {
      const auto _S_pi = __gnu_cxx::__const_pi(__w);
      const auto _S_pi_2 = __gnu_cxx::__const_pi_half(__w);
      if (__w < -(_S_pi_2 + _S_pi / _Tp{5})) // Choose exp'ly converging series.
	return __polylog_exp_negative_real_part(__s, std::complex<_Tp>(__w));
      else if (__gnu_cxx::__fp_is_zero(__w))
	return std::numeric_limits<_Tp>::infinity();
      else if (__w < _Tp{6}) // Arbitrary transition point...
	return __polylog_exp_neg(__s, std::complex<_Tp>(__w));
      else
	// FIXME: The series should terminate
	// after a finite number of terms.
	return __polylog_exp_asymp(static_cast<_Tp>(__s),
				   std::complex<_Tp>(__w));
    }

  /**
   * Return the polylog where s is a positive real value
   * and for complex argument.
   *
   * @param __s A positive real number.
   * @param __w the complex argument.
   * @return The value of the polylogarithm.
   */
  template<typename _Tp>
    std::complex<_Tp>
    __polylog_exp_real_pos(_Tp __s, std::complex<_Tp> __w)
    {
      const auto _S_2pi = __gnu_cxx::__const_2_pi(__s);
      const auto _S_pi = __gnu_cxx::__const_pi(__s);
      const auto _S_pi_2 = __gnu_cxx::__const_pi_half(__s);
      auto __rw = __w.real();
      auto __iw = __w.imag();
      if (__gnu_cxx::__fp_is_real(__w)
	  && __gnu_cxx::__fp_is_zero(std::remainder(__iw, _S_2pi)))
	{
	  if (__s > _Tp{1})
	    return __riemann_zeta(__s);
	  else
	    return std::numeric_limits<_Tp>::infinity();
	}
      if (__rw < -(_S_pi_2 + _S_pi/_Tp{5})) // Choose exp'ly converging series.
	return __polylog_exp_negative_real_part(__s, __w);
      if (__rw < _Tp{6}) // arbitrary transition point
	// The reductions of the imaginary part yield the same results
	// as Mathematica then.
	// Branch cuts??
	return __polylog_exp_pos(__s, __clamp_pi(__w));
      else
	// Wikipedia says that this is required for Wood's formula
	return __polylog_exp_asymp(__s, __clamp_0_m2pi(__w));
    }

  /**
   * Return the polylog where s is a positive real value and the argument
   * is real.
   *
   * @param __s  A positive real number tht does not reduce to an integer.
   * @param __w  The real argument w.
   * @return  The value of the polylogarithm.
   */
  template<typename _Tp>
    std::complex<_Tp>
    __polylog_exp_real_pos(_Tp __s, _Tp __w)
    {
      const auto _S_pi = __gnu_cxx::__const_pi(__s);
      const auto _S_pi_2 = __gnu_cxx::__const_pi_half(__s);
      if (__gnu_cxx::__fp_is_zero(__w))
	{
	  if (__s > _Tp{1})
	    return __riemann_zeta(__s);
	  else
	    return std::numeric_limits<_Tp>::infinity();
	}
      if (__w < -(_S_pi_2 + _S_pi / _Tp{5})) // Choose exp'ly converging series.
	return __polylog_exp_negative_real_part(__s, __w);
      if (__w < _Tp{6}) // arbitrary transition point
	return __polylog_exp_pos(__s, std::complex<_Tp>(__w));
      else
	return __polylog_exp_asymp(__s, std::complex<_Tp>(__w));
  }

  /**
   * Return the polylog where s is a negative real value
   * and for complex argument.
   * Now we branch depending on the properties of w in the specific functions
   *
   * @param __s  A negative real value that does not reduce
   *             to a negative integer.
   * @param __w  The complex argument.
   * @return  The value of the polylogarithm.
   */
  template<typename _Tp>
    std::complex<_Tp>
    __polylog_exp_real_neg(_Tp __s, std::complex<_Tp> __w)
    {
      const auto _S_pi = __gnu_cxx::__const_pi(__s);
      const auto _S_pi_2 = __gnu_cxx::__const_pi_half(__s);
      auto __rw = __w.real();
      auto __iw = __w.imag();
      if (__rw < -(_S_pi_2 + _S_pi/_Tp{5})) // Choose exp'ly converging series.
	return __polylog_exp_negative_real_part(__s, __w);
      else if (__rw < 6) // arbitrary transition point
	// The reductions of the imaginary part yield the same results
	// as Mathematica then.
	// Necessary to improve the speed of convergence.
	// Branch cuts??
	return __polylog_exp_neg(__s, __clamp_pi(__w));
      else
	// Wikipedia says that this is required for Wood's formula
	return __polylog_exp_asymp(__s, __clamp_0_m2pi(__w));
    }

  /**
   * Return the polylog where s is a negative real value and for real argument.
   * Now we branch depending on the properties of w in the specific functions.
   *
   * @param __s  A negative real value.
   * @param __w  A real argument.
   * @return  The value of the polylogarithm.
   */
  template<typename _Tp>
    std::complex<_Tp>
    __polylog_exp_real_neg(_Tp __s, _Tp __w)
    {
      const auto _S_pi = __gnu_cxx::__const_pi(__s);
      const auto _S_pi_2 = __gnu_cxx::__const_pi_half(__s);
      if (__w < -(_S_pi_2 + _S_pi / _Tp{5})) // Choose exp'ly converging series.
	return __polylog_exp_negative_real_part(__s, std::complex<_Tp>(__w));
      else if (__w < _Tp{6}) // arbitrary transition point
	return __polylog_exp_neg(__s, std::complex<_Tp>(__w));
      else
	return __polylog_exp_asymp(__s, std::complex<_Tp>(__w));
    }

  /**
   * This is the frontend function which calculates @f$ Li_s(e^w) @f$
   * First we branch into different parts depending on the properties of s.
   * This function is the same irrespective of a real or complex w,
   * hence the template parameter ArgType.
   *
   * @note: I *really* wish we could return a variant<Tp, std::complex<Tp>>.
   *
   * @param __s  The real order.
   * @param __w  The real or complex argument.
   * @return  The real or complex value of Li_s(e^w).
   */
  template<typename _Tp, typename _ArgType>
    __gnu_cxx::__promote_fp_t<std::complex<_Tp>, _ArgType>
    __polylog_exp(_Tp __s, _ArgType __w)
    {
      if (__isnan(__s) || __isnan(__w))
	return std::numeric_limits<_Tp>::quiet_NaN();
      else if (__s > _Tp{25}) // Cutoff chosen by some testing on the real axis.
	return __polylog_exp_negative_real_part(__s, __w);
      else
	{
	  const auto __p = __gnu_cxx::__fp_is_integer(__s);
	  if (__p)
	    { // The order s is an integer.
	      if (__p() >= 0)
		return __polylog_exp_int_pos(__p(), __w);
	      else
		return __polylog_exp_int_neg(__p(), __w);
	    }
	  else
	    {
	      if (__s > _Tp{0})
		return __polylog_exp_real_pos(__s, __w);
	      else
		return __polylog_exp_real_neg(__s, __w);
	    }
	}
    }

  /**
   * Return the polylog Li_s(x) for two real arguments.
   *
   * @param __s  The real index.
   * @param __x  The real argument.
   * @return The complex value of the polylogarithm.
   */
  template<typename _Tp>
    _Tp
    __polylog(_Tp __s, _Tp __x)
    {
      if (__isnan(__s) || __isnan(__x))
	return __gnu_cxx::__quiet_NaN(__s);
      else if (__gnu_cxx::__fp_is_zero(__x))
	return _Tp{0};
      else
	{
	  const auto __n = __gnu_cxx::__fp_is_integer(__s);
	  if (__n && __n() == 1)
	    return -std::log(_Tp{1} - __x);
	  else if (__n && __n() == 0)
	    return __x / (_Tp{1} - __x);
	  else if (__gnu_cxx::__fp_is_equal(__x, _Tp{-1}))
	    // Prevent blowups caused by reflecting the branch point.
	    return std::real(__polylog_exp(__s, _Tp{0})
				* (std::pow(_Tp{2}, _Tp{1} - __s) - _Tp{1}));
	  else if (__x < _Tp{0})
	    { // Use the reflection formula to access negative values.
	      auto __y = std::log(-__x);
	      return std::real(__polylog_exp(__s, _Tp{2} * __y)
				    * std::pow(_Tp{2}, _Tp{1} - __s)
			     - __polylog_exp(__s, __y));
	    }
	  else
	    {
	      auto __y = std::log(__x);
	      return std::real(__polylog_exp(__s, __y));
	    }
	}
    }

  /**
   * Return the polylog in those cases where we can calculate it.
   *
   * @param __s  The real index.
   * @param __w  The complex argument.
   * @return  The complex value of the polylogarithm.
   */
  template<typename _Tp>
    std::complex<_Tp>
    __polylog(_Tp __s, std::complex<_Tp> __w)
    {
      if (__isnan(__s) || __isnan(__w))
	return __gnu_cxx::__quiet_NaN(__s);
      else if (__gnu_cxx::__fp_is_real(__w))
	return __polylog(__s, std::real(__w));
      else
	return __polylog_exp(__s, std::log(__w));
    }

  /**
   * Return the Hurwitz Zeta function for real s and complex a.
   * This uses Jonquiere's identity:
   * @f[
   *    \frac{(i2\pi)^s}{\Gamma(s)}\zeta(a,1-s) = 
   *          Li_s(e^{i2\pi a}) + (-1)^s Li_s(e^{-i2\pi a})
   * @f]
   * @param __s The real argument
   * @param __a The complex parameter
   */
  template<typename _Tp>
    std::complex<_Tp>
    __hurwitz_zeta_polylog(_Tp __s, std::complex<_Tp> __a)
    {
      using _Cmplx = std::complex<_Tp>;
      const auto _S_pi = __gnu_cxx::__const_pi(__s);
      const auto _S_2pi = __gnu_cxx::__const_2_pi(__s);
      const auto _S_i2pi = _Cmplx{0, _S_2pi};
      const auto _S_pi_2 = __gnu_cxx::__const_pi_half(__s);
      if ((__a.imag() >= _Tp{0}
		&& (__a.real() >= _Tp{0} && __a.real() <  _Tp{1}))
       || (__a.imag() <  _Tp{0}
		&& (__a.real() >  _Tp{0} && __a.real() <= _Tp{1})))
	{
	  _Tp __t = _Tp{1} - __s;
	  auto __lpe = __polylog_exp(__t, _S_i2pi * __a);
	  /// @todo This __hurwitz_zeta_polylog prefactor is prone to overflow.
	  /// positive integer orders s?
	  auto __thing = std::exp(_Cmplx(_Tp{0}, -_S_pi_2 * __t));
	  return __gamma(__t)
	       * std::pow(_S_2pi, -__t)
	       * (__lpe * __thing + std::conj(__lpe * __thing));
	}
      else
	std::__throw_domain_error(__N("__hurwitz_zeta_polylog: Bad argument"));
    }

  /**
   * Return the Dirichlet eta function.
   * Currently, s must be real (complex type but negligible imaginary part.)
   * Otherwise std::domain_error is thrown.
   *
   * @param __s  The complex (but on-real-axis) argument.
   * @return  The complex Dirichlet eta function.
   * @throw  std::domain_error if the argument has a significant imaginary part.
   */
  template<typename _Tp>
    std::complex<_Tp>
    __dirichlet_eta(std::complex<_Tp> __s)
    {
      if (__isnan(__s))
	return __gnu_cxx::__quiet_NaN(std::imag(__s));
      else if (__gnu_cxx::__fp_is_real(__s))
	return -__polylog(__s.real(), _Tp{-1});
      else
	std::__throw_domain_error(__N("__dirichlet_eta: Bad argument"));
    }

  /**
   *  Return the Dirichlet eta function for real argument.
   *
   *  @param __s  The real argument.
   *  @return  The Dirichlet eta function.
   */
  template<typename _Tp>
    _Tp
    __dirichlet_eta(_Tp __s)
    {
      if (__isnan(__s))
	return __gnu_cxx::__quiet_NaN(__s);
      else
	return -std::real(__polylog(__s, _Tp{-1}));
    }

  /**
   * Return the Dirichlet beta function.
   * Currently, s must be real (complex type but negligible imaginary part.)
   * Otherwise std::domain_error is thrown.
   *
   * @param __s  The complex (but on-real-axis) argument.
   * @return  The Dirichlet Beta function of real argument.
   * @throw  std::domain_error if the argument has a significant imaginary part.
   */
  template<typename _Tp>
    _Tp
    __dirichlet_beta(std::complex<_Tp> __s)
    {
      const auto _S_i = std::complex<_Tp>{0, 1};
      if (__isnan(__s))
	return __gnu_cxx::__quiet_NaN(std::imag(__s));
      else if (__gnu_cxx::__fp_is_real(__s))
	return std::imag(__polylog(__s.real(), _S_i));
      else
	std::__throw_domain_error(__N("__dirichlet_beta: Bad argument."));
    }

  /**
   * Return the Dirichlet beta function for real argument.
   *
   * @param __s  The real argument.
   * @return  The Dirichlet Beta function of real argument.
   */
  template<typename _Tp>
    _Tp
    __dirichlet_beta(_Tp __s)
    {
      const auto _S_i = std::complex<_Tp>{0, 1};
      if (__isnan(__s))
	return std::numeric_limits<_Tp>::quiet_NaN();
      else
	return std::imag(__polylog(__s, _S_i));
    }

  /**
   *  Return the Dirichlet lambda function for real argument.
   *
   *  @param __s  The real argument.
   *  @return  The Dirichlet lambda function.
   */
  template<typename _Tp>
    _Tp
    __dirichlet_lambda(_Tp __s)
    {
      if (__isnan(__s))
	return std::numeric_limits<_Tp>::quiet_NaN();
      else
	return (__riemann_zeta(__s) + __dirichlet_eta(__s)) / _Tp{2};
    }

  /**
   * Return Clausen's function of integer order m and complex argument @c w.
   * The notation and connection to polylog is from Wikipedia
   *
   * @param __m  The non-negative integral order.
   * @param __w  The complex argument.
   * @return  The complex Clausen function.
   */
  template<typename _Tp>
    std::complex<_Tp>
    __clausen(unsigned int __m, std::complex<_Tp> __w)
    {
      const auto _S_i = std::complex<_Tp>{0, 1};
      auto __ple = __polylog_exp(_Tp(__m), _S_i * __w);
      if (__isnan(__w))
	return std::numeric_limits<_Tp>::quiet_NaN();
      else if (__m == 0)
	std::__throw_domain_error(__N("__clausen: Non-positive order"));
      else if (__m & 1)
	return __ple;
      else
	return _S_i * std::conj(__ple);
    }

  /**
   * Return Clausen's function of integer order m and real argument w.
   * The notation and connection to polylog is from Wikipedia
   *
   * @param __m  The integer order m >= 1.
   * @param __w  The real argument.
   * @return  The Clausen function.
   */
  template<typename _Tp>
    _Tp
    __clausen(unsigned int __m, _Tp __w)
    {
      const auto _S_i = std::complex<_Tp>{0, 1};
      auto __ple = __polylog_exp(_Tp(__m), _S_i * __w);
      if (__isnan(__w))
	return std::numeric_limits<_Tp>::quiet_NaN();
      else if (__m == 0)
	std::__throw_domain_error(__N("__clausen: Non-positive order"));
      else if (__m & 1)
	return std::real(__ple);
      else
	return std::imag(__ple);
    }

  /**
   * Return Clausen's sine sum Sl_m for positive integer order m
   * and complex argument w.
   * @see https://en.wikipedia.org/wiki/Clausen_function
   *
   * @param __m  The integer order m >= 1.
   * @param __w  The complex argument.
   * @return  The Clausen sine sum Sl_m(w),
   */
  template<typename _Tp>
    _Tp
    __clausen_s(unsigned int __m, std::complex<_Tp> __w)
    {
      const auto _S_i = std::complex<_Tp>{0, 1};
      auto __ple = __polylog_exp(_Tp(__m), _S_i * __w);
      if (__isnan(__w))
	return std::numeric_limits<_Tp>::quiet_NaN();
      else if (__m == 0)
	std::__throw_domain_error(__N("__clausen_s: Non-positive order"));
      else if (__m & 1)
	return std::imag(__ple);
      else
	return std::real(__ple);
    }

  /**
   * Return Clausen's sine sum Sl_m for positive integer order m
   * and real argument w.
   * @see https://en.wikipedia.org/wiki/Clausen_function
   *
   * @param __m  The integer order m >= 1.
   * @param __w  The complex argument.
   * @return  The Clausen sine sum Sl_m(w),
   */
  template<typename _Tp>
    _Tp
    __clausen_s(unsigned int __m, _Tp __w)
    {
      const auto _S_i = std::complex<_Tp>{0, 1};
      auto __ple = __polylog_exp(_Tp(__m), _S_i * __w);
      if (__isnan(__w))
	return std::numeric_limits<_Tp>::quiet_NaN();
      else if (__m == 0)
	std::__throw_domain_error(__N("__clausen_s: Non-positive order"));
      else if (__m & 1)
	return std::imag(__ple);
      else
	return std::real(__ple);
    }

  /**
   * Return Clausen's cosine sum Cl_m for positive integer order m
   * and complex argument w.
   * @see https://en.wikipedia.org/wiki/Clausen_function
   *
   * @param __m  The integer order m >= 1.
   * @param __w  The real argument.
   * @return  The Clausen cosine sum Cl_m(w),
   */
  template<typename _Tp>
    _Tp
    __clausen_c(unsigned int __m, std::complex<_Tp> __w)
    {
      const auto _S_i = std::complex<_Tp>{0, 1};
      auto __ple = __polylog_exp(_Tp(__m), _S_i * __w);
      if (__isnan(__w))
	return std::numeric_limits<_Tp>::quiet_NaN();
      else if (__m == 0)
	std::__throw_domain_error(__N("__clausen_c: Non-positive order"));
      else if (__m & 1)
	return std::real(__ple);
      else
	return std::imag(__ple);
    }

  /**
   * Return Clausen's cosine sum Cl_m for positive integer order m
   * and real argument w.
   * @see https://en.wikipedia.org/wiki/Clausen_function
   *
   * @param __m  The integer order m >= 1.
   * @param __w  The real argument.
   * @return  The real Clausen cosine sum Cl_m(w),
   */
  template<typename _Tp>
    _Tp
    __clausen_c(unsigned int __m, _Tp __w)
    {
      const auto _S_i = std::complex<_Tp>{0, 1};
      auto __ple = __polylog_exp(_Tp(__m), _S_i * __w);
      if (__isnan(__w))
	return std::numeric_limits<_Tp>::quiet_NaN();
      else if (__m == 0)
	std::__throw_domain_error(__N("__clausen_c: Non-positive order"));
      else if (__m & 1)
	return std::real(__ple);
      else
	return std::imag(__ple);
    }

  /**
   * Return the Fermi-Dirac integral of integer or real order s
   * and real argument x.
   * @see https://en.wikipedia.org/wiki/Clausen_function
   * @see http://dlmf.nist.gov/25.12.16
   *
   * @f[
   *    F_s(x) = \frac{1}{\Gamma(s+1)}\int_0^\infty \frac{t^s}{e^{t-x} + 1}dt
   *           = -Li_{s+1}(-e^x)
   * @f]
   *
   * @param __s  The order s > -1.
   * @param __x  The real argument.
   * @return  The real Fermi-Dirac cosine sum F_s(x),
   */
  template<typename _Sp, typename _Tp>
    _Tp
    __fermi_dirac(_Sp __s, _Tp __x)
    {
      const auto _S_i = std::complex<_Tp>{0, 1};
      const auto _S_pi = __gnu_cxx::__const_pi(__s);
      if (__isnan(__s) || __isnan(__x))
	return std::numeric_limits<_Tp>::quiet_NaN();
      else if (__s <= _Sp{-1})
	std::__throw_domain_error(__N("__fermi_dirac: "
				      "Order must be greater than -1"));
      else
	return -std::real(__polylog_exp(__s + _Sp{1}, __x + _S_i * _S_pi));
    }

  /**
   * Return the Bose-Einstein integral of integer or real order s
   * and real argument x.
   * @see https://en.wikipedia.org/wiki/Clausen_function
   * @see http://dlmf.nist.gov/25.12.16
   *
   * @f[
   *    G_s(x) = \frac{1}{\Gamma(s+1)}\int_0^\infty \frac{t^s}{e^{t-x} - 1}dt
   *           = Li_{s+1}(e^x)
   * @f]
   *
   * @param __s  The order s >= 0.
   * @param __x  The real argument.
   * @return  The real Fermi-Dirac cosine sum G_s(x),
   */
  template<typename _Sp, typename _Tp>
    _Tp
    __bose_einstein(_Sp __s, _Tp __x)
    {
      if (__isnan(__s) || __isnan(__x))
	return std::numeric_limits<_Tp>::quiet_NaN();
      else if (__s <= _Sp{0} && __x < _Tp{0})
	std::__throw_domain_error(__N("__bose_einstein: "
				      "Order must be greater than -1"));
      else
	return std::real(__polylog_exp(__s + _Sp{1}, __x));
    }

_GLIBCXX_END_NAMESPACE_VERSION
} // namespace __detail
} // namespace std

#endif // _GLIBCXX_BITS_SF_POLYLOG_TCC
