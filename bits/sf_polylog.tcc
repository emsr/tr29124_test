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


  /**
   * This class manages the termination of series.
   * Termination conditions involve both a maximum iteration count
   * and a relative precision.
   */
  template<typename _Tp>
    class _Terminator
    {
    private:

      using _Real = __num_traits_t<_Tp>;
      const std::size_t _M_max_iter;
      std::size_t _M_curr_iter;
      _Real _M_toler;

    public:

      _Terminator(std::size_t __max_iter, _Real __mul = _Real{1})
      : _M_max_iter(__max_iter), _M_curr_iter{0},
	_M_toler(std::abs(__mul) * std::numeric_limits<_Real>::epsilon())
      { }

      /// Return the current number of terms summed.
      std::size_t
      num_terms() const
      { return this->_M_curr_iter; }

      /// Detect if the sum should terminate either because the incoming term
      /// is small enough or the maximum number of terms has been reached.
      bool
      operator()(_Tp __term, _Tp __sum)
      {
	if (this->_M_curr_iter >= this->_M_max_iter
	    || ++this->_M_curr_iter == this->_M_max_iter)
	  return true;
	else if (std::abs(__term) < this->_M_toler * std::abs(__sum))
	  return true;
	else
	  return false;
      }
    };

  /**
   * This class manages the termination of asymptotic series.
   * In particular, this termination watches for the growth of the sequence
   * of terms to stop the series.
   *
   * Termination conditions involve both a maximum iteration count
   * and a relative precision.
   */
  template<typename _Tp>
    class _AsympTerminator
    {
    private:

      using _Real = __num_traits_t<_Tp>;
      const std::size_t _M_max_iter;
      std::size_t _M_curr_iter;
      _Real _M_toler;
      _Real _M_prev_term = std::numeric_limits<_Real>::max();
      bool _M_stop_asymp = false;

    public:

      _AsympTerminator(std::size_t __max_iter, _Real __mul = _Real{1})
      : _M_max_iter(__max_iter), _M_curr_iter{0},
	_M_toler(std::abs(__mul) * std::numeric_limits<_Real>::epsilon())
      { }

      /// Filter a term before applying it to the sum.
      _Tp
      operator<<(_Tp __term)
      {
	if (std::abs(__term) > this->_M_prev_term)
	  {
	    this->_M_stop_asymp = true;
	    return _Tp{0};
	  }
	else
	  return __term;
      }

      /// Return the current number of terms summed.
      std::size_t
      num_terms() const
      { return this->_M_curr_iter; }

      /// Detect if the sum should terminate either because the incoming term
      /// is small enough or because the terms are starting to grow or
      //  the maximum number of terms has been reached.
      bool
      operator()(_Tp __term, _Tp __sum)
      {
	if (this->_M_stop_asymp)
	  return true;
	else
	  {
	    const auto __aterm = std::abs(__term);
	    this->_M_stop_asymp = (__aterm > this->_M_prev_term);
	    this->_M_prev_term = __aterm;
	    if (this->_M_curr_iter >= this->_M_max_iter
	    || ++this->_M_curr_iter == this->_M_max_iter)
	      return true;
	    else if (__aterm < this->_M_toler * std::abs(__sum))
	      return true;
	    else if (this->_M_stop_asymp)
	      return true;
	    else
	      return false;
	  }
      }
    };

  template<typename _Tp>
    std::complex<_Tp>
    __clamp_pi(std::complex<_Tp> __z)
    {
      const auto _S_pi = __gnu_cxx::__const_pi(std::real(__z));
      const auto _S_i2pi = std::complex<_Tp>{0, _Tp{2} * _S_pi};
      while (__z.imag() > _S_pi)
	__z -= _S_i2pi;
      while (__z.imag() <= -_S_pi)
	__z += _S_i2pi;
      return __z;
    }

  template<typename _Tp>
    std::complex<_Tp>
    __clamp_0_m2pi(std::complex<_Tp> __z)
    {
      const auto _S_2pi = __gnu_cxx::__const_2_pi(std::real(__z));
      const auto _S_i2pi = std::complex<_Tp>{0, _S_2pi};
      while (__z.imag() > _Tp{0})
	__z = std::complex<_Tp>(__z.real(), __z.imag() - _S_2pi);
      while (__z.imag() <= -_S_2pi)
	__z = std::complex<_Tp>(__z.real(), __z.imag() + _S_2pi);
      return __z;
    }

  /**
   * This function treats the cases of positive integer index s
   * for complex argument w.
   *
   * @f[
   *   Li_s(e^w) = \sum_{k=0, k != s-1} \zeta(s-k) \frac{w^k}{k!}
   *             + \left[H_{s-1} - \log(-w)\right] \frac{w^{s-1}}{(s-1)!}
   * @f]
   * The radius of convergence is @f$ |w| < 2 \pi @f$.
   * Note that this series involves a @f$ \log(-x) @f$.
   * gcc and Mathematica differ in their implementation
   * of @f$ \log(e^{i\pi}) @f$:
   * gcc: @f$ \log(e^{+-i\pi}) = +-i\pi @f$
   * whereas Mathematica doesn't preserve the sign in this case:
   * @f$ \log(e^{+- i\pi}) = +i \pi @f$
   *
   * @param __s the positive integer index.
   * @param __w the argument.
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
	  __res += __fact * __riemann_zeta<_Tp>(__s - __k) * __wk;
	  __wk *= __w;
	  const auto __temp = _Tp{1} / _Tp(1 + __k);
	  __fact *= __temp;
	  __harmonicN += __temp;
	}
      // harmonicN now contains H_{s-1}.
      // fact should be 1/(s-1)!
      __res += (__harmonicN - std::log(-__w)) * __wk * __fact;
      __wk *= __w;
      __fact /= __s; // 1/s!
      __res -= __wk * __fact / _Tp{2};
      __wk *= __w;
      // Now comes the remainder of the series.
      const auto __pref = __wk / _S_pi / _S_2pi;
      __fact /= _Tp(__s + 1); // 1/(s+1)!
      // Subtract the zeroth order term.
      __res -= _S_pipio6 * __fact * __pref;
      __fact *= _Tp{2} / _Tp(__s + 2) * _Tp{3} / _Tp(__s + 3);
      const auto __wbar = __w / _S_2pi;
      const auto __w2 = -__wbar * __wbar;
      auto __w2k = __w2;
      auto __rzarg = _Tp{2};
      const unsigned int __maxit = 200;
      _Terminator<std::complex<_Tp>> __done(__maxit);
      while (true)
	{
	  __rzarg += _Tp{2};
	  const auto __rzeta = __riemann_zeta(__rzarg);
	  const auto __term = __pref * __fact * __rzeta * __w2k;
	  __res -= __term;
	  if (__done(__term, __res))
	    break;
	  __w2k *= __w2;
	  __fact *= _Tp(__rzarg) / _Tp(__s + __rzarg)
		  * _Tp(__rzarg + 1) / _Tp(__s + __rzarg + 1);
	}
      return __res;
    }

  /**
   * This function treats the cases of positive integer index s
   * for real argument w.
   *
   * This specialization is worthwhile to catch the differing behaviour
   * of log(x).
   * @f[
   *   Li_s(e^w) = \sum_{k=0, k != s-1} \zeta(s-k) \frac{w^k}{k!}
   *             + \left[H_{s-1} - \log(-w)\right] \frac{w^{s-1}}{(s-1)!}
   * @f]
   * The radius of convergence is @f$ |w| < 2 \pi @f$.
   * Note that this series involves a @f$ \log(-x) @f$.
   * gcc and Mathematica differ in their implementation
   * of @f$ \log(e^{i\pi}) @f$:
   * gcc: @f$ \log(e^{+-i\pi}) = +-i\pi @f$
   * whereas Mathematica doesn't preserve the sign in this case:
   * @f$ \log(e^{+- i\pi}) = +i\pi @f$
   *
   * @param __s the positive integer index.
   * @param __w the argument.
   * @return the value of the polylogarithm.
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
	  __res += __fact * __riemann_zeta<_Tp>(__s - __k) * __wk;
	  __wk *= __w;
	  const auto __temp = _Tp{1} / _Tp(1 + __k);
	  __fact *= __temp;
	  __harmonicN += __temp;
	}
      // harmonicN now contains H_{s-1}
      // fact should be 1/(s-1)!
      const auto __imagtemp = __fact * __wk
			    * (__harmonicN - std::log(std::complex<_Tp>(-__w)));
      __res += std::real(__imagtemp);
      __wk *= __w;
      __fact /= __s; // 1/s!
      __res -= __wk * __fact / _Tp{2};
      __wk *= __w;
      // Now comes the remainder of the series.
      const auto __pref = __wk / _S_pi / _S_2pi;
      __fact /= _Tp(__s + 1); // 1/(s+1)!
      // Subtract the zeroth order term.
      __res -= _S_pipio6 * __fact * __pref;
      __fact *= _Tp{2} / _Tp(__s + 2) * _Tp{3} / _Tp(__s + 3);
      const auto __wbar = __w / _S_2pi;
      const auto __w2 = -__wbar * __wbar;
      auto __w2k = __w2;
      auto __rzarg = _Tp{2};
      const unsigned int __maxit = 200;
      _Terminator<_Tp> __done(__maxit);
      while (true)
	{
	  __rzarg += _Tp{2};
	  const auto __rzeta = __riemann_zeta(__rzarg);
	  const auto __term = __pref * __fact * __rzeta * __w2k;
	  __res -= __term;
	  if (__done(__term, __res))
	    break;
	  __w2k *= __w2;
	  __fact *= _Tp(__rzarg) / _Tp(__s + __rzarg)
		  * _Tp(__rzarg + 1) / _Tp(__s + __rzarg + 1);
	}
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
      const auto __pref = _Tp{2} * std::pow(_S_2pi, -(_Tp{1} - __s));
      // Here we factor up the ratio of Gamma(1 - s + k)/k! .
      // This ratio should be well behaved even for large k in the series
      // afterwards
      // Note that we have a problem for large s.
      // Since s is negative we evaluate the Gamma Function
      // on the positive real axis where it is real.
      auto __gam = std::exp(__ls);

      const auto __phase = __polar_pi(_Tp{1}, __s / _Tp{2});
      const auto __cp = std::real(__phase);
      const auto __sp = std::imag(__phase);
      // Here we add the expression that would result from ignoring
      // the zeta function in the series.
      const auto __p = _S_2pi - _S_i * __w;
      const auto __q = _S_2pi + _S_i * __w;
      // This can be optimized for real values of w
      __res += _S_i * __gam * (std::conj(__phase) * std::pow(__p, __s - _Tp{1})
	     - __phase * std::pow(__q, __s - _Tp{1}));
      // The above expression is the result of
      // sum_k Gamma(1+k-s)/k! * sin(pi (s-k)/2) (w/2/pi)^k
      // Therefore we only need to sample values of zeta(n) on the real axis
      // that really differ from one
      std::complex<_Tp> __sum = __sp * __gam * __riemann_zeta_m_1(_Tp{1} - __s);
      unsigned int __j = 1;
      __gam *= (_Tp{1} - __s);
      constexpr unsigned int __maxit = 200;
      _Terminator<std::complex<_Tp>> __done(__maxit);
      while (true)
	{
	  const auto __rzarg = _Tp(1 + __j) - __s;
	  const auto __rz = __riemann_zeta_m_1(__rzarg);
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
	  const auto __term =  __w2k * (__gam * __sine * __rz);
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
   * Compute the polylogarithm for negative integer order.
   * @f[
   *   Li_{-p}(e^w) = p!(-w)^{-(p+1)}
   *     - \sum_{k=0}^{\infty} \frac{B_{p+2k+q+1}}{(p+2k+q+1)!}
   *                           \frac{(p+2k+q)!}{(2k+q)!}w^{2k+q}
   * @f]
   * where @f$ q = (p+1) mod 2 @f$.
   *
   * @param __n the negative integer index @f$ n = -p @f$.
   * @param __w the argument w.
   * @return the value of the polylogarithm.
   */
  template<typename _Tp>
    std::complex<_Tp>
    __polylog_exp_neg(int __n, std::complex<_Tp> __w)
    {
      const auto _S_inf = std::numeric_limits<_Tp>::infinity();
      if (__gnu_cxx::__fp_is_zero(__w))
	return std::complex<_Tp>{0};
      else if (__gnu_cxx::__fp_is_equal(__w, _Tp{1}))
	return std::complex<_Tp>{_S_inf, _Tp{0}};
      else
	{
	  const int __p = -__n;
	  const int __pp = 1 + __p;
	  const int __q = __p & 1 ? 0 : 1;
	  const auto __w2 = __w * __w;
	  auto __wp = __p & 1 ? std::complex<_Tp>{1} : __w;
	  unsigned int __2k = __q;
	  auto __gam = std::__detail::__factorial<_Tp>(__p + __2k);
	  const auto __pfact = std::__detail::__factorial<_Tp>(__p);
	  auto __res = __pfact * std::pow(-__w, _Tp(-__pp));
	  auto __sum = std::complex<_Tp>{};
	  constexpr unsigned int __maxit = 300;
	  _Terminator<std::complex<_Tp>> __done(__maxit);
	  while (true)
	    {
	      const auto __id = (__p + __2k + 1) / 2;
	      if (__id == std::__detail::_Num_Euler_Maclaurin_zeta)
		break;
	      const auto __term = __gam * __wp
			* _Tp(std::__detail::_S_Euler_Maclaurin_zeta[__id]);
	      __sum += __term;
	      if (__done(__term, __sum))
		break;
	      __gam *= _Tp(__p + __2k + 1) / _Tp(__2k + 1)
		     * _Tp(__p + __2k + 2) / _Tp(__2k + 2);
	      __wp *= __w2;
	      __2k += 2;
	    }
	  __res -= __sum;
	  return __res;
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
      const auto __phase = __polar_pi(_Tp{1}, __s / _Tp{2});
      const auto __cp = std::real(__phase);
      const auto __sp = std::imag(__phase);
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
      // Factor this out for now so we can compare with sum.
      __res /= __pref;
      // Now comes the remainder of the series
      unsigned int __j = 0;
      constexpr unsigned int __maxit = 100;
      _Terminator<std::complex<_Tp>> __done(__maxit);
      auto __wup = __w / _S_2pi;
      auto __wbark = std::pow(__wup, _Tp(__m + 1));
      // It is 1 < 2 - s + m < 2 => Gamma(2-s+m) will not overflow
      // Here we factor up the ratio of Gamma(1 - s + k) / k!.
      // This ratio should be well behaved even for large k
      auto __gam = __gamma(_Tp(2 + __m) - __s) * __fact;
      std::complex<_Tp> __sum{};
      while (true)
	{
	  const auto __idx = __m + 1 + __j;
	  const auto __zetaarg = _Tp(1 + __idx) - __s;
	  const auto __rz = __riemann_zeta(__zetaarg);
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
	  const auto __term = __wbark * __sine * __gam * __rz;
	  __wbark *= __wup;
	  __gam *= __zetaarg / _Tp(1 + __idx);
	  ++__j;
	  __sum += __term;
	  if (__done(__term, __res + __sum))
	    break;
	}
      __res += __sum;
      return __pref * __res;
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
   * For real w the imaginary part of the polylog is given by
   * @f$ Im(Li_s(e^w)) = -\pi w^{s-1}/\Gamma(s) @f$.
   * Check this relation for any benchmark that you use.
   *
   * @param __s the real index s.
   * @param __w the large complex argument w.
   * @return the value of the polylogarithm.
   */
  template<typename _Tp>
    std::complex<_Tp>
    __polylog_exp_asymp(_Tp __s, std::complex<_Tp> __w)
    {
      const auto _S_pi = __gnu_cxx::__const_pi(__s);
      // wgamma = w^{s-1} / Gamma(s)
      auto __wgamma = std::pow(__w, __s - _Tp{1}) * __gamma_reciprocal(__s);
      auto __res = std::complex<_Tp>(_Tp{0}, -_S_pi) * __wgamma;
      // wgamma = w^s / Gamma(s+1)
      __wgamma *= __w / __s;
      constexpr unsigned int __maxiter = 100;
      _AsympTerminator<std::complex<_Tp>> __done(__maxiter);
      // zeta(0) w^s / Gamma(s + 1)
      std::complex<_Tp> __oldterm = -_Tp{0.5L} * __wgamma;
      __res += _Tp{2} * __oldterm;
      std::complex<_Tp> __term;
      auto __wq = _Tp{1} / (__w * __w);
      int __k = 1;
      while (true)
	{
	  __wgamma *= __wq * (__s + _Tp(1 - 2 * __k)) * (__s + _Tp(2 - 2 * __k));
	  __term = __riemann_zeta<_Tp>(2 * __k) * __wgamma;
	  __res += __done << _Tp{2} * __term;
	  if (__done(_Tp{2} * __term, __res))
	    break;
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
   *   Li_s(e^z) = \sum_{k=1}^{\infty} e^{kz} k^{-s}
   * @f]
   *
   * @param __s is an arbitrary type, integral or float.
   * @param __w something with a negative real part.
   * @return the value of the polylogarithm.
   */
  template<typename _PowTp, typename _Tp>
    _Tp
    __polylog_exp_sum(_PowTp __s, _Tp __w)
    {
      auto __ew = std::exp(__w);
      const auto __up = __ew;
      auto __res = __ew;
      unsigned int __maxiter = 500;
      _Terminator<_Tp> __done(__maxiter);
      bool __terminate = false;
      unsigned int __k = 2;
      while (!__terminate)
	{
	  __ew *= __up;
	  _Tp __temp = std::pow(__k, __s); // This saves us a type conversion
	  const auto __term = __ew / __temp;
	  __res += __term;
	  __terminate = __done(__term, __res);
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
    __polylog_exp_pos_int(unsigned int __s, std::complex<_Tp> __w)
    {
      const auto _S_2pi = __gnu_cxx::__const_2_pi(std::real(__w));
      const auto _S_pi = __gnu_cxx::__const_pi(std::real(__w));
      const auto _S_pi_2 = __gnu_cxx::__const_pi_half(std::real(__w));
      const auto _S_max_asymp = _Tp{5};
      const auto __rw = __w.real();
      const auto __iw = __w.imag();
      if (__gnu_cxx::__fp_is_real(__w)
	  && __gnu_cxx::__fp_is_equal(std::remainder(__iw, _S_2pi), _Tp{0}))
	{
	  if (__s == 1)
	    return std::numeric_limits<_Tp>::infinity();
	  else
	    return __riemann_zeta<_Tp>(__s);
	}
      else if (0 == __s)
	{
	  const auto __t = std::exp(__w);
	  return __gnu_cxx::__fp_is_zero(_Tp{1} - __t)
	       ? std::numeric_limits<_Tp>::quiet_NaN()
	       : __t / (_Tp{1} - __t);
	}
      else if (1 == __s)
	{
	  const auto __t = std::exp(__w);
	  return __gnu_cxx::__fp_is_zero(_Tp{1} - __t)
	       ? std::numeric_limits<_Tp>::quiet_NaN()
	       : -std::log(_Tp{1} - __t);
	}
      else
	{
	  if (__rw < -(_S_pi_2 + _S_pi / _Tp{5}))
	    // Choose the exponentially converging series
	    return __polylog_exp_sum(__s, __w);
	  else if (__rw < _S_max_asymp)
	    // The transition point chosen here, is quite arbitrary
	    // and needs more testing.
	    // The reductions of the imaginary part yield the same results
	    // as Mathematica.
	    // Necessary to improve the speed of convergence
	    return __polylog_exp_pos(__s, __clamp_pi(__w));
	  else
	    // Wikipedia says that this is required for Wood's formula.
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
    __polylog_exp_pos_int(unsigned int __s, _Tp __w)
    {
      const auto _S_pi = __gnu_cxx::__const_pi(std::real(__w));
      const auto _S_pi_2 = __gnu_cxx::__const_pi_half(std::real(__w));
      const auto _S_max_asymp = _Tp{5};
      if (__gnu_cxx::__fp_is_zero(__w))
	{
	  if (__s == 1)
	    return std::numeric_limits<_Tp>::infinity();
	  else
	    return __riemann_zeta<_Tp>(__s);
	}
      else if (__s == 0)
	{
	  const auto __t = std::exp(__w);
	  return __gnu_cxx::__fp_is_zero(_Tp{1} - __t)
	       ? std::numeric_limits<_Tp>::infinity()
	       : __t / (_Tp{1} - __t);
	}
      else if (__s == 1)
	{
	  const auto __t = std::exp(__w);
	  return __gnu_cxx::__fp_is_zero(_Tp{1} - __t)
	       ? -std::numeric_limits<_Tp>::infinity()
	       : -std::log(_Tp{1} - __t);
	}
      else
	{
	  if (__w < -(_S_pi_2 + _S_pi / _Tp{5}))
	    // Choose the exponentially converging series
	    return __polylog_exp_sum(__s, __w);
	  else if (__w < _S_max_asymp)
	    return __polylog_exp_pos(__s, __w);
	  else
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
    __polylog_exp_neg_int(int __s, std::complex<_Tp> __w)
    {
      const auto _S_2pi = __gnu_cxx::__const_2_pi(std::real(__w));
      const auto _S_pi = __gnu_cxx::__const_pi(std::real(__w));
      const auto _S_pi_2 = __gnu_cxx::__const_pi_half(std::real(__w));
      const auto _S_max_asymp = _Tp{5};
      if ((((-__s) & 1) == 0) && __gnu_cxx::__fp_is_imag(__w))
	{
	  // Now s is odd and w on the unit-circle.
	  const auto __iw = imag(__w); // Get imaginary part.
	  const auto __rem = std::remainder(__iw, _S_2pi);
	  if (__gnu_cxx::__fp_is_equal(std::abs(__rem), _Tp{0.5L}))
	    // Due to: Li_{-n}(-1) + (-1)^n Li_{-n}(1/-1) = 0.
	    return _Tp{0};
	  else
	    // No asymptotic expansion available... check the reduction.
	    return __polylog_exp_neg(__s, std::complex<_Tp>(__w.real(), __rem));
	}
      else
	{
	  if (std::real(__w) < -(_S_pi_2 + _S_pi / _Tp{5}))
	    // Choose the exponentially converging series
	    return __polylog_exp_sum(__s, __w);
	  else if (std::real(__w) < _S_max_asymp)
	    // Arbitrary transition point...
	    // The reductions of the imaginary part yield the same results
	    // as Mathematica.
	    // Necessary to improve the speed of convergence.
	    return __polylog_exp_neg(__s, __clamp_pi(__w));
	  else
	    // Wikipedia says that this clamping is required for Wood's formula.
	    return __polylog_exp_asymp(_Tp(__s), __clamp_0_m2pi(__w));
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
    __polylog_exp_neg_int(int __s, _Tp __w)
    {
      const auto _S_pi = __gnu_cxx::__const_pi(__w);
      const auto _S_pi_2 = __gnu_cxx::__const_pi_half(__w);
      const auto _S_max_asymp = _Tp{5};
      if (__w < -(_S_pi_2 + _S_pi / _Tp{5}))
	// Choose exponentially converging series.
	return __polylog_exp_sum(__s, __w);
      else if (__gnu_cxx::__fp_is_zero(__w))
	return std::numeric_limits<_Tp>::infinity();
      else if (__w < _S_max_asymp)
	// Arbitrary transition point less than 2 pi.
	return __polylog_exp_neg(__s, std::complex<_Tp>(__w));
      else
	return __polylog_exp_asymp(_Tp(__s), std::complex<_Tp>(__w));
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
    __polylog_exp_pos_real(_Tp __s, std::complex<_Tp> __w)
    {
      const auto _S_2pi = __gnu_cxx::__const_2_pi(__s);
      const auto _S_pi = __gnu_cxx::__const_pi(__s);
      const auto _S_pi_2 = __gnu_cxx::__const_pi_half(__s);
      const auto _S_max_asymp = _Tp{5};
      const auto __rw = __w.real();
      const auto __iw = __w.imag();
      if (__gnu_cxx::__fp_is_real(__w)
	  && __gnu_cxx::__fp_is_zero(std::remainder(__iw, _S_2pi)))
	{
	  if (__gnu_cxx::__fp_is_equal(__s, _Tp{1}))
	    return std::numeric_limits<_Tp>::infinity();
	  else
	    return __riemann_zeta(__s);
	}
      if (__rw < -(_S_pi_2 + _S_pi / _Tp{5}))
        // Choose exponentially converging series.
	return __polylog_exp_sum(__s, __w);
      if (__rw < _S_max_asymp)
	// Arbitrary transition point.
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
    __polylog_exp_pos_real(_Tp __s, _Tp __w)
    {
      const auto _S_pi = __gnu_cxx::__const_pi(__s);
      const auto _S_pi_2 = __gnu_cxx::__const_pi_half(__s);
      const auto _S_max_asymp = _Tp{5};
      if (__gnu_cxx::__fp_is_zero(__w))
	{
	  if (__gnu_cxx::__fp_is_equal(__s, _Tp{1}))
	    return std::numeric_limits<_Tp>::infinity();
	  else
	    return __riemann_zeta(__s);
	}
      if (__w < -(_S_pi_2 + _S_pi / _Tp{5}))
	// Choose exponentially converging series.
	return __polylog_exp_sum(__s, __w);
      if (__w < _S_max_asymp)
	// Arbitrary transition point
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
    __polylog_exp_neg_real(_Tp __s, std::complex<_Tp> __w)
    {
      const auto _S_pi = __gnu_cxx::__const_pi(__s);
      const auto _S_pi_2 = __gnu_cxx::__const_pi_half(__s);
      const auto _S_max_asymp = _Tp{5};
      const auto __rw = __w.real();
      const auto __iw = __w.imag();
      if (__rw < -(_S_pi_2 + _S_pi / _Tp{5}))
	// Choose exponentially converging series.
	return __polylog_exp_sum(__s, __w);
      else if (__rw < _S_max_asymp)
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
    __polylog_exp_neg_real(_Tp __s, _Tp __w)
    {
      const auto _S_pi = __gnu_cxx::__const_pi(__s);
      const auto _S_pi_2 = __gnu_cxx::__const_pi_half(__s);
      const auto _S_max_asymp = _Tp{5};
      if (__w < -(_S_pi_2 + _S_pi / _Tp{5}))
	// Choose exponentially converging series.
	return __polylog_exp_sum(__s, __w);
      else if (__w < _S_max_asymp)
	// Arbitrary transition point
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
      else if (__s > _Tp{25})
	// Cutoff chosen by some testing on the real axis.
	return __polylog_exp_sum(__s, __w);
      else
	{
	  const auto __p = __gnu_cxx::__fp_is_integer(__s, _Tp{5});
	  if (__p)
	    { // The order s is an integer.
	      if (__p() >= 0)
		return __polylog_exp_pos_int(__p(), __w);
	      else
		return __polylog_exp_neg_int(__p(), __w);
	    }
	  else
	    {
	      if (__s > _Tp{0})
		return __polylog_exp_pos_real(__s, __w);
	      else
		return __polylog_exp_neg_real(__s, __w);
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
	  const auto __n = __gnu_cxx::__fp_is_integer(__s, _Tp{5});
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
	      const auto __y = std::log(-__x);
	      return std::real(__polylog_exp(__s, _Tp{2} * __y)
				    * std::pow(_Tp{2}, _Tp{1} - __s)
			     - __polylog_exp(__s, __y));
	    }
	  else
	    {
	      const auto __y = std::log(__x);
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
	  const auto __t = _Tp{1} - __s;
	  const auto __lpe = __polylog_exp(__t, _S_i2pi * __a);
	  /// @todo This __hurwitz_zeta_polylog prefactor is prone to overflow.
	  /// positive integer orders s?
	  const auto __thing = std::exp(_Cmplx(_Tp{0}, -_S_pi_2 * __t));
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
   * The Dirichlet eta function, in terms of the polylogarithm, is
   * @f[
   *   \renewcommand\Re{\operatorname{Re}}
   *   \renewcommand\Im{\operatorname{Im}}
   *   \eta(s) = -\Re{Li_s(-1)}
   * @f]
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
	return -__polylog(std::real(__s), _Tp{-1});
      else
	std::__throw_domain_error(__N("__dirichlet_eta: Bad argument"));
    }

  /**
   * Return the Dirichlet eta function for real argument.
   * The Dirichlet eta function, in terms of the polylogarithm, is
   * @f[
   *   \renewcommand\Re{\operatorname{Re}}
   *   \renewcommand\Im{\operatorname{Im}}
   *   \eta(s) = -\Re{Li_s(-1)}
   * @f]
   *
   * @param __s  The real argument.
   * @return  The Dirichlet eta function.
   */
  template<typename _Tp>
    _Tp
    __dirichlet_eta(_Tp __s)
    {
      if (__isnan(__s))
	return __gnu_cxx::__quiet_NaN(__s);
      else if (__s < _Tp{0})
	{
	  const auto __p = __gnu_cxx::__fp_is_integer(__s, _Tp{5});
	  if (__p && (__p() % 2 == 0))
	    return _Tp{0};
	  else
	    {
	      const auto _S_pi = __gnu_cxx::__const_pi(__s);
	      const auto __sc = _Tp{1} - __s;
	      const auto __p2 = std::pow(_Tp{2}, -__sc);
	      return _Tp{2} * (_Tp{1} - __p2) / (_Tp{1} - _Tp{2} * __p2)
		   * std::pow(_S_pi, -__sc) * __s * __sin_pi(__s / _Tp{2})
		   * __gamma(-__s) * __dirichlet_eta(__sc);
	    }
	}
      else
	return -std::real(__polylog(__s, _Tp{-1}));
    }

  /**
   * Return the Dirichlet beta function.
   * Currently, s must be real (complex type but negligible imaginary part.)
   * Otherwise std::domain_error is thrown.
   * The Dirichlet beta function, in terms of the polylogarithm, is
   * @f[
   *   \renewcommand\Re{\operatorname{Re}}
   *   \renewcommand\Im{\operatorname{Im}}
   *   \beta(s) = \Im{Li_s(i)}
   * @f]
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
   * The Dirichlet beta function, in terms of the polylogarithm, is
   * @f[
   *   \renewcommand\Re{\operatorname{Re}}
   *   \renewcommand\Im{\operatorname{Im}}
   *   \beta(s) = \Im{Li_s(i)}
   * @f]
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
   * Return the Dirichlet lambda function for real argument.
   * @f[
   *   \lambda(s) = \frac{1}{2}(\zeta(s) + \eta(s))
   * @f]
   *
   * @param __s  The real argument.
   * @return  The Dirichlet lambda function.
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
   * Return Clausen's function of integer order m and complex argument @c z.
   * The notation and connection to polylog is from Wikipedia
   *
   * @param __m  The non-negative integral order.
   * @param __z  The complex argument.
   * @return  The complex Clausen function.
   */
  template<typename _Tp>
    std::complex<_Tp>
    __clausen(unsigned int __m, std::complex<_Tp> __z)
    {
      if (__isnan(__z))
	return std::numeric_limits<_Tp>::quiet_NaN();
      else if (__m == 0)
	std::__throw_domain_error(__N("__clausen: Non-positive order"));
      else
	{
	  const auto _S_i = std::complex<_Tp>{0, 1};
	  const auto __ple = __polylog_exp(_Tp(__m), _S_i * __z);
	  if (__m & 1)
	    return __ple;
	  else
	    return _S_i * std::conj(__ple);
	}
    }

  /**
   * Return Clausen's function of integer order m and real argument x.
   * The notation and connection to polylog is from Wikipedia
   *
   * @param __m  The integer order m >= 1.
   * @param __x  The real argument.
   * @return  The Clausen function.
   */
  template<typename _Tp>
    _Tp
    __clausen(unsigned int __m, _Tp __x)
    {
      if (__isnan(__x))
	return std::numeric_limits<_Tp>::quiet_NaN();
      else if (__m == 0)
	std::__throw_domain_error(__N("__clausen: Non-positive order"));
      else
	{
	  const auto _S_i = std::complex<_Tp>{0, 1};
	  const auto __ple = __polylog_exp(_Tp(__m), _S_i * __x);
	  if (__m & 1)
	    return std::real(__ple);
	  else
	    return std::imag(__ple);
	}
    }

  /**
   * Return Clausen's sine sum Sl_m for positive integer order m
   * and complex argument z.
   * @see https://en.wikipedia.org/wiki/Clausen_function
   *
   * @param __m  The integer order m >= 1.
   * @param __z  The complex argument.
   * @return  The Clausen sine sum Sl_m(w),
   */
  template<typename _Tp>
    _Tp
    __clausen_sl(unsigned int __m, std::complex<_Tp> __z)
    {
      if (__isnan(__z))
	return std::numeric_limits<_Tp>::quiet_NaN();
      else if (__m == 0)
	std::__throw_domain_error(__N("__clausen_sl: Non-positive order"));
      else
	{
	  const auto _S_i = std::complex<_Tp>{0, 1};
	  const auto __ple = __polylog_exp(_Tp(__m), _S_i * __z);
	  if (__m & 1)
	    return std::imag(__ple);
	  else
	    return std::real(__ple);
	}
    }

  /**
   * Return Clausen's sine sum Sl_m for positive integer order m
   * and real argument x.
   * @see https://en.wikipedia.org/wiki/Clausen_function
   *
   * @param __m  The integer order m >= 1.
   * @param __x  The real argument.
   * @return  The Clausen sine sum Sl_m(w),
   */
  template<typename _Tp>
    _Tp
    __clausen_sl(unsigned int __m, _Tp __x)
    {
      if (__isnan(__x))
	return std::numeric_limits<_Tp>::quiet_NaN();
      else if (__m == 0)
	std::__throw_domain_error(__N("__clausen_sl: Non-positive order"));
      else
	{
	  const auto _S_i = std::complex<_Tp>{0, 1};
	  const auto __ple = __polylog_exp(_Tp(__m), _S_i * __x);
	  if (__m & 1)
	    return std::imag(__ple);
	  else
	    return std::real(__ple);
	}
    }

  /**
   * Return Clausen's cosine sum Cl_m for positive integer order m
   * and complex argument w.
   * @see https://en.wikipedia.org/wiki/Clausen_function
   *
   * @param __m  The integer order m >= 1.
   * @param __z  The complex argument.
   * @return  The Clausen cosine sum Cl_m(w),
   */
  template<typename _Tp>
    _Tp
    __clausen_cl(unsigned int __m, std::complex<_Tp> __z)
    {
      if (__isnan(__z))
	return std::numeric_limits<_Tp>::quiet_NaN();
      else if (__m == 0)
	std::__throw_domain_error(__N("__clausen_cl: Non-positive order"));
      else
	{
	  const auto _S_i = std::complex<_Tp>{0, 1};
	  const auto __ple = __polylog_exp(_Tp(__m), _S_i * __z);
	  if (__m & 1)
	    return std::real(__ple);
	  else
	    return std::imag(__ple);
	}
    }

  /**
   * Return Clausen's cosine sum Cl_m for positive integer order m
   * and real argument w.
   * @see https://en.wikipedia.org/wiki/Clausen_function
   *
   * @param __m  The integer order m >= 1.
   * @param __x  The real argument.
   * @return  The real Clausen cosine sum Cl_m(w),
   */
  template<typename _Tp>
    _Tp
    __clausen_cl(unsigned int __m, _Tp __x)
    {
      if (__isnan(__x))
	return std::numeric_limits<_Tp>::quiet_NaN();
      else if (__m == 0)
	std::__throw_domain_error(__N("__clausen_cl: Non-positive order"));
      else
	{
	  const auto _S_i = std::complex<_Tp>{0, 1};
	  const auto __ple = __polylog_exp(_Tp(__m), _S_i * __x);
	  if (__m & 1)
	    return std::real(__ple);
	  else
	    return std::imag(__ple);
	}
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
      if (__isnan(__s) || __isnan(__x))
	return std::numeric_limits<_Tp>::quiet_NaN();
      else if (__s <= _Sp{-1})
	std::__throw_domain_error(__N("__fermi_dirac: "
				      "Order must be greater than -1"));
      else
	{
	  const auto _S_i = std::complex<_Tp>{0, 1};
	  const auto _S_pi = __gnu_cxx::__const_pi(__s);
	  return -std::real(__polylog_exp(__s + _Sp{1}, __x + _S_i * _S_pi));
	}
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
				      "Order must be greater than 0"));
      else
	return std::real(__polylog_exp(__s + _Sp{1}, __x));
    }

_GLIBCXX_END_NAMESPACE_VERSION
} // namespace __detail
} // namespace std

#endif // _GLIBCXX_BITS_SF_POLYLOG_TCC
