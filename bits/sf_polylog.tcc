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

/** @file bits/sf_polylog.tcc
 *  This is an internal header file, included by other library headers.
 *  Do not attempt to use it directly. @headername{cmath}
 */

//
// ISO C++ 14882 TR29124: Mathematical Special Functions
//

// Written by Florian Goth.
//
// References:
// (1) David C. Wood, "The Computation of Polylogarithms."
//

#ifndef _GLIBCXX_BITS_SF_POLYLOG_TCC
#define _GLIBCXX_BITS_SF_POLYLOG_TCC 1

#include <ext/math_const.h>

namespace std _GLIBCXX_VISIBILITY(default)
{
// Implementation-space details.
namespace __detail
{
_GLIBCXX_BEGIN_NAMESPACE_VERSION

  /**
   * A function to reliably compare two floating point numbers.
   *
   * @param a the left hand side.
   * @param b the right hand side
   * @return returns true if a and b are equal to zero or differ only by max(a,b)* 5* eps
   */
  template<typename _Tp>
    bool
    __fpequal(const _Tp& __a, const _Tp& __b)
    {
      constexpr auto _S_eps = std::numeric_limits<_Tp>::epsilon();
      constexpr auto _S_tol = _Tp{5} * _S_eps;
      bool __retval = true;
      if ((__a != _Tp{0}) || (__b != _Tp{0}))
	// Looks mean, but is necessary that the next line has sense.
	__retval = (std::abs(__a - __b) < std::max(std::abs(__a), std::abs(__b)) * _S_tol);
      return __retval;
    }


  /**
   * A function to reliably test a complex number for realness.
   *
   * @param w  The complex argument.
   * @return  true if Im(w) is zero within 5 * epsilon, false otherwize.
   */
  template<typename _Tp>
    inline bool
    __fpreal(const std::complex<_Tp>& __w)
    {
      constexpr auto _S_eps = std::numeric_limits<_Tp>::epsilon();
      constexpr auto _S_tol = _Tp{5} * _S_eps;
      return __fpequal(std::imag(__w), _Tp{0});
    }

  // Specialize for real numbers.
  template<typename _Tp>
    inline bool
    __fpreal(const _Tp)
    { return true; }


  /**
   * A function to reliably test a complex number for imaginaryness [?].
   *
   * @param w  The complex argument.
   * @return  true if Re(w) is zero within 5 * epsilon, false otherwize.
   */
  template<typename _Tp>
    inline bool
    __fpimag(const std::complex<_Tp>& __w)
    {
      constexpr auto _S_eps = std::numeric_limits<_Tp>::epsilon();
      constexpr auto _S_tol = _Tp{5} * _S_eps;
      return __fpequal(std::real(__w), _Tp{0});
    }

  // Specialize for real numbers.
  template<typename _Tp>
    inline bool
    __fpimag(const _Tp)
    { return false; }


  template<typename _Tp>
    std::complex<_Tp>
    __clamp_pi(std::complex<_Tp> __w)
    {
      constexpr auto _S_pi = __gnu_cxx::__math_constants<_Tp>::__pi;
      constexpr auto _S_i2pi = std::complex<_Tp>{0, _Tp{2} * _S_pi};
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
      constexpr auto _S_2pi = _Tp{2} * __gnu_cxx::__math_constants<_Tp>::__pi;
      constexpr auto _S_i2pi = std::complex<_Tp>{0, _S_2pi};
      while (__w.imag() > _Tp{0})
	__w = std::complex<_Tp>(__w.real(), __w.imag() - _S_2pi);
      while (__w.imag() <= -_S_2pi)
	__w = std::complex<_Tp>(__w.real(), __w.imag() + _S_2pi);
      return __w;
    }

  /**
   * A function to calculate the values of zeta at even positive integers.
   * For values smaller than thirty a table is used.
   *
   * @param k an integer at which we evaluate the Riemann zeta function.
   * @return zeta(k)
   */
  template<typename _Tp = double> 
  inline _Tp
  evenzeta(unsigned int __k)
  {
    // The following constants were calculated with Mathematica 8
    constexpr _Tp
    __data[]
    {
     -0.50000000000000000000000000,
      1.6449340668482264364724152,
      1.0823232337111381915160037,
      1.0173430619844491397145179,
      1.0040773561979443393786852,
      1.0009945751278180853371460,
      1.0002460865533080482986380,
      1.0000612481350587048292585,
      1.0000152822594086518717326,
      1.0000038172932649998398565,
      1.0000009539620338727961132,
      1.0000002384505027277329900,
      1.0000000596081890512594796,
      1.0000000149015548283650412,
      1.0000000037253340247884571,   
    };
    constexpr auto __maxk = 2 * sizeof(__data) / sizeof(_Tp);
    if (__k < __maxk)
      return __data[__k / 2];
    else
      return std::__detail::__riemann_zeta(static_cast<_Tp>(__k));
  }

  /**
   * This function treats the cases of positive integer index s.
   *
   * @f[
   *   Li_s(e^w) = \sum_{k=0, k != s-1} \zeta(s-k) w^k/k! + (H_{s-1} - log(-w)) w^(s-1)/(s-1)!
   * @f]
   * The radius of convergence is |w| < 2 pi.
   * Note that this series involves a log(-x).
   * gcc and Mathematica differ in their implementation of \log(e^(i \pi)):
   * gcc: \log(e^(+- i * \pi)) = +- i \pi
   * whereas Mathematica doesn't preserve the sign in this case: \log(e^(+- i\pi)) = +i \pi
   *
   * @param s the index s.
   * @param w the argument w.
   * @return the value of the polylogarithm.
   */
  template<typename _Tp>
    std::complex<_Tp>
    __polylog_exp_pos(unsigned int __s, std::complex<_Tp> __w)
    { // positive integer s
      // Optimization possibility: s are positive integers
      constexpr auto _S_2pi = _Tp{2} * __gnu_cxx::__math_constants<_Tp>::__pi;
      constexpr auto _S_pi = __gnu_cxx::__math_constants<_Tp>::__pi;
      constexpr auto _S_pipio6 = __gnu_cxx::__math_constants<_Tp>::__pi_sqr_div_6;
      std::complex<_Tp> __res = std::__detail::__riemann_zeta(_Tp(__s));
      auto __wk = __w;
      auto __fac = _Tp{1};
      auto __harmonicN = _Tp{1}; // HarmonicNumber_1
      for (unsigned int __k = 1; __k <= __s - 2; ++__k)
	{
	  __res += __wk * __fac * std::__detail::__riemann_zeta(_Tp(__s - __k));
	  __wk *= __w;
	  _Tp __temp = _Tp{1}/(_Tp{1} + __k);
	  __fac *= __temp;
	  __harmonicN += __temp;
	}
      // harmonicN now contains H_{s-1}.
      // fac should be 1/(n-1)!
      __res += (__harmonicN - std::log(-__w)) * __wk * __fac;
      __wk *= __w;
      __fac /= __s;
      __res -= __wk * __fac / _Tp{2};
      __wk *= __w;
      // Now comes the remainder of the series.
      const auto __pref = __wk / _S_pi / _S_2pi;
      const unsigned int __maxit = 200;
      unsigned int __j = 1;
      bool __terminate = false;
      __fac /= (__s + _Tp{1}); // (1/(n+1)!)
      __res -= _S_pipio6 * __fac * __pref; //subtract the zeroth order term.
      // Remainder of series.
      __fac *= _Tp{3} * _Tp{2} / (__s + _Tp{2}) / (__s + _Tp{3});
      auto __upfac = -(__w / _S_2pi) * (__w / _S_2pi);
      auto __w2 = __upfac;
      while (!__terminate) // Assume uniform convergence.
	{
          auto __rzarg = 2 * __j + 2;
          //auto __rz = std::__detail::__riemann_zeta(rzarg);
	  auto __rz = evenzeta<_Tp>(__rzarg);
          auto __term = (__rz * __fac) * __w2;
          __w2 *= __upfac;
          __fac *= __rzarg / _Tp(__rzarg + __s)
	      * (__rzarg + 1) / _Tp(__rzarg + __s + 1);
          ++__j;
          __terminate = (__fpequal(std::abs(__res - __pref * __term), std::abs(__res)) || (__j > __maxit));
          __res -= __pref * __term;
	}
      return __res;
    }

  /**
   * This function treats the cases of positive integer index s for real w.
   *
   * This specialization is worthwhile to catch the differing behaviour of log(x).
   * @f[
   *   Li_s(e^w) = \sum_{k=0, k != s-1} \zeta(s-k) w^k/k! + (H_{s-1} - log(-w)) w^(s-1)/(s-1)!
   * @f]
   * The radius of convergence is |w| < 2 pi.
   * Note that this series involves a log(-x).
   * The use of evenzeta yields a speedup of about 2.5.
   * gcc and Mathematica differ in their implementation of \log(e^(i \pi)):
   * gcc: \log(e^(+- i \pi)) = +- i \pi
   * whereas Mathematica doesn't preserve the sign in this case: \log(e^(+- i * \pi)) = +i \pi
   *
   * @param s the index s.
   * @param w the argument w
   * @return the value of the Polylogarithm
   */
  template<typename _Tp>
    inline std::complex<_Tp>
    __polylog_exp_pos(unsigned int __s, _Tp __w)
    { // positive integer s
      // Optimization possibility: s are positive integers.
      constexpr auto _S_2pi = _Tp{2} * __gnu_cxx::__math_constants<_Tp>::__pi;
      constexpr auto _S_pi = __gnu_cxx::__math_constants<_Tp>::__pi;
      auto __res = std::__detail::__riemann_zeta(_Tp(__s));
      auto __wk = __w;
      auto __fac = _Tp{1};
      auto __harmonicN = _Tp{1}; // HarmonicNumber_1
      for (unsigned int __k = 1; __k <= __s - 2; ++__k)
	{
	  __res += __wk * __fac * std::__detail::__riemann_zeta(_Tp(__s - __k));
	  __wk *= __w;
	  auto __temp = _Tp{1} / (_Tp{1} + __k);
	  __fac *= __temp;
	  __harmonicN += __temp;
	}
      // HarmonicN now contains H_{s-1}
      // fac should be 1/(n-1)!
      auto __imagtemp = __fac * __wk * (__harmonicN - std::log(std::complex<_Tp>(-__w, _Tp{0})));
      __res += real(__imagtemp);
      __wk *= __w;
      __fac /= __s;
      __res -= __wk * __fac / _Tp{2};
      __wk *= __w;
      // Now comes the remainder of the series.
      const auto __pref = __wk / _S_pi / _S_2pi;
      const unsigned int __maxit = 200;
      unsigned int __j = 1;
      bool __terminate = false;
      __fac /= (__s + _Tp{1}); // (1/(n+1)!)
      // Subtract the zeroth order term.
      __res -= _S_pi * _S_pi / _Tp{6} * __fac * __pref;
      // Remainder of series.
      __fac *= _Tp{3} * _Tp{2} / (__s + _Tp{2}) / (__s + _Tp{3});
      auto __upfac = -(__w / _S_2pi) * (__w / _S_2pi);
      auto __w2 = __upfac;
      while (!__terminate) // Assume convergence
	{
	  auto __rzarg = _Tp(2 * __j + 2);
	  auto __rz = evenzeta<_Tp>(__rzarg);
	  auto __term = __rz * __fac * __w2;
	  __w2 *= __upfac;
	  __fac *= __rzarg / (__rzarg + __s) * (__rzarg + _Tp{1}) / (__rzarg + __s + _Tp{1});
	  ++__j;
	  __terminate = (__fpequal(std::abs(__res - __pref * __term), std::abs(__res)) || (__j > __maxit));
	  __res -= __pref * __term;
	}
      return std::complex<_Tp>(__res, std::imag(__imagtemp));
    }

  /**
   * This function treats the cases of negative real index s.
   * Theoretical convergence is present for |w| < 2*pi.
   * We use an optimized version of
   * @f[
   *   Li_s(e^w) = \Gamma(1-s)(-w)^(s-1) + (2\pi)^(-s)/pi A_p(w)
   * @f]
   * @f[
   *   A_p(w) = \sum_k \Gamma(1+k-s)/k!\Sin(\pi/2*(s-k))(w/2/\pi)^k\zeta(1+k-s)
   * @f]
   * @param s  The index s.
   * @param w  The Argument w.
   * @return  The value of the polylogarithm.
   */
  template<typename _Tp>
    inline std::complex<_Tp>
    __polylog_exp_neg(_Tp __s, std::complex<_Tp> __w)
    {
      constexpr auto _S_i = std::complex<_Tp>{0, 1};
      constexpr auto _S_2pi = _Tp{2} * __gnu_cxx::__math_constants<_Tp>::__pi;
      constexpr auto _S_pi_2 = __gnu_cxx::__math_constants<_Tp>::__pi_half;
      // Basic general loop, but s is a negative quantity here
      // FIXME Large s makes problems.
      // The series should be rearrangeable so that we only need the ratio Gamma(1-s)/(2 pi)^s
      auto __ls = std::lgamma(_Tp{1} - __s);
      auto __res = std::exp(__ls - (_Tp{1} - __s) * std::log(-__w));
      const auto __wup = __w / _S_2pi;
      auto __w2 = __wup;
      auto __pref = _Tp{2} * std::pow(_S_2pi, -(_Tp{1} - __s));
      // here we factor up the ratio of Gamma(1 - s + k)/k! .
      // This ratio should be well behaved even for large k in the series afterwards
      // Note that we have a problem for large s
      // Since s is negative we evaluate the Gamma Function on the positive real axis where it is real.
      auto __gam = std::exp(__ls);

      auto __phase = std::polar(_Tp{1}, _S_pi_2 * __s);
      auto __cp = std::real(__phase);
      auto __sp = std::imag(__phase);
      // Here we add the expression that would result from ignoring the zeta function in the series.
      std::complex<_Tp> __expis(__cp, __sp);
      auto __p = _S_2pi - _S_i * __w;
      auto __q = _S_2pi + _S_i * __w;
      // This can be optimized for real values of w
      __res += _S_i * __gam * (std::conj(__expis) * std::pow(__p, __s - _Tp{1})
             - __expis * std::pow(__q, __s - _Tp{1}));
      // The above expression is the result of sum_k Gamma(1+k-s) /k! * sin(pi /2* (s-k)) * (w/2/pi)^k
      // Therefore we only need to sample values of zeta(n) on the real axis that really differ from one
      __res += __pref * (__sp * __gam * (std::__detail::__riemann_zeta(_Tp{1} - __s) - _Tp{1}));
      constexpr unsigned int __maxit = 200;
      unsigned int __j = 1;
      bool __terminate = false;
      __gam *= (_Tp{1} - __s);
      while (!__terminate) // Assume convergence
	{
	  auto __rzarg = _Tp(1 + __j) - __s;
	  auto __rz = std::__detail::__riemann_zeta_m_1(__rzarg);
	  _Tp __sine;
	  // Save the repeated recalculation of the sines
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
	  auto __term =  __w2 * (__gam * __sine * __rz);
	  __w2 *= __wup;
	  ++__j;
	  __gam  *= __rzarg / (__j); // equal to 1/(j+1) since we have incremented j in the line above
	  __terminate = (__fpequal(std::abs(__res + __pref * __term), std::abs(__res)) || (__j > __maxit));
	  __res += __pref * __term;
	}
      return __res;
    }

  /**
   * This function treats the cases of negative integer index s
   * which are multiples of two.
   *
   * In that case the sine occuring in the expansion occasionally
   *  takes on the value zero.
   * We use that to provide an optimized series for p = 2n:
   *
   * In the template parameter sigma we transport whether p = 4k (sigma = 1)
   * or p = 4k + 2  (sigma = -1)
   * @f[
   *   Li_p(e^w) = Gamma(1-p) (-w)^{p-1} - A_p(w) - sigma * B_p(w)
   * @f]
   * with
   * @f[
   *   A_p(w) = 2 (2\pi)^(p-1) (-p)! / (2 \pi)^(-p/2)
   *           (1 + w^2/(4 pi^2))^{-1/2 + p/2}
   *                          cos((1 - p) ArcTan(2 pi/ w))
   * @f]
   * and 
   * @f[
   *   B_p(w) = - 2 (2 pi)^(p-1) \sum \limits_{k = 0}^\infty \Gamma(2 + 2k - p)
   *           / (2k+1)! (-1)^k (w/2\pi)^(2k+1) (Zeta(2 + 2k - p) - _Tp{1})
   * @f]
   * This is suitable for |w| < 2 pi
   * The original series is (This might be worthwhile if we use
   * the already present table of the Bernoullis)
   * @f[
   *   Li_p(e^w) = Gamma(1-p) (-w)^{p-1} - sigma (2 \pi)^p / pi
   *              \sum \limits_{k = 0}^\infty
   *               \Gamma(2 + 2k - p)/ (2k+1)! (-1)^k (w/2\pi)^(2k+1)
   *             \zeta(2 + 2k - p)
   * @f]
   *
   * @param n the index n = 4k.
   * @param w The Argument w
   * @return the value of the Polylogarithm.
   */
  template<typename _Tp, int __sigma>
    inline std::complex<_Tp>
    __polylog_exp_neg_even(unsigned int __n, std::complex<_Tp> __w)
    {
      constexpr auto _S_2pi = _Tp{2} * __gnu_cxx::__math_constants<_Tp>::__pi;
      constexpr auto _S_pi = __gnu_cxx::__math_constants<_Tp>::__pi;
      const auto __np = 1 + __n;
      auto __lnp = std::lgamma(_Tp(__np));
      auto __res = std::exp(__lnp - _Tp(__np) * std::log(-__w));
      auto __wup = __w / _S_2pi;
      auto __wq = __wup * __wup;
      auto __pref = _Tp{2} * std::pow(_S_2pi, -_Tp(1 + __n));
      // Subtract the expression A_p(w)
      __res -= std::exp(__lnp - _Tp{0.5L} * __np * std::log(_Tp{1} + __wq))
	     * __pref * std::cos(_Tp(__np) * std::atan(_Tp{1} / __wup));
      unsigned int __k = 0;
      bool __terminate = false;
      constexpr unsigned int __maxit = 300;
      auto __gam = std::tgamma(_Tp(2 + __n));
      if (__sigma != 1)
	__pref = -__pref;
      while (!__terminate)
	{
	  auto __term = __gam * __riemann_zeta_m_1<_Tp>(2 * __k + 2 + __n)
		      * __wup;
	  __gam *= - _Tp(2 * __k + 2 + __n + 1) / _Tp(2 * __k + 2 + 1)
		 * _Tp(2 * __k + 2 + __n) / _Tp(2 * __k + 1 + 1);
	  __wup *= __wq;
	  __terminate = (__k > __maxit)
		     || __fpequal(std::abs(__res - __pref * __term),
				  std::abs(__res));
	  __res -= __pref * __term;
	  ++__k;
	}
    return __res;
  }

  /**
   * This function treats the cases of negative integer index s which are odd.
   *
   * In that case the sine occuring in the expansion occasionally vanishes.
   * We use that to provide an optimized series for p = 1 + 2k:
   * In the template parameter sigma we transport whether p = 1 + 4k (sigma = 1)
   * or p = 3 + 4k  (sigma = -1)
   *
   * @f[
   *   Li_p(e^w) = Gamma(1-p) * (-w)^{p-1} + sigma * A_p(w) - sigma * B_p(w)
   * @f]
   * with
   * @f[
   *   A_p(w) = 2 (2\pi)^(p-1) * Gamma(1-p) (1 + w^2/(4 pi^2))^{-1/2 + p/2} cos((1 - p) ArcTan(2 pi/ w))
   * @f]
   * and 
   * @f[
   *   B_p(w) = 2 (2 pi)^(p-1) * \sum \limits_{k = 0}^\infty \Gamma(1 + 2k - p)/ (2k)! (-w^2/4/\pi^2)^k (Zeta(1 + 2k - p) - _Tp{1})
   * @f]
   * This is suitable for |w| < 2 pi .
   * The use of evenzeta gives a speedup of about 50
   * The original series is (This might be worthwhile if we use the already present table of the Bernoullis)
   * @f[
   *   Li_p(e^w) = Gamma(1-p) * (-w)^{p-1} - sigma *2*(2 pi)^(p-1) * \sum \limits_{k = 0}^\infty \Gamma(1 + 2k - p)/ (2k)! (-1)^k (w/2/\pi)^(2k) Zeta(1 + 2k - p)
   * @f]
   *
   * @param n the index n = 4k.
   * @param w The Argument w.
   * @return The value of the Polylogarithm.
   */
  template<typename _Tp, int __sigma>
    inline std::complex<_Tp>
    __polylog_exp_neg_odd(unsigned int __n, std::complex<_Tp> __w)
    {
      constexpr auto _S_pi = __gnu_cxx::__math_constants<_Tp>::__pi;
      const unsigned int __np = 1 + __n;
      auto __lnp = std::lgamma(_Tp(__np));
      auto __res = std::exp(__lnp - _Tp(__np) * std::log(-__w));
      constexpr auto __itp = _Tp{1} / (_Tp{2} * _S_pi);
      auto __wq = -__w * __itp * __w * __itp;
      auto __pref = _Tp{2} * std::pow(__itp, _Tp(__np));
      // Subtract the expression A_p(w)
      __res += std::exp(__lnp - _Tp{0.5L} * __np * std::log(_Tp{1} - __wq))
	     * __pref * std::cos(_Tp(__np)
			* std::atan(std::complex<_Tp>{2 * _S_pi} / __w));
      if (__sigma != 1)
	__pref = -__pref;
      bool __terminate = false;
      constexpr unsigned int __maxit = 300;
      _Tp __gam = std::exp(__lnp);
      // zeroth order
      __res -= __pref * __gam * (evenzeta<_Tp>(__np) - _Tp{1});
      unsigned int __k = 0;
      auto __wup = __wq;
      while (!__terminate)
	{
	  auto __zk = 2 * __k;
	  __gam *= _Tp(__zk + __np) / _Tp(1 + __zk)
		 * _Tp(1 + __zk + __np) / _Tp(__zk + 2);
	  auto __term = (__gam * __riemann_zeta_m_1<_Tp>(__zk + 2 + __np)) * __wup;
	  __wup *= __wq;
	  __terminate = __k > __maxit
		     || __fpequal(std::abs(__res - __pref * __term),
				  std::abs(__res));
	  __res -= __pref * __term;
	  ++__k;
	}
      return __res;
  }

  /**
   * This function treats the cases of negative integer index s and branches accordingly
   *
   * @param s the integer index s.
   * @param w The Argument w
   * @return The value of the Polylogarithm evaluated by a suitable function.
   */
  template<typename _Tp>
    inline std::complex<_Tp> __polylog_exp_neg(int __s, std::complex<_Tp> __w)
    { // negative integer __s
      const auto __n = -__s;
      switch (__n % 4)
      {
      case 0:
	return __polylog_exp_neg_even<_Tp, 1>(__n, __w);
      case 1:
	return __polylog_exp_neg_odd<_Tp, 1>(__n, __w);
      case 2:
	return __polylog_exp_neg_even<_Tp, -1>(__n, __w);
      case 3:
	return __polylog_exp_neg_odd<_Tp, -1>(__n, __w);
	break;
      }
    }

  /**
   * This function treats the cases of positive real index s.
   *
   * The defining series is
   * @f[
   *   Li_s(e^w) = A_s(w) + B_s(w) + \Gamma(1-s)(-w)^(s-1)
   * @f]
   * with
   * @f[
   *   A_s(w) = \sum_{k=0}^{m} \zeta(s-k)w^k/k!
   * @f]
   * @f[
   *   B_s(w) = \sum_{k=m+1}^\infty \sin(\pi/2(s-k))
   *             \Gamma(1-s+k)\zeta(1-s+k) (w/2/\pi)^k/k!
   * @f]
   *
   * @param s the positive real index s.
   * @param w The argument w.
   * @return the value of the polylogarithm.
   */
  template<typename _Tp>
    inline std::complex<_Tp>
    __polylog_exp_pos(_Tp __s, std::complex<_Tp> __w)
    { // positive s
      constexpr auto _S_2pi = _Tp{2} * __gnu_cxx::__math_constants<_Tp>::__pi;
      constexpr auto _S_pi = __gnu_cxx::__math_constants<_Tp>::__pi;
      constexpr auto _S_pi_2 = __gnu_cxx::__math_constants<_Tp>::__pi_half;
      std::complex<_Tp> __res = std::__detail::__riemann_zeta(__s);
      auto __wk = __w;
      auto __phase = std::polar(_Tp{1}, _S_pi_2 * __s);
      auto __cp = std::real(__phase);
      auto __sp = std::imag(__phase);
      // This is \Gamma(1-s)(-w)^(s-1)
      __res += _S_pi / (_Tp{2} * __sp * __cp)
	  * std::exp(-std::lgamma(__s) + (__s - _Tp{1}) * std::log(-__w));
      auto __fac = _Tp{1};
      const auto __m = static_cast<unsigned int>(std::floor(__s));
      for (unsigned int __k = 1; __k <= __m; ++__k)
	{
          __res += __wk * __fac
		 * std::__detail::__riemann_zeta(__s - _Tp(__k));
          __wk *= __w;
          __fac /= _Tp(1 + __k);
	}
      // fac should now be 1/(m+1)!
      const auto __pref = _Tp{2} * std::pow(_S_2pi, __s - _Tp{1});
      // Now comes the remainder of the series
      constexpr unsigned int __maxit = 100;
      unsigned int __j = 0;
      bool __terminate = false;
      auto __wup = __w / _S_2pi;
      auto __w2 = std::pow(__wup, _Tp(__m + 1));
      // It is 1 < 2 - s + m < 2 => Gamma(2-s+m) will not overflow
      // Here we factor up the ratio of Gamma(1 - s + k) / k!.
      // This ratio should be well behaved even for large k
      auto __gam = std::tgamma(_Tp(2 + __m) - __s) * __fac;
      while (!__terminate)
	{ // FIXME: optimize.
	  auto __idx = __m + 1 + __j;
	  auto __zetaarg = _Tp(1 + __idx) - __s;
	  auto __rz = std::__detail::__riemann_zeta(__zetaarg);
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
          auto __term = __w2 * __sine * __gam * __rz;
          __w2 *= __wup;
          __gam *= __zetaarg / _Tp(1 + __idx);
          ++__j;
          __terminate = (__fpequal(std::abs(__res + __pref * __term),
				   std::abs(__res))
		     || (__j > __maxit));
          __res += __pref * __term;
	}
      return __res;
    }

  /**
   * This function implements the asymptotic series for the polylog.
   * It is given by 2 \sum \limits_{k=0}^\infty \zeta(2k) w^{s-2k}/Gamma(s-2k+1) -i \pi w^(s-1)/Gamma(s)
   * for Re(w) >> 1
   *
   * Don't check this against Mathematica 8.
   * For real u the imaginary part of the polylog is given by Im(Li_s(e^u)) = - \pi u^{s-1}/Gamma(s)
   * Check this relation for any benchmark that you use.
   * The use of evenzeta leads to a speedup of about 1000.
   *
   * @param s the index s.
   * @param w the large argument w.
   * @return the value of the polylogarithm.
   */
  template<typename _Tp>
    std::complex<_Tp>
    __polylog_exp_asymp(const _Tp __s, std::complex<_Tp> __w)
    { // asymptotic expansion
      constexpr auto _S_pi = __gnu_cxx::__math_constants<_Tp>::__pi;
      // wgamma = w^(s-1) / Gamma(s)
      auto __wgamma = std::exp((__s - _Tp{1}) * std::log(__w) - std::lgamma(__s));
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
          __wgamma *= __wq * (__s + _Tp(1 - 2 * __k)) * (__s + _Tp(2 - 2 * __k));
	  __term = evenzeta<_Tp>(2 * __k) * __wgamma;
          if (std::abs(__term) > std::abs(__oldterm))
	    __terminate = true; // Termination due to failure of asymptotic expansion
          if (__fpequal(std::abs(__res + _Tp{2} * __term), std::abs(__res)))
	    __terminate = true; // Precision goal reached.
          if (__k > __maxiter)
	    __terminate = true; // Stop the iteration somewhen
          if (!__terminate)
            {
              __res += _Tp{2} * __term;
              __oldterm = __term;
              ++__k;
            }
	}
      return __res;
  }

  /**
   * Theoretical convergence for Re(w) < 0.
   *
   * Seems to beat the other expansions for Re(w) < -pi/2 - pi/5.
   * Note that this is an implementation of the basic series:
   * @f[
   *   Li_s(e^z) = \sum_{k=1} e^(k*z) * k^(-s)
   * @f]
   *
   * @param s is an arbitrary type, Integer or float.
   * @param w something with a negative real part.
   * @return the value of the polylogarithm.
   */
  template<typename _PowTp, typename _Tp>
    inline _Tp
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
	  __terminate = (__fpequal(std::abs(__res + __term), std::abs(__res))) || (__k > __maxiter);
	  __res += __term;
	  ++__k;
	}
      return __res;
  }

  /**
   * Here s is a positive integer and the function descends
   * into the different kernels depending on w.
   *
   * @param s a positive integer.
   * @param w an arbitrary complex number.
   * @return The value of the polylogarithm.
   */
  template<typename _Tp>
    std::complex<_Tp>
    __polylog_exp_int_pos(const unsigned int __s, std::complex<_Tp> __w)
    {
      constexpr auto _S_2pi = _Tp{2} * __gnu_cxx::__math_constants<_Tp>::__pi;
      constexpr auto _S_pi = __gnu_cxx::__math_constants<_Tp>::__pi;
      constexpr auto _S_pi_2 = __gnu_cxx::__math_constants<_Tp>::__pi_half;
      auto __rw = __w.real();
      auto __iw = __w.imag();
      if (__fpreal(__w) && __fpequal(std::remainder(__iw, _S_2pi), _Tp{0}))
	{
	  if (__s > 1)
	    return std::__detail::__riemann_zeta(_Tp(__s));
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
	    // The transition point chosen here, is quite arbitrary and needs more testing.
	    // The reductions of the imaginary part yield the same results as Mathematica.
	    // Necessary to improve the speed of convergence
	    return __polylog_exp_pos(__s, __clamp_pi(__w));
	  else
	    // Wikipedia says that this is required for Wood's formula.
	    // FIXME: The series should terminate after a finite number of terms.
	    return __polylog_exp_asymp(static_cast<_Tp>(__s), __clamp_0_m2pi(__w));
	}
    }

  /**
   * Here s is a positive integer and the function descends
   * into the different kernels depending on w.
   *
   * @param s a positive integer
   * @param w an arbitrary real argument w
   * @return the value of the polylogarithm.
   */
  template<typename _Tp>
    std::complex<_Tp>
    __polylog_exp_int_pos(const unsigned int __s, _Tp __w)
    {
      constexpr auto _S_pi = __gnu_cxx::__math_constants<_Tp>::__pi;
      constexpr auto _S_pi_2 = __gnu_cxx::__math_constants<_Tp>::__pi_half;
      if (__fpequal(__w, _Tp{0}))
	{
	  if (__s > 1)
	    return std::__detail::__riemann_zeta(_Tp(__s));
	  else
	    return std::numeric_limits<_Tp>::infinity();
	}
      else if (__s == 0)
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
	    return __polylog_exp_negative_real_part(__s, std::complex<_Tp>(__w));
          else if (__w < _Tp{6})
            // The transition point chosen here, is quite arbitrary and needs more testing.
	    return __polylog_exp_pos(__s, __w);
	  else
	    // FIXME: The series should terminate after a finite number of terms.
	    return __polylog_exp_asymp(static_cast<_Tp>(__s), std::complex<_Tp>(__w));
	}
    }

  /**
   * This treats the case where s is a negative integer.
   *
   * @param s a negative integer.
   * @param w an arbitrary complex number
   * @return the value of the polylogarith,.
   */
  template<typename _Tp>
    std::complex<_Tp>
    __polylog_exp_int_neg(const int __s, std::complex<_Tp> __w)
    {
      constexpr auto _S_2pi = _Tp{2} * __gnu_cxx::__math_constants<_Tp>::__pi;
      constexpr auto _S_pi = __gnu_cxx::__math_constants<_Tp>::__pi;
      constexpr auto _S_pi_2 = __gnu_cxx::__math_constants<_Tp>::__pi_half;
      if ((((-__s) & 1) == 0) && __fpequal(std::real(__w), _Tp{0}))
	{
	  // Now s is odd and w on the unit-circle.
	  auto __iw = imag(__w);  //get imaginary part
	  auto __rem = std::remainder(__iw, _S_2pi);
	  if (__fpequal(std::abs(__rem), _Tp{0.5L}))
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
	    // The reductions of the imaginary part yield the same results as Mathematica.
	    // Necessary to improve the speed of convergence.
	    return __polylog_exp_neg(__s, __clamp_pi(__w));
	  else
	    // Wikipedia says that this clamping is required for Wood's formula.
	    // FIXME: The series should terminate after a finite number of terms.
	    return __polylog_exp_asymp(static_cast<_Tp>(__s), __clamp_0_m2pi(__w));
	}
    }

  /**
   * This treats the case where s is a negative integer and w is a real.
   *
   * @param s a negative integer.
   * @param w the argument.
   * @return the value of the polylogarithm.
   */
  template<typename _Tp>
    inline std::complex<_Tp>
    __polylog_exp_int_neg(const int __s, _Tp __w)
    {
      constexpr auto _S_pi = __gnu_cxx::__math_constants<_Tp>::__pi;
      constexpr auto _S_pi_2 = __gnu_cxx::__math_constants<_Tp>::__pi_half;
      if (__w < -(_S_pi_2 + _S_pi / _Tp{5})  ) // Choose the exponentially converging series
	return __polylog_exp_negative_real_part(__s, std::complex<_Tp>(__w));
      else if (__fpequal(__w, _Tp{0}))
	return std::numeric_limits<_Tp>::infinity();
      else if (__w < _Tp{6}) // Arbitrary transition point...
	return __polylog_exp_neg(__s, std::complex<_Tp>(__w));
      else
	// FIXME: The series should terminate after a finite number of terms.
	return __polylog_exp_asymp(static_cast<_Tp>(__s), std::complex<_Tp>(__w));
    }

  /**
   * Return the polylog where s is a positive real value and for complex argument.
   *
   * @param s A positive real number.
   * @param w the complex argument.
   * @return The value of the polylogarithm.
   */
  template<typename _Tp>
    std::complex<_Tp>
    __polylog_exp_real_pos(_Tp __s, std::complex<_Tp> __w)
    {
      constexpr auto _S_2pi = _Tp{2} * __gnu_cxx::__math_constants<_Tp>::__pi;
      constexpr auto _S_pi = __gnu_cxx::__math_constants<_Tp>::__pi;
      constexpr auto _S_pi_2 = __gnu_cxx::__math_constants<_Tp>::__pi_half;
      auto __rw = __w.real();
      auto __iw = __w.imag();
      if (__fpreal(__w) && __fpequal(std::remainder(__iw, _S_2pi), _Tp{0}))
	{
	  if (__s > _Tp{1})
	    return std::__detail::__riemann_zeta(__s);
	  else
	    return std::numeric_limits<_Tp>::infinity();
	}
      if (__rw < -(_S_pi_2 + _S_pi/_Tp{5})) // Choose the exponentially converging series
	return __polylog_exp_negative_real_part(__s, __w);
      if (__rw < _Tp{6}) // arbitrary transition point
	// The reductions of the imaginary part yield the same results as Mathematica then.
	// Branch cuts??
	return __polylog_exp_pos(__s, __clamp_pi(__w));
      else
	// Wikipedia says that this is required for Wood's formula
	return __polylog_exp_asymp(__s, __clamp_0_m2pi(__w));
    }

  /**
   * Return the polylog where s is a positive real value and the argument is real.
   *
   * @param s  A positive real number tht does not reduce to an integer.
   * @param w  The real argument w.
   * @return  The value of the polylogarithm.
   */
  template<typename _Tp>
    std::complex<_Tp>
    __polylog_exp_real_pos(_Tp __s, _Tp __w)
    {
      constexpr auto _S_pi = __gnu_cxx::__math_constants<_Tp>::__pi;
      constexpr auto _S_pi_2 = __gnu_cxx::__math_constants<_Tp>::__pi_half;
      if (__fpequal(__w, _Tp{0}))
	{
          if (__s > _Tp{1})
	    return std::__detail::__riemann_zeta(__s);
	  else
	    return std::numeric_limits<_Tp>::infinity();
	}
      if (__w < -(_S_pi_2 + _S_pi / _Tp{5})) // Choose the exponentially converging series
	return __polylog_exp_negative_real_part(__s, __w);
      if (__w < _Tp{6}) // arbitrary transition point
	return __polylog_exp_pos(__s, std::complex<_Tp>(__w));
      else
	return __polylog_exp_asymp(__s, std::complex<_Tp>(__w));
  }

  /**
   * Return the polylog where s is a negative real value and for complex argument.
   * Now we branch depending on the properties of w in the specific functions
   *
   * @param s  A negative real value that does not reduce to a negative integer.
   * @param w  The complex argument.
   * @return  The value of the polylogarithm.
   */
  template<typename _Tp>
    std::complex<_Tp>
    __polylog_exp_real_neg(_Tp __s, std::complex<_Tp> __w)
    {
      constexpr auto _S_pi = __gnu_cxx::__math_constants<_Tp>::__pi;
      constexpr auto _S_pi_2 = __gnu_cxx::__math_constants<_Tp>::__pi_half;
      auto __rw = __w.real();
      auto __iw = __w.imag();
      if (__rw < -(_S_pi_2 + _S_pi/_Tp{5})) // Choose the exponentially converging series
	return __polylog_exp_negative_real_part(__s, __w);
      else if (__rw < 6) // arbitrary transition point
	// The reductions of the imaginary part yield the same results as Mathematica then.
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
   * @param s  A negative real value.
   * @param w  A real argument.
   * @return  The value of the polylogarithm.
   */
  template<typename _Tp>
    std::complex<_Tp>
    __polylog_exp_real_neg(_Tp __s, _Tp __w)
    {
      constexpr auto _S_pi = __gnu_cxx::__math_constants<_Tp>::__pi;
      constexpr auto _S_pi_2 = __gnu_cxx::__math_constants<_Tp>::__pi_half;
      if (__w < -(_S_pi_2 + _S_pi / _Tp{5})) // Choose the exponentially converging series
	return __polylog_exp_negative_real_part(__s, std::complex<_Tp>(__w));
      else if (__w < _Tp{6}) // arbitrary transition point
	return __polylog_exp_neg(__s, std::complex<_Tp>(__w));
      else
	return __polylog_exp_asymp(__s, std::complex<_Tp>(__w));
    }

  /**
   * This is the frontend function which calculates Li_s(e^w)
   * First we branch into different parts depending on the properties of s.
   * This function is the same irrespective of a real or complex w,
   * hence the template parameter ArgType.
   *
   * @note: I *really* wish we could return a variant<Tp, std::complex<Tp>>.
   *
   * @param s  The real order.
   * @param w  The real or complex argument.
   * @return  The real or complex value of Li_s(e^w).
   */
  template<typename _Tp, typename ArgType>
    __gnu_cxx::__promote_num_t<std::complex<_Tp>, ArgType>
    __polylog_exp(_Tp __s, ArgType __w)
    {
      if (__isnan(__s) || __isnan(__w))
	return std::numeric_limits<_Tp>::quiet_NaN();
      else if (__s > _Tp{25}) // Cutoff chosen by some testing on the real axis.
	return __polylog_exp_negative_real_part(__s, __w);
      else if (__fpequal<_Tp>(std::rint(__s), __s))
	{
          // In this branch of the if statement, s is an integer
          int __p = int(std::lrint(__s));
	  if (__p > 0)
	    return __polylog_exp_int_pos(__p, __w);
	  else
	    return __polylog_exp_int_neg(__p, __w);
	}
      else
	{
	  if (__s > _Tp{0})
	    return __polylog_exp_real_pos(__s, __w);
	  else
	    return __polylog_exp_real_neg(__s, __w);
	}
    }

  /**
   * Return the polylog Li_s(x) for two real arguments.
   *
   * @param s  The real index.
   * @param x  The real argument.
   * @return The complex value of the polylogarithm.
   */
  template<typename _Tp>
    _Tp
    __polylog(_Tp __s, _Tp __x)
    {
      if (__isnan(__s) || __isnan(__x))
	return std::numeric_limits<_Tp>::quiet_NaN();
      else if (__fpequal(__x, _Tp{0}))
	return _Tp{0}; // According to Mathematica
      else if (__x < _Tp{0})
	{ // Use the reflection formula to access negative values.
          auto __xp = -__x;
          auto __y = std::log(__xp);
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

  /**
   * Return the polylog in those cases where we can calculate it.
   *
   * @param s  The real index.
   * @param w  The complex argument.
   * @return  The complex value of the polylogarithm.
   */
  template<typename _Tp>
    inline std::complex<_Tp>
    __polylog(_Tp __s, std::complex<_Tp> __w)
    {
      if (__isnan(__s) || __isnan(__w))
	return std::numeric_limits<_Tp>::quiet_NaN();
      else if (__fpequal(std::imag(__w), _Tp{0}))
	return __polylog(__s, std::real(__w));
      else
	return __polylog_exp(__s, std::log(__w));
    }

  /**
   * Return the Dirichlet eta function.
   * Currently, w must be real (complex type but negligible imaginary part.)
   * Otherwise std::domain_error is thrown.
   *
   * @param w  The complex (but on-real-axis) argument.
   * @return  The complex Dirichlet eta function.
   * @throw  std::domain_error if the argument has a significant imaginary part.
   */
  template<typename _Tp>
    inline std::complex<_Tp>
    __dirichlet_eta(std::complex<_Tp> __w)
    {
      if (__isnan(__w))
	return std::numeric_limits<_Tp>::quiet_NaN();
      else if (__fpequal(std::imag(__w), _Tp{0}))
	return -__polylog(__w.real(), _Tp{-1});
      else
	std::__throw_domain_error(__N("__dirichlet_eta: Bad argument"));
    }

  /**
   *  Return the Dirichlet eta function for real argument.
   *
   *  @param w  The real argument.
   *  @return  The Dirichlet eta function.
   */
  template<typename _Tp>
    inline _Tp
    __dirichlet_eta(_Tp __w)
    {
      if (__isnan(__w))
	return std::numeric_limits<_Tp>::quiet_NaN();
      else
	return -std::real(__polylog(__w, _Tp{-1}));
    }

  /**
   * Return the Dirichlet beta function.
   * Currently, w must be real (complex type but negligible imaginary part.)
   * Otherwise std::domain_error is thrown.
   *
   * @param w  The complex (but on-real-axis) argument.
   * @return  The Dirichlet Beta function of real argument.
   * @throw  std::domain_error if the argument has a significant imaginary part.
   */
  template<typename _Tp>
    inline _Tp
    __dirichlet_beta(std::complex<_Tp> __w)
    {
      constexpr auto _S_i = std::complex<_Tp>{0, 1};
      if (__isnan(__w))
	return std::numeric_limits<_Tp>::quiet_NaN();
      else if (__fpequal(std::imag(__w), _Tp{0}))
	return std::imag(__polylog(__w.real(), _S_i));
      else
	std::__throw_domain_error(__N("__dirichlet_beta: Bad argument."));
    }

  /**
   * Return the Dirichlet beta function for real argument.
   *
   * @param w  The real argument.
   * @return  The Dirichlet Beta function of real argument.
   */
  template<typename _Tp>
    inline _Tp
    __dirichlet_beta(_Tp __w)
    {
      constexpr auto _S_i = std::complex<_Tp>{0, 1};
      if (__isnan(__w))
	return std::numeric_limits<_Tp>::quiet_NaN();
      else
	return std::imag(__polylog(__w, _S_i));
    }

  /**
   * Return Clausen's function of integer order m and complex argument w.
   * The notation and connection to polylog is from Wikipedia
   *
   * @param w  The complex argument.
   * @return  The complex Clausen function.
   */
  template<typename _Tp>
    inline std::complex<_Tp>
    __clausen(unsigned int __m, std::complex<_Tp> __w)
    {
      constexpr auto _S_i = std::complex<_Tp>{0, 1};
      auto __ple = __polylog_exp(_Tp(__m), _S_i * __w);
      if (__isnan(__w))
	return std::numeric_limits<_Tp>::quiet_NaN();
      else if (__m == 0)
	std::__throw_domain_error(__N("__clausen: Non-positive order"));
      else if (__m & 1)
	return __ple;
      else
	return std::conj(__ple);
    }

  /**
   * Return Clausen's sine sum Sl for positive integer order m
   * and complex argument w.
   * The notation and connection to polylog from Wikipedia
   *
   * @param w  The complex argument.
   * @return  The Clausen sine sum S_m(w),
   */
  template<typename _Tp>
    inline _Tp
    __clausen_s(unsigned int __m, std::complex<_Tp> __w)
    {
      constexpr auto _S_i = std::complex<_Tp>{0, 1};
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
   * Return Clausen's cosine sum Cl for positive integer order m
   * and complex argument w.
   * The notation and connection to polylog from Wikipedia
   *
   * @param w  The complex argument.
   * @return  The Clausen cosine sum C_m(w),
   */
  template<typename _Tp>
    inline _Tp
    __clausen_c(unsigned int __m, std::complex<_Tp> __w)
    {
      constexpr auto _S_i = std::complex<_Tp>{0, 1};
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

_GLIBCXX_END_NAMESPACE_VERSION
} // namespace __detail
} // namespace std

#endif // _GLIBCXX_BITS_SF_POLYLOG_TCC