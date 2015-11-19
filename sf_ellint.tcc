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

/** @file bits/sf_ellint.tcc
 *  This is an internal header file, included by other library headers.
 *  Do not attempt to use it directly. @headername{cmath}
 */

//
// ISO C++ 14882 TR29124: Mathematical Special Functions
//

// Written by Edward Smith-Rowland.
//
// References:
//   (1)  B. C. Carlson Numer. Math. 33, 1 (1979)
//   (2)  B. C. Carlson, Special Functions of Applied Mathematics (1977)
//   (3)  The Gnu Scientific Library, http://www.gnu.org/software/gsl
//   (4)  Numerical Recipes in C, 2nd ed, by W. H. Press, S. A. Teukolsky,
//        W. T. Vetterling, B. P. Flannery, Cambridge University Press
//        (1992), pp. 261-269

#ifndef _GLIBCXX_BITS_SF_ELLINT_TCC
#define _GLIBCXX_BITS_SF_ELLINT_TCC 1

namespace std _GLIBCXX_VISIBILITY(default)
{
// Implementation-space details.
namespace __detail
{
_GLIBCXX_BEGIN_NAMESPACE_VERSION

  /**
   *   @brief Return the Carlson elliptic function @f$ R_F(x,y,z) @f$
   *          of the first kind.
   *
   *   The Carlson elliptic function of the first kind is defined by:
   *   @f[
   *       R_F(x,y,z) = \frac{1}{2} \int_0^\infty
   *                 \frac{dt}{(t + x)^{1/2}(t + y)^{1/2}(t + z)^{1/2}}
   *   @f]
   *
   *   @param  __x  The first of three symmetric arguments.
   *   @param  __y  The second of three symmetric arguments.
   *   @param  __z  The third of three symmetric arguments.
   *   @return  The Carlson elliptic function of the first kind.
   */
  template<typename _Tp>
    _Tp
    __ellint_rf(_Tp __x, _Tp __y, _Tp __z)
    {
      using _Val = __num_traits_t<_Tp>;
      const _Val __r = std::numeric_limits<_Val>::epsilon();
      _Tp __xt = __x;
      _Tp __yt = __y;
      _Tp __zt = __z;
      _Tp __a0 = (__x + __y + __z) / _Val(3);
      _Val __q = std::pow( _Val(3) * __r, -_Val(1) / _Val(6) )
	       * std::max(std::abs(__a0 - __z),
			  std::max(std::abs(__a0 - __x),
				   std::abs(__a0 - __y)));
      _Tp __a = __a0;
      _Val __f = _Val(1);

      while (true)
	{
	  _Tp __lambda = std::sqrt(__xt) * std::sqrt(__yt)
		       + std::sqrt(__yt) * std::sqrt(__zt)
		       + std::sqrt(__zt) * std::sqrt(__xt);
	  __a = (__a + __lambda) / _Val(4);
	  __xt = (__xt + __lambda) / _Val(4);
	  __yt = (__yt + __lambda) / _Val(4);
	  __zt = (__zt + __lambda) / _Val(4);
	  __f *= _Val(4);
	  if (__q < __f * std::abs(__a))
	    {
	      _Tp __xf = (__a0 - __x) / (__f * __a);
	      _Tp __yf = (__a0 - __y) / (__f * __a);
	      _Tp __zf = -(__xf + __yf);
	      _Tp __e2 = __xf * __yf - __zf * __zf;
	      _Tp __e3 = __xf * __yf * __zf;
	      return (_Val(1)
		    - __e2 / _Val(10)
		    + __e3 / _Val(14)
		    + __e2 * __e2 / _Val(24)
		    - _Val(3) * __e2 * __e3 / _Val(44)) / std::sqrt(__a);
	    }
	}

      return _Tp{0};
    }

  /**
   *   @brief  Return the Carlson elliptic function
   *           @f$ R_C(x,y) = R_F(x,y,y) @f$ where @f$ R_F(x,y,z) @f$
   *           is the Carlson elliptic function of the first kind.
   *
   *   The Carlson elliptic function is defined by:
   *   @f[
   *       R_C(x,y) = \frac{1}{2} \int_0^\infty
   *                 \frac{dt}{(t + x)^{1/2}(t + y)}
   *   @f]
   *
   *   Based on Carlson's algorithms:
   *   -  B. C. Carlson Numer. Math. 33, 1 (1979)
   *   -  B. C. Carlson, Special Functions of Applied Mathematics (1977)
   *   -  Numerical Recipes in C, 2nd ed, pp. 261-269,
   *      by Press, Teukolsky, Vetterling, Flannery (1992)
   *
   *   @param  __x  The first argument.
   *   @param  __y  The second argument.
   *   @return  The Carlson elliptic function.
   */
  template<typename _Tp>
    _Tp
    __ellint_rc(_Tp __x, _Tp __y)
    {
      using _Val = __num_traits_t<_Tp>;
      if (std::imag(__y) == _Val(0) && std::real(__y) < _Val(0))
	return std::sqrt(__x / (__x - __y)) * __ellint_rc(__x - __y, -__y);
      const _Val __r = std::numeric_limits<_Val>::epsilon();
      _Tp __xt = __x;
      _Tp __yt = __y;
      _Tp __a0 = (__x + _Val(2) * __y) / _Val(3);
      _Val __q = std::pow( _Val(3) * __r, -_Val(1) / _Val(8) )
	       * std::abs(__a0 - __x);
      _Tp __a = __a0;
      _Val __f = _Val(1);

      while (true)
	{
	  _Tp __lambda = _Val(2) * std::sqrt(__xt) * std::sqrt(__yt) + __yt;
	  __a = (__a + __lambda) / _Val(4);
	  __xt = (__xt + __lambda) / _Val(4);
	  __yt = (__yt + __lambda) / _Val(4);
	  __f *= _Val(4);
	  if (__q < __f * std::abs(__a))
	    {
	      _Tp __s = (__y - __a0) / (__f * __a);
	      return (_Val(1) + __s * __s * (_Val(3) / _Val(10)
		    + __s * (_Val(1) / _Val(7)
		    + __s * (_Val(3) / _Val(8)
		    + __s * (_Val(9) / _Val(22)
		    + __s * (_Val(159) / _Val(208)
		    + __s * (_Val(9) / _Val(8)))))))) / std::sqrt(__a);
	    }
	}

      return _Tp{0};
    }

  /**
   *   @brief  Return the Carlson elliptic function @f$ R_J(x,y,z,p) @f$
   *           of the third kind.
   *
   *   The Carlson elliptic function of the third kind is defined by:
   *   @f[
   *       R_J(x,y,z,p) = \frac{3}{2} \int_0^\infty
   *       \frac{dt}{(t + x)^{1/2}(t + y)^{1/2}(t + z)^{1/2}(t + p)}
   *   @f]
   *
   *   Based on Carlson's algorithms:
   *   -  B. C. Carlson Numer. Math. 33, 1 (1979)
   *   -  B. C. Carlson, Special Functions of Applied Mathematics (1977)
   *   -  Numerical Recipes in C, 2nd ed, pp. 261-269,
   *      by Press, Teukolsky, Vetterling, Flannery (1992)
   *
   *   @param  __x  The first of three symmetric arguments.
   *   @param  __y  The second of three symmetric arguments.
   *   @param  __z  The third of three symmetric arguments.
   *   @param  __p  The fourth argument.
   *   @return  The Carlson elliptic function of the fourth kind.
   */
  template<typename _Tp>
    _Tp
    __ellint_rj(_Tp __x, _Tp __y, _Tp __z, _Tp __p)
    {
      using _Val = __num_traits_t<_Tp>;
      const _Val __r = std::numeric_limits<_Val>::epsilon();
      _Tp __xt = __x;
      _Tp __yt = __y;
      _Tp __zt = __z;
      _Tp __pt = __p;
      _Tp __a0 = (__x + __y + __z + _Val(2) * __p) / _Val(5);
      _Tp __delta = (__p - __x) * (__p - __y) * (__p - __z);
      _Val __q = std::pow(__r / _Val(4), -_Val(1) / _Val(6))
	       * std::max(std::abs(__a0 - __z),
			  std::max(std::abs(__a0 - __x),
				   std::max(std::abs(__a0 - __y), std::abs(__a0 - __p))));
      _Tp __a = __a0;
      _Val __f = _Val(1);
      _Val __fe = _Val(1);
      _Tp __sum = _Tp{};

      while (true)
	{
	  _Tp __xroot = std::sqrt(__xt);
	  _Tp __yroot = std::sqrt(__yt);
	  _Tp __zroot = std::sqrt(__zt);
	  _Tp __proot = std::sqrt(__pt);
	  _Tp __lambda = __xroot * __yroot
		       + __yroot * __zroot
		       + __zroot * __xroot;
	  __a = (__a + __lambda) / _Val(4);
	  __xt = (__xt + __lambda) / _Val(4);
	  __yt = (__yt + __lambda) / _Val(4);
	  __zt = (__zt + __lambda) / _Val(4);
	  __pt = (__pt + __lambda) / _Val(4);
	  _Tp __d = (__proot + __xroot) * (__proot + __yroot) * (__proot + __zroot);
	  _Tp __e = __delta / (__fe * __d * __d);
	  __sum += __ellint_rc(_Tp{1}, _Tp{1} + __e) / (__f * __d);
	  __f *= _Val(4);
	  __fe *= _Val(64);
	  if (__q < __f * std::abs(__a))
	    {
	      _Tp __xf = (__a0 - __x) / (__f * __a);
	      _Tp __yf = (__a0 - __y) / (__f * __a);
	      _Tp __zf = (__a0 - __z) / (__f * __a);
	      _Tp __xyz = __xf * __yf * __zf;
	      _Tp __pf = -(__xf + __yf + __zf) / _Val(2);
	      _Tp __pp = __pf * __pf;
	      _Tp __ppp = __pp * __pf;
	      _Tp __e2 = __xf * __yf + __yf * __zf + __zf * __xf - _Val(3) * __pp;
	      _Tp __e3 = __xyz + _Val(2) * __e2 * __pf + _Tp{4} * __ppp;
	      _Tp __e4 = (_Val(2) * __xyz + __e2 * __pf + _Val(3) * __ppp) * __pf;
	      _Tp __e5 = __xyz * __pp;
	      return (_Val(1) - _Val(3) * __e2 / _Val(14)
			      + __e3 / _Val(6)
			      + _Val(9) * __e2 * __e2 / _Val(88)
			      - _Val(3) * __e4 / _Val(22)
			      - _Val(9) * __e2 * __e3 / _Val(52)
			      + _Val(3) * __e5 / _Val(26)) / __f / __a / std::sqrt(__a)
			      + _Val(6) * __sum;
	    }
	}

      return _Tp{0};
    }

  /**
   *   @brief  Return the Carlson elliptic function of the second kind
   *           @f$ R_D(x,y,z) = R_J(x,y,z,z) @f$ where
   *           @f$ R_J(x,y,z,p) @f$ is the Carlson elliptic function
   *           of the third kind.
   *
   *   The Carlson elliptic function of the second kind is defined by:
   *   @f[
   *       R_D(x,y,z) = \frac{3}{2} \int_0^\infty
   *                 \frac{dt}{(t + x)^{1/2}(t + y)^{1/2}(t + z)^{3/2}}
   *   @f]
   *
   *   Based on Carlson's algorithms:
   *   -  B. C. Carlson Numer. Math. 33, 1 (1979)
   *   -  B. C. Carlson, Special Functions of Applied Mathematics (1977)
   *   -  Numerical Recipes in C, 2nd ed, pp. 261-269,
   *      by Press, Teukolsky, Vetterling, Flannery (1992)
   *
   *   @param  __x  The first of two symmetric arguments.
   *   @param  __y  The second of two symmetric arguments.
   *   @param  __z  The third argument.
   *   @return  The Carlson elliptic function of the second kind.
   */
  template<typename _Tp>
    _Tp
    __ellint_rd(_Tp __x, _Tp __y, _Tp __z)
    {
      using _Val = __num_traits_t<_Tp>;
      const _Val __r = std::numeric_limits<_Val>::epsilon();
      _Tp __xt = __x;
      _Tp __yt = __y;
      _Tp __zt = __z;
      _Tp __a0 = (__x + __y + _Val(3) * __z) / _Val(5);
      _Val __q = std::pow(__r / _Val(4), -_Val(1) / _Val(6))
	       * std::max(std::abs(__a0 - __z),
			  std::max(std::abs(__a0 - __x),
			  std::abs(__a0 - __y)));
      _Tp __a = __a0;
      _Val __f = _Val(1);
      _Tp __sum = _Tp{};

      while (true)
	{
	  _Tp __lambda = std::sqrt(__xt) * std::sqrt(__yt)
		       + std::sqrt(__yt) * std::sqrt(__zt)
		       + std::sqrt(__zt) * std::sqrt(__xt);
	  __sum += _Val(1) / __f / std::sqrt(__zt) / (__zt + __lambda);
	  __a = (__a + __lambda) / _Val(4);
	  __xt = (__xt + __lambda) / _Val(4);
	  __yt = (__yt + __lambda) / _Val(4);
	  __zt = (__zt + __lambda) / _Val(4);
	  __f *= _Val(4);
	  if (__q < __f * std::abs(__a))
	    {
	      _Tp __xf = (__a0 - __x) / (__f * __a);
	      _Tp __yf = (__a0 - __y) / (__f * __a);
	      _Tp __zf = -(__xf + __yf) / _Val(3);
	      _Tp __zz = __zf * __zf;
	      _Tp __xy = __xf * __yf;
	      _Tp __e2 = __xy - _Val(6) * __zz;
	      _Tp __e3 = (_Val(3) * __xy - _Val(8) * __zz) * __zf;
	      _Tp __e4 = _Val(3) * (__xy - __zz) * __zz;
	      _Tp __e5 = __xy * __zf * __zz;
	      return (_Val(1)
		    - _Val(3) * __e2 / _Val(14)
		    + __e3 / _Val(6)
		    + _Val(9) * __e2 * __e2 / _Val(88)
		    - _Val(3) * __e4 / _Val(22)
		    - _Val(9) * __e2 * __e3 / _Val(52)
		    + _Val(3) * __e5 / _Val(26)) / __f / __a / std::sqrt(__a)
		    + _Val(3) * __sum;
	    }
	}

      return _Tp{0};
    }

  template<typename _Tp>
    _Tp
    __comp_ellint_rf(_Tp __x, _Tp __y)
    {
      using _Val = __num_traits_t<_Tp>;
      const _Val __r = std::numeric_limits<_Val>::epsilon();
      const _Val __tolfact = _Val(2.7L) * std::sqrt(__r);
      __x = std::sqrt(__x);
      __y = std::sqrt(__y);
      while (true)
	{
	  _Tp __xt = __x;
	  __x = (__x + __y) / _Tp{2};
	  __y = std::sqrt(__xt) * std::sqrt(__y);
	  if (std::abs(__x - __y) < __tolfact * std::abs(__x))
	    return _Val(__gnu_cxx::__math_constants<_Tp>::__pi) / (__x + __y);
	}
    }

  template<typename _Tp>
    _Tp
    __comp_ellint_rg(_Tp __x, _Tp __y);

  /**
   *   @brief  Return the symmetric Carlson elliptic function of the second kind
   *           @f$ R_G(x,y,z) @f$.
   *
   *   The Carlson symmetric elliptic function of the second kind is defined by:
   *   @f[
   *       R_G(x,y,z) = \frac{1}{4} \int_0^\infty
   *                 dt t [(t + x)(t + y)(t + z)]^{-1/2}
   *                 (\frac{x}{t + x} + \frac{y}{t + y} + \frac{z}{t + z})
   *   @f]
   *
   *   Based on Carlson's algorithms:
   *   -  B. C. Carlson Numer. Math. 33, 1 (1979)
   *   -  B. C. Carlson, Special Functions of Applied Mathematics (1977)
   *   -  Numerical Recipes in C, 2nd ed, pp. 261-269,
   *      by Press, Teukolsky, Vetterling, Flannery (1992)
   *
   *   @param  __x  The first of three symmetric arguments.
   *   @param  __y  The second of three symmetric arguments.
   *   @param  __z  The third of three symmetric arguments.
   *   @return  The Carlson symmetric elliptic function of the second kind.
   */

  template<typename _Tp>
    _Tp
    __ellint_rg(_Tp __x, _Tp __y, _Tp __z)
    {
      using _Val = __num_traits_t<_Tp>;
      if (__z == _Tp{})
	{
	  if (__x == _Tp{})
	    return std::sqrt(__y);
	  else if (__y == _Tp{})
	    return std::sqrt(__x);
	  else
	    return __comp_ellint_rg(__x, __y);
	}
      else if (__x == _Tp{})
	{
	  if (__y == _Tp{})
	    return std::sqrt(__z);
	  else if (__z == _Tp{})
	    return std::sqrt(__y);
	  else
	    return __comp_ellint_rg(__y, __z);
	}
      else if (__y == _Tp{})
	{
	  if (__z == _Tp{})
	    return std::sqrt(__x);
	  else if (__x == _Tp{})
	    return std::sqrt(__z);
	  else
	    return __comp_ellint_rg(__z, __x);
	}
      else
	return (__z * __ellint_rf(__x, __y, __z)
	      - (__x - __z) * (__y - __z) * __ellint_rd(__x, __y, __z) / _Val(3)
	      + (std::sqrt(__x) * std::sqrt(__y) / std::sqrt(__z))) / _Val(2);
    }

  template<typename _Tp>
    _Tp
    __comp_ellint_rg(_Tp __x, _Tp __y)
    {
      using _Val = __num_traits_t<_Tp>;
      const _Val __r = std::numeric_limits<_Val>::epsilon();
      const _Val __tolfact = _Val(2.7L) * std::sqrt(__r);
      _Tp __xt = std::sqrt(__x);
      _Tp __yt = std::sqrt(__y);
      const _Tp __a = (__xt + __yt) / _Val(2);
      _Tp __sum = _Tp{};
      _Val __sf = _Val(1) / _Val(2);
      while (true)
	{
	  _Tp __xtt = __xt;
	  __xt = (__xt + __yt) / _Tp{2};
	  __yt = std::sqrt(__xtt) * std::sqrt(__yt);
	  _Tp __del = __xt - __yt;
	  if (std::abs(__del) < __tolfact * std::abs(__xt))
	    return (__a * __a - __sum)
		 * _Val(__gnu_cxx::__math_constants<_Tp>::__pi)
		 / (__xt + __yt) / _Val(2);
	  __sum += __sf * __del * __del;
	  __sf *= _Val(2);
	}
    }

  /**
   *   @brief  Return the complete elliptic integral of the first kind
   *           @f$ K(k) @f$ using the Carlson formulation.
   *
   *   The complete elliptic integral of the first kind is defined as
   *   @f[
   *     K(k) = F(k,\pi/2) = \int_0^{\pi/2}\frac{d\theta}
   *                                           {\sqrt{1 - k^2 sin^2\theta}}
   *   @f]
   *   where @f$ F(k,\phi) @f$ is the incomplete elliptic integral of the
   *   first kind.
   *
   *   @param  __k  The argument of the complete elliptic function.
   *   @return  The complete elliptic function of the first kind.
   */
  template<typename _Tp>
    _Tp
    __comp_ellint_1(_Tp __k)
    {
      constexpr auto _S_nan = std::numeric_limits<_Tp>::quiet_NaN();
      if (__isnan(__k))
	return _S_nan;
      else if (std::abs(__k) == _Tp{1})
	return _S_nan;
      else
	return __ellint_rf(_Tp{0}, _Tp{1} - __k * __k, _Tp{1});
    }

  /**
   *   @brief  Return the incomplete elliptic integral of the first kind
   *           @f$ F(k,\phi) @f$ using the Carlson formulation.
   *
   *   The incomplete elliptic integral of the first kind is defined as
   *   @f[
   *     F(k,\phi) = \int_0^{\phi}\frac{d\theta}
   *                                   {\sqrt{1 - k^2 sin^2\theta}}
   *   @f]
   *
   *   @param  __k  The argument of the elliptic function.
   *   @param  __phi  The integral limit argument of the elliptic function.
   *   @return  The elliptic function of the first kind.
   */
  template<typename _Tp>
    _Tp
    __ellint_1(_Tp __k, _Tp __phi)
    {
      constexpr auto _S_nan = __gnu_cxx::__math_constants<_Tp>::__NaN;
      constexpr auto _S_pi = __gnu_cxx::__math_constants<_Tp>::__pi;
      if (__isnan(__k) || __isnan(__phi))
	return _S_nan;
      else if (std::abs(__k) > _Tp{1})
	std::__throw_domain_error(__N("__ellint_1: bad argument"));
      else
	{
	  //  Reduce phi to -pi/2 < phi < +pi/2.
	  const int __n = std::floor(__phi / _S_pi + _Tp{0.5L});
	  const _Tp __phi_red = __phi - __n * _S_pi;

	  const _Tp __s = std::sin(__phi_red);
	  const _Tp __c = std::cos(__phi_red);

	  const _Tp __F = __s
			* __ellint_rf(__c * __c,
				_Tp{1} - __k * __k * __s * __s, _Tp{1});

	  if (__n == 0)
	    return __F;
	  else
	    return __F + _Tp{2} * __n * __comp_ellint_1(__k);
	}
    }

  /**
   *   @brief  Return the complete elliptic integral of the second kind
   *           @f$ E(k) @f$ using the Carlson formulation.
   *
   *   The complete elliptic integral of the second kind is defined as
   *   @f[
   *     E(k,\pi/2) = \int_0^{\pi/2}\sqrt{1 - k^2 sin^2\theta}
   *   @f]
   *
   *   @param  __k  The argument of the complete elliptic function.
   *   @return  The complete elliptic function of the second kind.
   */
  template<typename _Tp>
    _Tp
    __comp_ellint_2(_Tp __k)
    {
      constexpr auto _S_nan = std::numeric_limits<_Tp>::quiet_NaN();
      if (__isnan(__k))
	return _S_nan;
      else if (std::abs(__k) == 1)
	return _Tp{1};
      else if (std::abs(__k) > _Tp{1})
	std::__throw_domain_error(__N("__comp_ellint_2: bad argument"));
      else
	{
	  const _Tp __kk = __k * __k;

	  return __ellint_rf(_Tp{0}, _Tp{1} - __kk, _Tp{1})
	       - __kk * __ellint_rd(_Tp{0}, _Tp{1} - __kk, _Tp{1}) / _Tp{3};
	}
    }

  /**
   *   @brief  Return the incomplete elliptic integral of the second kind
   *           @f$ E(k,\phi) @f$ using the Carlson formulation.
   *
   *   The incomplete elliptic integral of the second kind is defined as
   *   @f[
   *     E(k,\phi) = \int_0^{\phi} \sqrt{1 - k^2 sin^2\theta}
   *   @f]
   *
   *   @param  __k  The argument of the elliptic function.
   *   @param  __phi  The integral limit argument of the elliptic function.
   *   @return  The elliptic function of the second kind.
   */
  template<typename _Tp>
    _Tp
    __ellint_2(_Tp __k, _Tp __phi)
    {
      constexpr auto _S_nan = __gnu_cxx::__math_constants<_Tp>::__NaN;
      constexpr auto _S_pi = __gnu_cxx::__math_constants<_Tp>::__pi;
      if (__isnan(__k) || __isnan(__phi))
	return _S_nan;
      else if (std::abs(__k) > _Tp{1})
	std::__throw_domain_error(__N("__ellint_2: bad argument"));
      else
	{
	  //  Reduce phi to -pi/2 < phi < +pi/2.
	  const int __n = std::floor(__phi / _S_pi + _Tp{0.5L});
	  const _Tp __phi_red = __phi - __n * _S_pi;

	  const _Tp __kk = __k * __k;
	  const _Tp __s = std::sin(__phi_red);
	  const _Tp __ss = __s * __s;
	  const _Tp __sss = __ss * __s;
	  const _Tp __c = std::cos(__phi_red);
	  const _Tp __cc = __c * __c;

	  const _Tp __E = __s
			* __ellint_rf(__cc, _Tp{1} - __kk * __ss, _Tp{1})
			- __kk * __sss
			* __ellint_rd(__cc, _Tp{1} - __kk * __ss, _Tp{1})
			/ _Tp{3};

	  if (__n == 0)
	    return __E;
	  else
	    return __E + _Tp{2} * __n * __comp_ellint_2(__k);
	}
    }

  /**
   *   @brief Return the complete elliptic integral of the third kind
   *          @f$ \Pi(k,\nu) = \Pi(k,\nu,\pi/2) @f$ using the
   *          Carlson formulation.
   *
   *   The complete elliptic integral of the third kind is defined as
   *   @f[
   *     \Pi(k,\nu) = \int_0^{\pi/2}
   *                   \frac{d\theta}
   *                 {(1 - \nu \sin^2\theta)\sqrt{1 - k^2 \sin^2\theta}}
   *   @f]
   *
   *   @param  __k  The argument of the elliptic function.
   *   @param  __nu  The second argument of the elliptic function.
   *   @return  The complete elliptic function of the third kind.
   */
  template<typename _Tp>
    _Tp
    __comp_ellint_3(_Tp __k, _Tp __nu)
    {
      constexpr auto _S_nan = __gnu_cxx::__math_constants<_Tp>::__NaN;
      if (__isnan(__k) || __isnan(__nu))
	return _S_nan;
      else if (__nu == _Tp{1})
	return __gnu_cxx::__math_constants<_Tp>::__inf;
      else if (std::abs(__k) > _Tp{1})
	std::__throw_domain_error(__N("__comp_ellint_3: bad argument"));
      else
	{
	  const _Tp __kk = __k * __k;

	  return __ellint_rf(_Tp{0}, _Tp{1} - __kk, _Tp{1})
	       - __nu
	       * __ellint_rj(_Tp{0}, _Tp{1} - __kk, _Tp{1}, _Tp{1} + __nu)
	       / _Tp{3};
	}
    }

  /**
   *   @brief Return the incomplete elliptic integral of the third kind
   *          @f$ \Pi(k,\nu,\phi) @f$ using the Carlson formulation.
   *
   *   The incomplete elliptic integral of the third kind is defined as
   *   @f[
   *     \Pi(k,\nu,\phi) = \int_0^{\phi}
   *                       \frac{d\theta}
   *                            {(1 - \nu \sin^2\theta)
   *                             \sqrt{1 - k^2 \sin^2\theta}}
   *   @f]
   *
   *   @param  __k  The argument of the elliptic function.
   *   @param  __nu  The second argument of the elliptic function.
   *   @param  __phi  The integral limit argument of the elliptic function.
   *   @return  The elliptic function of the third kind.
   */
  template<typename _Tp>
    _Tp
    __ellint_3(_Tp __k, _Tp __nu, _Tp __phi)
    {
      constexpr auto _S_nan = __gnu_cxx::__math_constants<_Tp>::__NaN;
      constexpr auto _S_pi = __gnu_cxx::__math_constants<_Tp>::__pi;
      if (__isnan(__k) || __isnan(__nu) || __isnan(__phi))
	return _S_nan;
      else if (std::abs(__k) > _Tp{1})
	std::__throw_domain_error(__N("__ellint_3: bad argument"));
      else
	{
	  //  Reduce phi to -pi/2 < phi < +pi/2.
	  const int __n = std::floor(__phi / _S_pi + _Tp{0.5L});
	  const _Tp __phi_red = __phi - __n * _S_pi;

	  const _Tp __kk = __k * __k;
	  const _Tp __s = std::sin(__phi_red);
	  const _Tp __ss = __s * __s;
	  const _Tp __sss = __ss * __s;
	  const _Tp __c = std::cos(__phi_red);
	  const _Tp __cc = __c * __c;

	  const _Tp __Pi = __s
			 * __ellint_rf(__cc, _Tp{1} - __kk * __ss, _Tp{1})
			 - __nu * __sss
			 * __ellint_rj(__cc, _Tp{1} - __kk * __ss, _Tp{1},
				       _Tp{1} + __nu * __ss) / _Tp{3};

	  if (__n == 0)
	    return __Pi;
	  else
	    return __Pi + _Tp{2} * __n * __comp_ellint_3(__k, __nu);
	}
    }

_GLIBCXX_END_NAMESPACE_VERSION
} // namespace __detail
} // namespace std

#endif // _GLIBCXX_BITS_SF_ELLINT_TCC

