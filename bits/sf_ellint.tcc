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
// (1)  B. C. Carlson Numer. Math. 33, 1 (1979)
// (2)  B. C. Carlson, Special Functions of Applied Mathematics (1977)
// (3)  The Gnu Scientific Library, http://www.gnu.org/software/gsl
// (4)  Numerical Recipes in C, 2nd ed, by W. H. Press, S. A. Teukolsky,
//      W. T. Vetterling, B. P. Flannery, Cambridge University Press
//      (1992), pp. 261-269
// (5)  Toshio Fukushima, Elliptic functions and elliptic integrals for
//      celestial mechanics and dynamical astronomy

#ifndef _GLIBCXX_BITS_SF_ELLINT_TCC
#define _GLIBCXX_BITS_SF_ELLINT_TCC 1

#pragma GCC system_header

#include <complex>
#include <ext/math_const.h>

namespace std _GLIBCXX_VISIBILITY(default)
{
// Implementation-space details.
namespace __detail
{
_GLIBCXX_BEGIN_NAMESPACE_VERSION

  /**
   * @brief  Return the Carlson elliptic function
   * 	     @f$ R_C(x,y) = R_F(x,y,y) @f$ where @f$ R_F(x,y,z) @f$
   * 	     is the Carlson elliptic function of the first kind.
   *
   * The Carlson elliptic function is defined by:
   * @f[
   * 	 R_C(x,y) = \frac{1}{2} \int_0^\infty
   * 		   \frac{dt}{(t + x)^{1/2}(t + y)}
   * @f]
   *
   * Based on Carlson's algorithms:
   * -  B. C. Carlson Numer. Math. 33, 1 (1979)
   * -  B. C. Carlson, Special Functions of Applied Mathematics (1977)
   * -  Numerical Recipes in C, 2nd ed, pp. 261-269,
   * 	by Press, Teukolsky, Vetterling, Flannery (1992)
   *
   * @param  __x  The first argument.
   * @param  __y  The second argument.
   * @return  The Carlson elliptic function.
   */
  template<typename _Tp>
    _Tp
    __ellint_rc(_Tp __x, _Tp __y)
    {
      using _Real = __num_traits_t<_Tp>;
      constexpr auto _S_NaN = __gnu_cxx::__quiet_NaN<_Real>();
      constexpr auto _S_min = __gnu_cxx::__min<_Real>();
      constexpr auto _S_max = __gnu_cxx::__max<_Real>();
      constexpr auto _S_eps = __gnu_cxx::__epsilon<_Real>();
      constexpr auto _S_lolim = _Real(5) * _S_min;
      constexpr auto _S_uplim = _S_max / _Real(5);

      if (__isnan(__x) || __isnan(__y))
	return _S_NaN;
      else if (std::imag(__x) == _Real{} && std::real(__x) < _Real{})
	std::__throw_domain_error(__N("__ellint_rc: argument less than zero"));
      else if (std::abs(__x) + std::abs(__y) < _S_lolim)
        std::__throw_domain_error(__N("__ellint_rc: arguments too small"));
      else if (std::imag(__y) == _Real{0} && std::real(__y) < _Real{0})
	{
	  if (std::abs(__x) == _Real{0})
	    return _Tp{};
	  else
	    return std::sqrt(__x / (__x - __y)) * __ellint_rc(__x - __y, -__y);
	}
      else
	{
	  auto __xt = __x;
	  auto __yt = __y;
	  auto _A0 = (__x + _Real{2} * __y) / _Real{3};
	  auto _Q = std::pow(_Real{3} * _S_eps, -_Real{1} / _Real{8})
		  * std::abs(_A0 - __x);
	  auto _A = _A0;
	  auto __f = _Real{1};

	  while (true)
	    {
	      auto __lambda = _Real{2} * std::sqrt(__xt) * std::sqrt(__yt)
			    + __yt;
	      _A = (_A + __lambda) / _Real{4};
	      __xt = (__xt + __lambda) / _Real{4};
	      __yt = (__yt + __lambda) / _Real{4};
	      __f *= _Real{4};
	      if (_Q < __f * std::abs(_A))
		{
		  auto __s = (__y - _A0) / (__f * _A);
		  return (_Real{1} + __s * __s * (_Real{3} / _Real{10}
			+ __s * (_Real{1} / _Real{7}
			+ __s * (_Real{3} / _Real{8}
			+ __s * (_Real{9} / _Real{22}
			+ __s * (_Real{159} / _Real{208}
			+ __s * (_Real{9} / _Real{8}))))))) / std::sqrt(_A);
		}
	    }

	  return _Tp{};
	}
    }

  /**
   * @brief  Return the Carlson elliptic function of the second kind
   * 	     @f$ R_D(x,y,z) = R_J(x,y,z,z) @f$ where
   * 	     @f$ R_J(x,y,z,p) @f$ is the Carlson elliptic function
   * 	     of the third kind.
   *
   * The Carlson elliptic function of the second kind is defined by:
   * @f[
   * 	 R_D(x,y,z) = \frac{3}{2} \int_0^\infty
   * 		   \frac{dt}{(t + x)^{1/2}(t + y)^{1/2}(t + z)^{3/2}}
   * @f]
   *
   * Based on Carlson's algorithms:
   * -  B. C. Carlson Numer. Math. 33, 1 (1979)
   * -  B. C. Carlson, Special Functions of Applied Mathematics (1977)
   * -  Numerical Recipes in C, 2nd ed, pp. 261-269,
   * 	by Press, Teukolsky, Vetterling, Flannery (1992)
   *
   * @param  __x  The first of two symmetric arguments.
   * @param  __y  The second of two symmetric arguments.
   * @param  __z  The third argument.
   * @return  The Carlson elliptic function of the second kind.
   */
  template<typename _Tp>
    _Tp
    __ellint_rd(_Tp __x, _Tp __y, _Tp __z)
    {
      using _Real = __num_traits_t<_Tp>;
      constexpr auto _S_NaN = __gnu_cxx::__quiet_NaN<_Real>();
      constexpr auto _S_min = __gnu_cxx::__min<_Real>();
      constexpr auto _S_max = __gnu_cxx::__max<_Real>();
      constexpr auto _S_eps = __gnu_cxx::__epsilon<_Real>();
      constexpr auto _S_lolim = _Real(5) * _S_min;
      constexpr auto _S_uplim = _S_max / _Real(5);

      if (__isnan(__x) || __isnan(__y) || __isnan(__z))
	return _S_NaN;
      else if ((std::imag(__x) == _Real{} && std::real(__x) < _Real{})
	    || (std::imag(__y) == _Real{} && std::real(__y) < _Real{})
	    || (std::imag(__z) == _Real{} && std::real(__z) < _Real{}))
        std::__throw_domain_error(__N("__ellint_rd: argument less than zero"));
      else if (std::abs(__x) + std::abs(__y) < _S_lolim
	    || std::abs(__z) < _S_lolim)
	std::__throw_domain_error(__N("__ellint_rd: arguments too small"));
      else
	{
	  auto __xt = __x;
	  auto __yt = __y;
	  auto __zt = __z;
	  auto _A0 = (__x + __y + _Real{3} * __z) / _Real{5};
	  auto _Q = std::pow(_S_eps / _Real{4}, -_Real{1} / _Real{6})
		  * std::max(std::abs(_A0 - __z),
			     std::max(std::abs(_A0 - __x),
			     std::abs(_A0 - __y)));
	  auto _A = _A0;
	  auto __f = _Real{1};
	  auto __sum = _Tp{};

	  while (true)
	    {
	      auto __lambda = std::sqrt(__xt) * std::sqrt(__yt)
			    + std::sqrt(__yt) * std::sqrt(__zt)
			    + std::sqrt(__zt) * std::sqrt(__xt);
	      __sum += _Real{1} / __f / std::sqrt(__zt) / (__zt + __lambda);
	      _A = (_A + __lambda) / _Real{4};
	      __xt = (__xt + __lambda) / _Real{4};
	      __yt = (__yt + __lambda) / _Real{4};
	      __zt = (__zt + __lambda) / _Real{4};
	      __f *= _Real{4};
	      if (_Q < __f * std::abs(_A))
		{
		  auto _Xi = (_A0 - __x) / (__f * _A);
		  auto _Yi = (_A0 - __y) / (__f * _A);
		  auto _Zi = -(_Xi + _Yi) / _Real{3};
		  auto _ZZ = _Zi * _Zi;
		  auto _XY = _Xi * _Yi;
		  auto _E2 = _XY - _Real{6} * _ZZ;
		  auto _E3 = (_Real{3} * _XY - _Real{8} * _ZZ) * _Zi;
		  auto _E4 = _Real{3} * (_XY - _ZZ) * _ZZ;
		  auto _E5 = _XY * _Zi * _ZZ;
		  return (_Real{1}
			- _Real{3} * _E2 / _Real{14}
			+ _E3 / _Real{6}
			+ _Real{9} * _E2 * _E2 / _Real{88}
			- _Real{3} * _E4 / _Real{22}
			- _Real{9} * _E2 * _E3 / _Real{52}
			+ _Real{3} * _E5 / _Real{26}) / __f / _A / std::sqrt(_A)
			+ _Real{3} * __sum;
		}
	    }

	  return _Tp{};
	}
    }

  template<typename _Tp>
    _Tp
    __comp_ellint_rf(_Tp __x, _Tp __y)
    {
      using _Real = __num_traits_t<_Tp>;
      constexpr auto _S_NaN = __gnu_cxx::__quiet_NaN<_Real>();
      constexpr auto _S_eps = __gnu_cxx::__epsilon<_Real>();
      constexpr auto _S_pi = __gnu_cxx::__math_constants<_Real>::__pi;
      const auto _S_tolfact = _Real{2.7L} * __gnu_cxx::__sqrt_eps<_Real>();

      if (__isnan(__x) || __isnan(__y))
	return _S_NaN;
      else
	{
	  __x = std::sqrt(__x);
	  __y = std::sqrt(__y);
	  while (true)
	    {
	      auto __xt = __x;
	      __x = (__x + __y) / _Tp{2};
	      __y = std::sqrt(__xt) * std::sqrt(__y);
	      if (std::abs(__x - __y) < _S_tolfact * std::abs(__x))
		return _S_pi / (__x + __y);
	    }
	}
    }

  /**
   * @brief Return the Carlson elliptic function @f$ R_F(x,y,z) @f$
   * 	    of the first kind.
   *
   * The Carlson elliptic function of the first kind is defined by:
   * @f[
   * 	 R_F(x,y,z) = \frac{1}{2} \int_0^\infty
   * 		   \frac{dt}{(t + x)^{1/2}(t + y)^{1/2}(t + z)^{1/2}}
   * @f]
   *
   * @param  __x  The first of three symmetric arguments.
   * @param  __y  The second of three symmetric arguments.
   * @param  __z  The third of three symmetric arguments.
   * @return  The Carlson elliptic function of the first kind.
   */
  template<typename _Tp>
    _Tp
    __ellint_rf(_Tp __x, _Tp __y, _Tp __z)
    {
      using _Real = __num_traits_t<_Tp>;
      constexpr auto _S_NaN = __gnu_cxx::__quiet_NaN<_Real>();
      constexpr auto _S_min = __gnu_cxx::__min<_Real>();
      constexpr auto _S_max = __gnu_cxx::__max<_Real>();
      constexpr auto _S_eps = __gnu_cxx::__epsilon<_Real>();
      constexpr auto _S_lolim = _Real(5) * _S_min;
      constexpr auto _S_uplim = _S_max / _Real(5);

      if (__isnan(__x) || __isnan(__y) || __isnan(__z))
	return _S_NaN;
      else if (std::imag(__x) == _Real{} && std::real(__x) < _Real{}
	    || std::imag(__y) == _Real{} && std::real(__y) < _Real{}
	    || std::imag(__z) == _Real{} && std::real(__z) < _Real{})
        std::__throw_domain_error(__N("__ellint_rf: argument less than zero"));
      else if (std::abs(__x) + std::abs(__y) < _S_lolim
	    || std::abs(__x) + std::abs(__z) < _S_lolim
	    || std::abs(__y) + std::abs(__z) < _S_lolim)
        std::__throw_domain_error(__N("Argument too small in __ellint_rf"));

      if (std::abs(__z) < _S_eps)
        return __comp_ellint_rf(__x, __y);
      else if (std::abs(__z - __y) < _S_eps)
	return __ellint_rc(__x, __y);
      else
	{
	  auto __xt = __x;
	  auto __yt = __y;
	  auto __zt = __z;
	  auto _A0 = (__x + __y + __z) / _Real{3};
	  auto _Q = std::pow(_Real{3} * _S_eps, -_Real{1} / _Real{6})
		  * std::max(std::abs(_A0 - __z),
			     std::max(std::abs(_A0 - __x),
				      std::abs(_A0 - __y)));
	  auto _A = _A0;
	  auto __f = _Real{1};

	  while (true)
	    {
	      auto __lambda = std::sqrt(__xt) * std::sqrt(__yt)
			    + std::sqrt(__yt) * std::sqrt(__zt)
			    + std::sqrt(__zt) * std::sqrt(__xt);
	      _A = (_A + __lambda) / _Real{4};
	      __xt = (__xt + __lambda) / _Real{4};
	      __yt = (__yt + __lambda) / _Real{4};
	      __zt = (__zt + __lambda) / _Real{4};
	      __f *= _Real{4};
	      if (_Q < __f * std::abs(_A))
		{
		  auto _Xi = (_A0 - __x) / (__f * _A);
		  auto _Yi = (_A0 - __y) / (__f * _A);
		  auto _Zi = -(_Xi + _Yi);
		  auto _E2 = _Xi * _Yi - _Zi * _Zi;
		  auto _E3 = _Xi * _Yi * _Zi;
		  return (_Real{1}
			- _E2 / _Real{10}
			+ _E3 / _Real{14}
			+ _E2 * _E2 / _Real{24}
			- _Real{3} * _E2 * _E3 / _Real{44}) / std::sqrt(_A);
		}
	    }

	  return _Tp{};
	}
    }

  template<typename _Tp>
    _Tp
    __comp_ellint_rg(_Tp __x, _Tp __y)
    {
      using _Real = __num_traits_t<_Tp>;
      constexpr auto _S_NaN = __gnu_cxx::__quiet_NaN<_Real>();
      constexpr auto _S_eps = __gnu_cxx::__epsilon<_Real>();
      constexpr auto _S_pi = __gnu_cxx::__math_constants<_Real>::__pi;
      const auto _S_tolfact = _Real{2.7L} * __gnu_cxx::__sqrt_eps<_Real>();

      if (__isnan(__x) || __isnan(__y))
	return _S_NaN;
      else if (__x == _Tp{} && __y == _Tp{})
	return _Tp{};
      else if (__x == _Tp{})
	return std::sqrt(__y) / _Real{2};
      else if (__y == _Tp{})
	return std::sqrt(__x) / _Real{2};
      else
	{
	  auto __xt = std::sqrt(__x);
	  auto __yt = std::sqrt(__y);
	  const auto _A = (__xt + __yt) / _Real{2};
	  auto __sum = _Tp{};
	  auto __sf = _Real{1} / _Real{2};
	  while (true)
	    {
	      auto __xtt = __xt;
	      __xt = (__xt + __yt) / _Tp{2};
	      __yt = std::sqrt(__xtt) * std::sqrt(__yt);
	      auto __del = __xt - __yt;
	      if (std::abs(__del) < _S_tolfact * std::abs(__xt))
		return (_A * _A - __sum) * _S_pi / (__xt + __yt) / _Real{2};
	      __sum += __sf * __del * __del;
	      __sf *= _Real{2};
	    }
	}
    }

  /**
   * @brief  Return the symmetric Carlson elliptic function of the second kind
   * 	     @f$ R_G(x,y,z) @f$.
   *
   * The Carlson symmetric elliptic function of the second kind is defined by:
   * @f[
   * 	 R_G(x,y,z) = \frac{1}{4} \int_0^\infty
   * 		   dt t [(t + x)(t + y)(t + z)]^{-1/2}
   * 		   (\frac{x}{t + x} + \frac{y}{t + y} + \frac{z}{t + z})
   * @f]
   *
   * Based on Carlson's algorithms:
   * -  B. C. Carlson Numer. Math. 33, 1 (1979)
   * -  B. C. Carlson, Special Functions of Applied Mathematics (1977)
   * -  Numerical Recipes in C, 2nd ed, pp. 261-269,
   * 	by Press, Teukolsky, Vetterling, Flannery (1992)
   *
   * @param  __x  The first of three symmetric arguments.
   * @param  __y  The second of three symmetric arguments.
   * @param  __z  The third of three symmetric arguments.
   * @return  The Carlson symmetric elliptic function of the second kind.
   */

  template<typename _Tp>
    _Tp
    __ellint_rg(_Tp __x, _Tp __y, _Tp __z)
    {
      using _Real = __num_traits_t<_Tp>;
      constexpr auto _S_NaN = __gnu_cxx::__quiet_NaN<_Real>();

      if (__isnan(__x) || __isnan(__y) || __isnan(__z))
	return _S_NaN;
      else if (__z == _Tp{})
	return __comp_ellint_rg(__x, __y);
      else if (__x == _Tp{})
	return __comp_ellint_rg(__y, __z);
      else if (__y == _Tp{})
	return __comp_ellint_rg(__z, __x);
      else
	//return (__z * __ellint_rf(__x, __y, __z)
	//     - (__x - __z) * (__y - __z) * __ellint_rd(__x, __y, __z) / _Real{3}
	//     + (std::sqrt(__x) * std::sqrt(__y) / std::sqrt(__z))) / _Real{2};
	// There is a symmetric version that is less subject to cancellation loss
	// when the arguments are real:
	return (__x * (__y + __z) * __ellint_rd(__y, __z, __x)
	      + __y * (__z + __x) * __ellint_rd(__z, __x, __y)
	      + __z * (__x + __y) * __ellint_rd(__x, __y, __z)) / _Tp{6};
    }

  /**
   * @brief  Return the Carlson elliptic function @f$ R_J(x,y,z,p) @f$
   * 	     of the third kind.
   *
   * The Carlson elliptic function of the third kind is defined by:
   * @f[
   * 	 R_J(x,y,z,p) = \frac{3}{2} \int_0^\infty
   * 	 \frac{dt}{(t + x)^{1/2}(t + y)^{1/2}(t + z)^{1/2}(t + p)}
   * @f]
   *
   * Based on Carlson's algorithms:
   * -  B. C. Carlson Numer. Math. 33, 1 (1979)
   * -  B. C. Carlson, Special Functions of Applied Mathematics (1977)
   * -  Numerical Recipes in C, 2nd ed, pp. 261-269,
   * 	by Press, Teukolsky, Vetterling, Flannery (1992)
   *
   * @param  __x  The first of three symmetric arguments.
   * @param  __y  The second of three symmetric arguments.
   * @param  __z  The third of three symmetric arguments.
   * @param  __p  The fourth argument.
   * @return  The Carlson elliptic function of the fourth kind.
   */
  template<typename _Tp>
    _Tp
    __ellint_rj(_Tp __x, _Tp __y, _Tp __z, _Tp __p)
    {
      using _Real = __num_traits_t<_Tp>;
      constexpr auto _S_NaN = __gnu_cxx::__quiet_NaN<_Real>();
      constexpr auto _S_min = __gnu_cxx::__min<_Real>();
      constexpr auto _S_max = __gnu_cxx::__max<_Real>();
      constexpr auto _S_eps = __gnu_cxx::__epsilon<_Real>();
      constexpr auto _S_lolim = _Real(5) * _S_min;
      constexpr auto _S_uplim = _S_max / _Real(5);

      if (__isnan(__x) || __isnan(__y) || __isnan(__z) || __isnan(__p))
	return _S_NaN;
      else if (std::imag(__x) == _Real{} && std::real(__x) < _Real{}
	    || std::imag(__y) == _Real{} && std::real(__y) < _Real{}
	    || std::imag(__z) == _Real{} && std::real(__z) < _Real{})
        std::__throw_domain_error(__N("__ellint_rj: argument less than zero"));
      else if (std::abs(__x) + std::abs(__y) < _S_lolim
	    || std::abs(__x) + std::abs(__z) < _S_lolim
	    || std::abs(__y) + std::abs(__z) < _S_lolim
	    || std::abs(__p) < _S_lolim)
        std::__throw_domain_error(__N("__ellint_rj: argument too small"));
      else if (std::abs(__p - __z) < _S_eps)
	return __ellint_rd(__x, __y, __z);
      else
	{
	  auto __xt = __x;
	  auto __yt = __y;
	  auto __zt = __z;
	  auto __pt = __p;
	  auto _A0 = (__x + __y + __z + _Real{2} * __p) / _Real{5};
	  auto __delta = (__p - __x) * (__p - __y) * (__p - __z);
	  auto _Q = std::pow(_S_eps / _Real{4}, -_Real{1} / _Real{6})
		  * std::max(std::abs(_A0 - __z),
		      std::max(std::abs(_A0 - __x),
			std::max(std::abs(_A0 - __y), std::abs(_A0 - __p))));
	  auto _A = _A0;
	  auto __f = _Real{1};
	  auto __fe = _Real{1};
	  auto __sum = _Tp{};

	  while (true)
	    {
	      auto __xroot = std::sqrt(__xt);
	      auto __yroot = std::sqrt(__yt);
	      auto __zroot = std::sqrt(__zt);
	      auto __proot = std::sqrt(__pt);
	      auto __lambda = __xroot * __yroot
			    + __yroot * __zroot
			    + __zroot * __xroot;
	      _A = (_A + __lambda) / _Real{4};
	      __xt = (__xt + __lambda) / _Real{4};
	      __yt = (__yt + __lambda) / _Real{4};
	      __zt = (__zt + __lambda) / _Real{4};
	      __pt = (__pt + __lambda) / _Real{4};
	      auto __d = (__proot + __xroot)
		       * (__proot + __yroot)
		       * (__proot + __zroot);
	      auto _E = __delta / (__fe * __d * __d);
	      __sum += __ellint_rc(_Tp{1}, _Tp{1} + _E) / (__f * __d);
	      __f *= _Real{4};
	      __fe *= _Real{64};
	      if (_Q < __f * std::abs(_A))
		{
		  auto _Xi = (_A0 - __x) / (__f * _A);
		  auto _Yi = (_A0 - __y) / (__f * _A);
		  auto _Zi = (_A0 - __z) / (__f * _A);
		  auto _XYZ = _Xi * _Yi * _Zi;
		  auto _Pi = -(_Xi + _Yi + _Zi) / _Real{2};
		  auto _PP = _Pi * _Pi;
		  auto _PPP = _PP * _Pi;
		  auto _E2 = _Xi * _Yi
			   + _Yi * _Zi
			   + _Zi * _Xi
			   - _Real{3} * _PP;
		  auto _E3 = _XYZ + _Real{2} * _E2 * _Pi + _Tp{4} * _PPP;
		  auto _E4 = _Pi
			   * (_Real{2} * _XYZ + _E2 * _Pi + _Real{3} * _PPP);
		  auto _E5 = _XYZ * _PP;
		  return (_Real{1} - _Real{3} * _E2 / _Real{14}
			+ _E3 / _Real{6}
			+ _Real{9} * _E2 * _E2 / _Real{88}
			- _Real{3} * _E4 / _Real{22}
			- _Real{9} * _E2 * _E3 / _Real{52}
			+ _Real{3} * _E5 / _Real{26}) / __f / _A / std::sqrt(_A)
			+ _Real{6} * __sum;
		}
	    }

	  return _Tp{};
	}
    }

  /**
   * @brief  Return the complete elliptic integral of the first kind
   * @f$ K(k) @f$ using the Carlson formulation.
   *
   * The complete elliptic integral of the first kind is defined as
   * @f[
   *   K(k) = F(k,\pi/2) = \int_0^{\pi/2}\frac{d\theta}
   * 					     {\sqrt{1 - k^2 sin^2\theta}}
   * @f]
   * where @f$ F(k,\phi) @f$ is the incomplete elliptic integral of the
   * first kind.
   *
   * @param  __k  The modulus of the complete elliptic function.
   * @return  The complete elliptic function of the first kind.
   */
  template<typename _Tp>
    _Tp
    __comp_ellint_1(_Tp __k)
    {
      using _Real = __num_traits_t<_Tp>;
      constexpr auto _S_NaN = __gnu_cxx::__quiet_NaN<_Real>();

      if (__isnan(__k))
	return _S_NaN;
      else if (std::abs(__k) == _Real{1})
	return _S_NaN;
      else
	return __comp_ellint_rf(_Tp{1} - __k * __k, _Tp{1});
    }

  /**
   * @brief  Return the incomplete elliptic integral of the first kind
   * @f$ F(k,\phi) @f$ using the Carlson formulation.
   *
   * The incomplete elliptic integral of the first kind is defined as
   * @f[
   *   F(k,\phi) = \int_0^{\phi}\frac{d\theta}
   * 				     {\sqrt{1 - k^2 sin^2\theta}}
   * @f]
   *
   * @param  __k  The argument of the elliptic function.
   * @param  __phi  The integral limit argument of the elliptic function.
   * @return  The elliptic function of the first kind.
   */
  template<typename _Tp>
    _Tp
    __ellint_1(_Tp __k, _Tp __phi)
    {
      using _Real = __num_traits_t<_Tp>;
      constexpr auto _S_NaN = __gnu_cxx::__quiet_NaN<_Real>();
      constexpr auto _S_pi = __gnu_cxx::__math_constants<_Real>::__pi;

      if (__isnan(__k) || __isnan(__phi))
	return _S_NaN;
      else if (std::abs(__k) > _Real{1})
	std::__throw_domain_error(__N("__ellint_1: bad argument"));
      else
	{
	  // Reduce phi to -pi/2 < phi < +pi/2.
	  const int __n = std::floor(__phi / _S_pi + _Real{0.5L});
	  const auto __phi_red = __phi - __n * _S_pi;

	  const auto __s = std::sin(__phi_red);
	  const auto __c = std::cos(__phi_red);

	  const auto __F = __s
			 * __ellint_rf(__c * __c,
				       _Real{1} - __k * __k * __s * __s, _Tp{1});

	  if (__n == 0)
	    return __F;
	  else
	    return __F + _Tp{2} * __n * __comp_ellint_1(__k);
	}
    }

  /**
   * @brief  Return the complete elliptic integral of the second kind
   * 	     @f$ E(k) @f$ using the Carlson formulation.
   *
   * The complete elliptic integral of the second kind is defined as
   * @f[
   *   E(k,\pi/2) = \int_0^{\pi/2}\sqrt{1 - k^2 sin^2\theta}
   * @f]
   *
   * @param  __k  The modulus of the complete elliptic function.
   * @return  The complete elliptic function of the second kind.
   */
  template<typename _Tp>
    _Tp
    __comp_ellint_2(_Tp __k)
    {
      using _Real = __num_traits_t<_Tp>;
      constexpr auto _S_NaN = __gnu_cxx::__quiet_NaN<_Real>();

      if (__isnan(__k))
	return _S_NaN;
      else if (std::abs(__k) == _Real{1})
	return _Tp{1};
      else if (std::abs(__k) > _Real{1})
	std::__throw_domain_error(__N("__comp_ellint_2: bad argument"));
      else
	{
	  const auto __kk = __k * __k;

	  return __ellint_rf(_Real{0}, _Real{1} - __kk, _Real{1})
	       - __kk * __ellint_rd(_Real{0}, _Real{1} - __kk, _Real{1}) / _Tp{3};
	}
    }

  /**
   * @brief  Return the incomplete elliptic integral of the second kind
   * @f$ E(k,\phi) @f$ using the Carlson formulation.
   *
   * The incomplete elliptic integral of the second kind is defined as
   * @f[
   *   E(k,\phi) = \int_0^{\phi} \sqrt{1 - k^2 sin^2\theta}
   * @f]
   *
   * @param  __k  The argument of the elliptic function.
   * @param  __phi  The integral limit argument of the elliptic function.
   * @return  The elliptic function of the second kind.
   */
  template<typename _Tp>
    _Tp
    __ellint_2(_Tp __k, _Tp __phi)
    {
      using _Real = __num_traits_t<_Tp>;
      constexpr auto _S_NaN = __gnu_cxx::__quiet_NaN<_Real>();
      constexpr auto _S_pi = __gnu_cxx::__math_constants<_Real>::__pi;

      if (__isnan(__k) || __isnan(__phi))
	return _S_NaN;
      else if (std::abs(__k) > _Real{1})
	std::__throw_domain_error(__N("__ellint_2: bad argument"));
      else
	{
	  // Reduce phi to -pi/2 < phi < +pi/2.
	  const int __n = std::floor(std::abs(__phi) / _S_pi + _Real{0.5L});
	  const auto __phi_red = __phi - __n * _S_pi;

	  const auto __kk = __k * __k;
	  const auto __s = std::sin(__phi_red);
	  const auto __ss = __s * __s;
	  const auto __sss = __ss * __s;
	  const auto __c = std::cos(__phi_red);
	  const auto __cc = __c * __c;

	  const auto _E = __s
		        * __ellint_rf(__cc, _Tp{1} - __kk * __ss, _Tp{1})
		        - __kk * __sss
		        * __ellint_rd(__cc, _Tp{1} - __kk * __ss, _Tp{1})
		        / _Tp{3};

	  if (__n == 0)
	    return _E;
	  else
	    return _E + _Tp{2} * __n * __comp_ellint_2(__k);
	}
    }

  /**
   * @brief Return the complete elliptic integral of the third kind
   * @f$ \Pi(k,\nu) = \Pi(k,\nu,\pi/2) @f$ using the
   * Carlson formulation.
   *
   * The complete elliptic integral of the third kind is defined as
   * @f[
   *   \Pi(k,\nu) = \int_0^{\pi/2}
   * 		     \frac{d\theta}
   * 		   {(1 - \nu \sin^2\theta)\sqrt{1 - k^2 \sin^2\theta}}
   * @f]
   *
   * @param  __k  The argument of the elliptic function.
   * @param  __nu  The second argument of the elliptic function.
   * @return  The complete elliptic function of the third kind.
   */
  template<typename _Tp>
    _Tp
    __comp_ellint_3(_Tp __k, _Tp __nu)
    {
      using _Real = __num_traits_t<_Tp>;
      constexpr auto _S_NaN = __gnu_cxx::__quiet_NaN<_Real>();

      if (__isnan(__k) || __isnan(__nu))
	return _S_NaN;
      else if (__nu == _Tp{1})
	return __gnu_cxx::__infinity<_Real>();
      else if (std::abs(__k) > _Real{1})
	std::__throw_domain_error(__N("__comp_ellint_3: bad argument"));
      else
	{
	  const auto __kk = __k * __k;

	  return __ellint_rf(_Tp{0}, _Tp{1} - __kk, _Tp{1})
	       - __nu
	       * __ellint_rj(_Tp{0}, _Tp{1} - __kk, _Tp{1}, _Tp{1} + __nu)
	       / _Tp{3};
	}
    }

  /**
   * @brief Return the incomplete elliptic integral of the third kind
   * @f$ \Pi(k,\nu,\phi) @f$ using the Carlson formulation.
   *
   * The incomplete elliptic integral of the third kind is defined as
   * @f[
   *   \Pi(k,\nu,\phi) = \int_0^{\phi}
   * 			 \frac{d\theta}
   * 			      {(1 - \nu \sin^2\theta)
   * 			       \sqrt{1 - k^2 \sin^2\theta}}
   * @f]
   *
   * @param  __k  The argument of the elliptic function.
   * @param  __nu  The second argument of the elliptic function.
   * @param  __phi  The integral limit argument of the elliptic function.
   * @return  The elliptic function of the third kind.
   */
  template<typename _Tp>
    _Tp
    __ellint_3(_Tp __k, _Tp __nu, _Tp __phi)
    {
      using _Real = __num_traits_t<_Tp>;
      constexpr auto _S_NaN = __gnu_cxx::__quiet_NaN<_Real>();
      constexpr auto _S_pi = __gnu_cxx::__math_constants<_Real>::__pi;

      if (__isnan(__k) || __isnan(__nu) || __isnan(__phi))
	return _S_NaN;
      else if (std::abs(__k) > _Real{1})
	std::__throw_domain_error(__N("__ellint_3: bad argument"));
      else
	{
	  // Reduce phi to -pi/2 < phi < +pi/2.
	  const int __n = std::floor(std::real(__phi) / _S_pi + _Real{0.5L});
	  const auto __phi_red = __phi - __n * _S_pi;

	  const auto __kk = __k * __k;
	  const auto __s = std::sin(__phi_red);
	  const auto __ss = __s * __s;
	  const auto __sss = __ss * __s;
	  const auto __c = std::cos(__phi_red);
	  const auto __cc = __c * __c;

	  const auto _Pi = __s
			 * __ellint_rf(__cc, _Tp{1} - __kk * __ss, _Tp{1})
			 - __nu * __sss
			 * __ellint_rj(__cc, _Tp{1} - __kk * __ss, _Tp{1},
				       _Tp{1} + __nu * __ss) / _Tp{3};

	  if (__n == 0)
	    return _Pi;
	  else
	    return _Pi + _Tp{2} * __n * __comp_ellint_3(__k, __nu);
	}
    }

  /**
   * Return the Legendre elliptic integral D.
   */
  template<typename _Tp>
    _Tp
    __ellint_d(_Tp __k, _Tp __phi)
    {
      using _Real = __num_traits_t<_Tp>;
      constexpr auto _S_NaN = __gnu_cxx::__quiet_NaN<_Real>();

      if (__isnan(__k) || __isnan(__phi))
	return _S_NaN;
      else if (std::abs(__k) > _Real{1})
	std::__throw_domain_error(__N("__ellint_d: bad argument"));
      else
	{
	  auto __sinphi = std::sin(__phi);
	  auto __sinphi2 = __sinphi * __sinphi;
	  auto __k2 = __k * __k;
	  auto __arg1 = _Tp{1} - __sinphi2;
	  auto __arg2 = _Tp{1} - __k2 * __sinphi2;
	  return __sinphi * __sinphi2 * __ellint_rd(__arg1, __arg2, _Tp{1}) / _Tp{3};
	}
    }

  /**
   * Return the complete Legendre elliptic integral D.
   */
  template<typename _Tp>
    _Tp
    __comp_ellint_d(_Tp __k)
    {
      using _Real = __num_traits_t<_Tp>;
      constexpr auto _S_NaN = __gnu_cxx::__quiet_NaN<_Real>();

      if (__isnan(__k))
	return _S_NaN;
      else
	return __ellint_rd(_Tp{0}, _Tp{1} - __k * __k, _Tp{1}) / _Tp{3};
    }

  /**
   * Return the Bulirsch elliptic integrals of the first kind.
   */
  template<typename _Tp>
    _Tp
    __ellint_el1(_Tp __x, _Tp __k_c)
    {
      using _Real = __num_traits_t<_Tp>;
      constexpr auto _S_NaN = __gnu_cxx::__quiet_NaN<_Real>();

      if (__isnan(__x) || __isnan(__k_c))
	return _S_NaN;
      else
	{
	  auto __x2 = __x * __x;
	  auto __k2_c = __k_c * __k_c;
	  auto __arg2 = _Tp{1} + __k2_c * __x2;
	  return __x * __ellint_rf(_Tp{1}, __arg2, _Tp{1} + __x2);
	}
    }

  /**
   * Return the Bulirsch elliptic integrals of the second kind.
   */
  template<typename _Tp>
    _Tp
    __ellint_el2(_Tp __x, _Tp __k_c, _Tp __a, _Tp __b)
    {
      using _Real = __num_traits_t<_Tp>;
      constexpr auto _S_NaN = __gnu_cxx::__quiet_NaN<_Real>();

      if (__isnan(__x) || __isnan(__k_c) || __isnan(__a) || __isnan(__b))
	return _S_NaN;
      else
	{
	  auto __x2 = __x * __x;
	  auto __x3 = __x * __x;
	  auto __k2_c = __k_c * __k_c;
	  auto __arg2 = _Tp{1} + __k2_c * __x2;
	  auto __arg3 = _Tp{1} + __x2;
	  return __a * __x * __ellint_rf(_Tp{1}, __arg2, __arg3)
               + (__b - __a) * __x3
		 * __ellint_rd(_Tp{1}, __arg2, __arg3) / _Tp{3};
	}
    }

  /**
   * Return the Bulirsch elliptic integrals of the third kind.
   */
  template<typename _Tp>
    _Tp
    __ellint_el3(_Tp __x, _Tp __k_c, _Tp __p)
    {
      using _Real = __num_traits_t<_Tp>;
      constexpr auto _S_NaN = __gnu_cxx::__quiet_NaN<_Real>();

      if (__isnan(__x) || __isnan(__k_c) || __isnan(__p))
	return _S_NaN;
      else
	{
	  auto __x2 = __x * __x;
	  auto __x3 = __x * __x;
	  auto __k2_c = __k_c * __k_c;
	  auto __arg2 = _Tp{1} + __k2_c * __x2;
	  auto __arg3 = _Tp{1} + __x2;
	  auto __arg4 = _Tp{1} + __p * __x2;
	  return __x * __ellint_rf(_Tp{1}, __arg2, __arg3)
	       + (_Tp{1} - __p) * __x3
		 * __ellint_rj(_Tp{1}, __arg2, __arg3, __arg4) / _Tp{3};
	}
    }

  /**
   * Return the Bulirsch complete elliptic integrals.
   */
  template<typename _Tp>
    _Tp
    __ellint_cel(_Tp __k_c, _Tp __p, _Tp __a, _Tp __b)
    {
      using _Real = __num_traits_t<_Tp>;
      constexpr auto _S_NaN = __gnu_cxx::__quiet_NaN<_Real>();

      if (__isnan(__k_c) || __isnan(__p) || __isnan(__a) || __isnan(__b))
	return _S_NaN;
      else
	{
	  auto __k2_c = __k_c * __k_c;
	  return __a * __ellint_rf(_Tp{0}, __k2_c, _Tp{1})
	       + (__b - __p * __a)
		 * __ellint_rj(_Tp{0}, __k2_c, _Tp{1}, __p) / _Tp{3};
	}
    }

  /**
   * Return the Jacobi zeta function.
   */
  template<typename _Tp>
    _Tp
    __jacobi_zeta(_Tp __k, _Tp __phi)
    {
      using _Real = __num_traits_t<_Tp>;
      constexpr auto _S_NaN = __gnu_cxx::__quiet_NaN<_Real>();
      constexpr auto _S_pi_2 = __gnu_cxx::__math_constants<_Real>::__pi_half;
      constexpr auto _S_eps = __gnu_cxx::__epsilon<_Real>();

      if (__isnan(__k) || __isnan(__phi))
	return _S_NaN;
      else if (std::abs(__k) > _Real{1})
	std::__throw_domain_error(__N("__jacobi_zeta: bad argument"));
      else if (std::abs(__k) < _S_eps)
	return _Tp{0};
      else if (std::abs(__k - _Tp{1}) < _S_eps)
	return std::sin(__phi);
      else if (std::real(__phi) < _Real{0})
	return -__jacobi_zeta(__k, -__phi);
      else if (std::abs(__phi) < _S_eps || std::abs(__phi - _S_pi_2) < _S_eps)
	return _Tp{0};
      else
	{
	  auto __mc = __k * __k;
	  auto __cosphi = std::cos(__phi);
	  auto __sinphi = std::sin(__phi);
	  auto __m = _Tp{1} - __mc;
	  auto __arg4 = _Tp{1} - __mc  * __sinphi * __sinphi;
	  return __mc * __cosphi * __sinphi * std::sqrt(__arg4)
	       * __ellint_rj(_Tp{0}, __m, _Tp{1}, __arg4)
	       / (_Tp{3} * __comp_ellint_1(__k));
	}
    }

  /**
   * Return the Heuman lambda function.
   */
  template<typename _Tp>
    _Tp
    __heuman_lambda(_Tp __k, _Tp __phi)
    {
      using _Real = __num_traits_t<_Tp>;
      constexpr auto _S_NaN = __gnu_cxx::__quiet_NaN<_Real>();
      constexpr auto _S_pi = __gnu_cxx::__math_constants<_Real>::__pi;
      constexpr auto _S_pi_2 = __gnu_cxx::__math_constants<_Real>::__pi_half;
      constexpr auto _S_eps = __gnu_cxx::__epsilon<_Real>();

      if (__isnan(__k) || __isnan(__phi))
	return _S_NaN;
      else if (std::abs(__k) > _Tp{1})
	std::__throw_domain_error(__N("__heuman_lambda: bad argument"));
      else if (std::abs(std::abs(__k) - _Tp{1}) < _S_eps)
	return __phi / _S_pi_2;
      else if (std::abs(__k) < _S_eps)
	return std::sin(__phi);
      else if (std::real(__phi) < _Real{0})
	return -__heuman_lambda(__k, -__phi);
      else if (std::abs(__phi - _S_pi_2) < _Tp{5} * _S_eps)
	return _Tp{1};
      else if (std::abs(__phi) < _Tp{5} * _S_eps)
	return _Tp{0};
      else if (std::abs(__phi) < _S_pi_2)
	{
	  auto __mc = __k * __k;
	  auto __m = _Tp{1} - __mc;
	  auto __cosphi = std::cos(__phi);
	  auto __sinphi = std::sin(__phi);
	  auto _Delta2 = _Tp{1} - __m * __sinphi * __sinphi;
	  if (std::abs(_Delta2) < _Real{0})
	    _Delta2 = _Tp{0};
	  auto __fact = _Tp{2} * __m * __cosphi * __sinphi
		      / (_S_pi * std::sqrt(_Delta2));
	  auto __arg4 = _Tp{1} - __mc / _Delta2;
	  auto __fact2 = __mc / (_Tp{3} * _Delta2);

	  return __fact * (__ellint_rf(_Tp{0}, __m, _Tp{1})
			+ __fact2 * __ellint_rj(_Tp{0}, __m, _Tp{1}, __arg4));
	}
      else
	{
	  auto __kc = std::sqrt(_Tp{1} - __k * __k);
	  return __ellint_1(__kc, __phi) / __comp_ellint_1(__kc)
	       + __comp_ellint_1(__k) * __jacobi_zeta(__kc, __phi) / _S_pi_2;
	}
    }

_GLIBCXX_END_NAMESPACE_VERSION
} // namespace __detail
} // namespace std

#endif // _GLIBCXX_BITS_SF_ELLINT_TCC

