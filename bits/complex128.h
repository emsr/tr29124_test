// -*- C++ -*- header.

// Copyright (C) 2016 Free Software Foundation, Inc.
//
// This file is part of the GNU ISO C++ Library.  This library is free
// software; you can redistribute it and/or modify it under the
// terms of the GNU General Public License as published by the
// Free Software Foundation; either version 3, or (at your option)
// any later version.

// This library is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.

// Under Section 7 of GPL version 3, you are granted additional
// permissions described in the GCC Runtime Library Exception, version
// 3.1, as published by the Free Software Foundation.

// You should have received a copy of the GNU General Public License and
// a copy of the GCC Runtime Library Exception along with this program;
// see the files COPYING3 and COPYING.RUNTIME respectively.  If not, see
// <http://www.gnu.org/licenses/>.

/** @file bits/complex128.h
 *  This is an internal header file, included by other library headers.
 *  Do not attempt to use it directly. @headername{xxxxx}
 */

#ifndef _GLIBCXX_BITS_COMPLEX128_H
#define _GLIBCXX_BITS_COMPLEX128_H 1

#pragma GCC system_header

#if !defined(__STRICT_ANSI__) && defined(_GLIBCXX_USE_FLOAT128)

#include <bits/float128.h>
#include <complex>

namespace std _GLIBCXX_VISIBILITY(default)
{
_GLIBCXX_BEGIN_NAMESPACE_VERSION

  /// complex<__float128> specialization
  template<>
    struct complex<__float128>
    {
      typedef __float128 value_type;

      // From quadmath.h
      //typedef _Complex float __attribute__((mode(TC))) __complex128;
      typedef __complex128 _ComplexT;

      _GLIBCXX_CONSTEXPR
      complex(_ComplexT __z)
      : _M_value(__z)
      { }

      _GLIBCXX_CONSTEXPR
      complex(__float128 __r = 0.0Q, 
	      __float128 __i = 0.0Q)
#if __cplusplus >= 201103L
      : _M_value{ __r, __i }
      { }
#else
      {
	__real__ _M_value = __r;
	__imag__ _M_value = __i;
      }
#endif

      _GLIBCXX_CONSTEXPR
      complex(const complex<float>& __z)
      : _M_value(__z.__rep())
      { }

      _GLIBCXX_CONSTEXPR
      complex(const complex<double>& __z)
      : _M_value(__z.__rep())
      { }

      _GLIBCXX_CONSTEXPR
      complex(const complex<long double>& __z)
      : _M_value(__z.__rep())
      { }

#if __cplusplus >= 201103L
      // _GLIBCXX_RESOLVE_LIB_DEFECTS
      // DR 387. std::complex over-encapsulated.
      __attribute ((__abi_tag__ ("cxx11")))
      constexpr __float128 
      real() const
      { return __real__ _M_value; }

      __attribute ((__abi_tag__ ("cxx11")))
      constexpr __float128 
      imag() const
      { return __imag__ _M_value; }
#else
      __float128& 
      real()
      { return __real__ _M_value; }

      const __float128& 
      real() const
      { return __real__ _M_value; }

      __float128& 
      imag()
      { return __imag__ _M_value; }

      const __float128& 
      imag()
      const { return __imag__ _M_value; }
#endif

      // _GLIBCXX_RESOLVE_LIB_DEFECTS
      // DR 387. std::complex over-encapsulated.
      void 
      real(__float128 __val)
      { __real__ _M_value = __val; }

      void 
      imag(__float128 __val)
      { __imag__ _M_value = __val; }

      complex&
      operator=(__float128 __r)
      {
	_M_value = __r;
	return *this;
      }

      complex&
      operator+=(__float128 __r)
      {
	_M_value += __r;
	return *this;
      }

      complex&
      operator-=(__float128 __r)
      {
	_M_value -= __r;
	return *this;
      }

      complex&
      operator*=(__float128 __r)
      {
	_M_value *= __r;
	return *this;
      }

      complex&
      operator/=(__float128 __r)
      {
	_M_value /= __r;
	return *this;
      }

      // The compiler knows how to do this efficiently
      // complex& operator=(const complex&);

      template<typename _Tp>
        complex&
        operator=(const complex<_Tp>& __z)
	{
	  __real__ _M_value = __z.real();
	  __imag__ _M_value = __z.imag();
	  return *this;
	}

      template<typename _Tp>
        complex&
	operator+=(const complex<_Tp>& __z)
	{
	  __real__ _M_value += __z.real();
	  __imag__ _M_value += __z.imag();
	  return *this;
	}

      template<typename _Tp>
        complex&
	operator-=(const complex<_Tp>& __z)
	{
	  __real__ _M_value -= __z.real();
	  __imag__ _M_value -= __z.imag();
	  return *this;
	}

      template<typename _Tp>
        complex&
	operator*=(const complex<_Tp>& __z)
	{
	  _ComplexT __t;
	  __real__ __t = __z.real();
	  __imag__ __t = __z.imag();
	  _M_value *= __t;
	  return *this;
	}

      template<typename _Tp>
        complex&
	operator/=(const complex<_Tp>& __z)
	{
	  _ComplexT __t;
	  __real__ __t = __z.real();
	  __imag__ __t = __z.imag();
	  _M_value /= __t;
	  return *this;
	}

      _GLIBCXX_CONSTEXPR _ComplexT
      __rep() const
      { return _M_value; }

    private:

      _ComplexT _M_value;
    };

  // @todo Ctors from larger types are marked explicit in the smaller classes.
  //inline _GLIBCXX_CONSTEXPR
  //complex<float>::complex(const complex<__float128>& __z)
  //: _M_value(__z.__rep()) { }

  //inline _GLIBCXX_CONSTEXPR
  //complex<double>::complex(const complex<__float128>& __z)
  //: _M_value(__z.__rep()) { }

  //inline _GLIBCXX_CONSTEXPR
  //complex<long double>::complex(const complex<__float128>& __z)
  //: _M_value(__z.__rep()) { }

} // namespace std

  inline __float128
  cabsq(const std::complex<__float128>& z) _GLIBCXX_USE_NOEXCEPT
  { return cabsq(z.__rep()); }

  inline __float128
  cargq(const std::complex<__float128>& z) _GLIBCXX_USE_NOEXCEPT
  { return cargq(z.__rep()); }

  inline __float128
  cimagq(const std::complex<__float128>& z) _GLIBCXX_USE_NOEXCEPT
  { return cimagq(z.__rep()); }

  inline __float128
  crealq(const std::complex<__float128>& z) _GLIBCXX_USE_NOEXCEPT
  { return crealq(z.__rep()); }

  inline __complex128
  cacosq(const std::complex<__float128>& z) _GLIBCXX_USE_NOEXCEPT
  { return cacosq(z.__rep()); }

  inline __complex128
  cacoshq(const std::complex<__float128>& z) _GLIBCXX_USE_NOEXCEPT
  { return cacoshq(z.__rep()); }

  inline __complex128
  casinq(const std::complex<__float128>& z) _GLIBCXX_USE_NOEXCEPT
  { return casinq(z.__rep()); }

  inline __complex128
  casinhq(const std::complex<__float128>& z) _GLIBCXX_USE_NOEXCEPT
  { return casinhq(z.__rep()); }

  inline __complex128
  catanq(const std::complex<__float128>& z) _GLIBCXX_USE_NOEXCEPT
  { return catanq(z.__rep()); }

  inline __complex128
  catanhq(const std::complex<__float128>& z) _GLIBCXX_USE_NOEXCEPT
  { return catanhq(z.__rep()); }

  inline __complex128
  ccosq(const std::complex<__float128>& z) _GLIBCXX_USE_NOEXCEPT
  { return ccosq(z.__rep()); }

  inline __complex128
  ccoshq(const std::complex<__float128>& z) _GLIBCXX_USE_NOEXCEPT
  { return ccoshq(z.__rep()); }

  inline __complex128
  cexpq(const std::complex<__float128>& z) _GLIBCXX_USE_NOEXCEPT
  { return cexpq(z.__rep()); }

  inline __complex128
  cexpiq(__float128 x) _GLIBCXX_USE_NOEXCEPT
  { return cexpiq(x); }

  inline __complex128
  clogq(const std::complex<__float128>& z) _GLIBCXX_USE_NOEXCEPT
  { return clogq(z.__rep()); }

  inline __complex128
  clog10q(const std::complex<__float128>& z) _GLIBCXX_USE_NOEXCEPT
  { return clog10q(z.__rep()); }

  inline __complex128
  conjq(const std::complex<__float128>& z) _GLIBCXX_USE_NOEXCEPT
  { return conjq(z.__rep()); }

  inline __complex128
  cpowq(const std::complex<__float128>& z, const std::complex<__float128>& w)
  _GLIBCXX_USE_NOEXCEPT
  { return cpowq(z.__rep(), w.__rep()); }

  inline __complex128
  cprojq(const std::complex<__float128>& z) _GLIBCXX_USE_NOEXCEPT
  { return cprojq(z.__rep()); }

  inline __complex128
  csinq(const std::complex<__float128>& z) _GLIBCXX_USE_NOEXCEPT
  { return csinq(z.__rep()); }

  inline __complex128
  csinhq(const std::complex<__float128>& z) _GLIBCXX_USE_NOEXCEPT
  { return csinhq(z.__rep()); }

  inline __complex128
  csqrtq(const std::complex<__float128>& z) _GLIBCXX_USE_NOEXCEPT
  { return csqrtq(z.__rep()); }

  inline __complex128
  ctanq(const std::complex<__float128>& z) _GLIBCXX_USE_NOEXCEPT
  { return ctanq(z.__rep()); }

  inline __complex128
  ctanhq(const std::complex<__float128>& z) _GLIBCXX_USE_NOEXCEPT
  { return ctanhq(z.__rep()); };

namespace std
{

  inline __float128
  abs(const __complex128& __z) _GLIBCXX_USE_NOEXCEPT
  { return cabsq(__z); }

  inline __float128
  arg(const __complex128& __z) _GLIBCXX_USE_NOEXCEPT
  { return cargq(__z); }

  inline __float128
  imag(const __complex128& __z) _GLIBCXX_USE_NOEXCEPT
  { return cimagq(__z); }

  inline __float128
  real(const __complex128& __z) _GLIBCXX_USE_NOEXCEPT
  { return crealq(__z); }

  inline __complex128
  acos(const __complex128& __z) _GLIBCXX_USE_NOEXCEPT
  { return cacosq(__z); }

  inline __complex128
  acosh(const __complex128& __z) _GLIBCXX_USE_NOEXCEPT
  { return cacoshq(__z); }

  inline __complex128
  asin(const __complex128& __z) _GLIBCXX_USE_NOEXCEPT
  { return casinq(__z); }

  inline __complex128
  asinh(const __complex128& __z) _GLIBCXX_USE_NOEXCEPT
  { return casinhq(__z); }

  inline __complex128
  atan(const __complex128& __z) _GLIBCXX_USE_NOEXCEPT
  { return catanq(__z); }

  inline __complex128
  atanh(const __complex128& __z) _GLIBCXX_USE_NOEXCEPT
  { return catanhq(__z); }

  inline __complex128
  cos(const __complex128& __z) _GLIBCXX_USE_NOEXCEPT
  { return ccosq(__z); }

  inline __complex128
  cosh(const __complex128& __z) _GLIBCXX_USE_NOEXCEPT
  { return ccoshq(__z); }

  inline __complex128
  exp(const __complex128& __z) _GLIBCXX_USE_NOEXCEPT
  { return cexpq(__z); }

  inline __complex128
  expi(__float128 __x) _GLIBCXX_USE_NOEXCEPT
  { return cexpiq(__x); }

  inline __complex128
  log(const __complex128& __z) _GLIBCXX_USE_NOEXCEPT
  { return clogq(__z); }

  inline __complex128
  log10(const __complex128& __z) _GLIBCXX_USE_NOEXCEPT
  { return clog10q(__z); }

  inline __complex128
  conj(const __complex128& __z) _GLIBCXX_USE_NOEXCEPT
  { return conjq(__z); }

  inline __complex128
  pow(const __complex128& __z, const __complex128& __w) _GLIBCXX_USE_NOEXCEPT
  { return cpowq(__z, __w); }

  inline __complex128
  proj(const __complex128& __z) _GLIBCXX_USE_NOEXCEPT
  { return cprojq(__z); }

  inline __complex128
  sin(const __complex128& __z) _GLIBCXX_USE_NOEXCEPT
  { return csinq(__z); }

  inline __complex128
  sinh(const __complex128& __z) _GLIBCXX_USE_NOEXCEPT
  { return csinhq(__z); }

  inline __complex128
  sqrt(const __complex128& __z) _GLIBCXX_USE_NOEXCEPT
  { return csqrtq(__z); }

  inline __complex128
  tan(const __complex128& __z) _GLIBCXX_USE_NOEXCEPT
  { return ctanq(__z); }

  inline __complex128
  tanh(const __complex128& __z) _GLIBCXX_USE_NOEXCEPT
  { return ctanhq(__z); }

#if _GLIBCXX_USE_C99_COMPLEX

  inline __float128
  __complex_abs(const __complex128& __z)
  { return cabsq(__z); }

  inline __float128
  __complex_arg(const __complex128& __z)
  { return cargq(__z); }

  inline __complex128
  __complex_cos(const __complex128& __z)
  { return ccosq(__z); }

  inline __complex128
  __complex_cosh(const __complex128& __z)
  { return ccoshq(__z); }

  inline __complex128
  __complex_exp(const __complex128& __z)
  { return cexpq(__z); }

  inline __complex128
  __complex_log(const __complex128& __z)
  { return clogq(__z); }

  inline __complex128
  __complex_sin(const __complex128& __z)
  { return csinq(__z); }

  inline __complex128
  __complex_sinh(const __complex128& __z)
  { return csinhq(__z); }      

  inline __complex128
  __complex_sqrt(const __complex128& __z)
  { return csqrtq(__z); }

  inline __complex128
  __complex_tan(const __complex128& __z)
  { return ctanq(__z); }

  inline __complex128
  __complex_tanh(const __complex128& __z)
  { return ctanhq(__z); }

  inline __complex128
  __complex_pow(const __complex128& __x, const __complex128& __y)
  { return cpowq(__x, __y); }

#endif  

_GLIBCXX_END_NAMESPACE_VERSION

#if __cplusplus > 201103L

inline namespace literals {
inline namespace complex_literals {
_GLIBCXX_BEGIN_NAMESPACE_VERSION

  std::complex<__float128>
  operator""iq(const char* __str)
  { return complex<__float128>(0.0Q, strtoflt128(__str, 0)); }

_GLIBCXX_END_NAMESPACE_VERSION
} // inline namespace complex_literals
} // inline namespace literals

#endif // C++14

} // namespace std

#endif // __STRICT_ANSI__ && _GLIBCXX_USE_FLOAT128

#endif // _GLIBCXX_BITS_COMPLEX128_H
