
// Copyright (C) 2017-2019 Free Software Foundation, Inc.
// Copyright (C) 2020-2022 Edward M. Smith-Rowland
//
// This program is free software; you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation; either version 3 of the License, or (at
// your option) any later version.

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

/** @file ext/complex128.h
 *  This is an internal header file, included by other library headers.
 *  Do not attempt to use it directly. @headername{xxxxx}
 */

#ifndef COMPLEX128_H
#define COMPLEX128_H 1

#include <emsr/float128.h>

#ifdef EMSR_HAVE_FLOAT128

#include <emsr/float128_io.h>

typedef _Complex float __attribute__((mode(TC))) __complex128;

namespace std
{

  // Forward declarations.?
  //template<typename _Tp> class complex;
  //template<> class complex<float>;
  //template<> class complex<double>;
  //template<> class complex<long double>;
  //template<> class complex<__float128>;

  /// complex<__float128> specialization
  template<>
    struct complex<__float128>
    {
      typedef __float128 value_type;

      // From quadmath.h
      //typedef _Complex float __attribute__((mode(TC))) __complex128;
      typedef __complex128 _ComplexT;

      constexpr
      complex(_ComplexT z)
      : _M_value(z)
      { }

      constexpr
      complex(__float128 r = 0.0Q, 
	      __float128 i = 0.0Q)
#if __cplusplus >= 201103L
      : _M_value{ r, i }
      { }
#else
      {
	__real__ _M_value = r;
	__imag__ _M_value = i;
      }
#endif

      constexpr
      complex(const complex<float>& z)
      : _M_value(z.__rep())
      { }

      constexpr
      complex(const complex<double>& z)
      : _M_value(z.__rep())
      { }

      constexpr
      complex(const complex<long double>& z)
      : _M_value(z.__rep())
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
      operator=(__float128 r)
      {
	_M_value = r;
	return *this;
      }

      complex&
      operator+=(__float128 r)
      {
	_M_value += r;
	return *this;
      }

      complex&
      operator-=(__float128 r)
      {
	_M_value -= r;
	return *this;
      }

      complex&
      operator*=(__float128 r)
      {
	_M_value *= r;
	return *this;
      }

      complex&
      operator/=(__float128 r)
      {
	_M_value /= r;
	return *this;
      }

      // The compiler knows how to do this efficiently
      // complex& operator=(const complex&);

      template<typename _Tp>
        complex&
        operator=(const complex<_Tp>& z)
	{
	  __real__ _M_value = z.real();
	  __imag__ _M_value = z.imag();
	  return *this;
	}

      template<typename _Tp>
        complex&
	operator+=(const complex<_Tp>& z)
	{
	  __real__ _M_value += z.real();
	  __imag__ _M_value += z.imag();
	  return *this;
	}

      template<typename _Tp>
        complex&
	operator-=(const complex<_Tp>& z)
	{
	  __real__ _M_value -= z.real();
	  __imag__ _M_value -= z.imag();
	  return *this;
	}

      template<typename _Tp>
        complex&
	operator*=(const complex<_Tp>& z)
	{
	  _ComplexT __t;
	  __real__ __t = z.real();
	  __imag__ __t = z.imag();
	  _M_value *= __t;
	  return *this;
	}

      template<typename _Tp>
        complex&
	operator/=(const complex<_Tp>& z)
	{
	  _ComplexT __t;
	  __real__ __t = z.real();
	  __imag__ __t = z.imag();
	  _M_value /= __t;
	  return *this;
	}

      constexpr _ComplexT
      __rep() const
      { return _M_value; }

    private:

      _ComplexT _M_value;
    };

  inline constexpr
  complex<float>::complex(const complex<__float128>& z)
  : _M_value(z.__rep()) { }

  inline constexpr
  complex<double>::complex(const complex<__float128>& z)
  : _M_value(z.__rep()) { }

  inline constexpr
  complex<long double>::complex(const complex<__float128>& z)
  : _M_value(z.__rep()) { }

  template istream& operator>>(istream&, complex<__float128>&);
  template ostream& operator<<(ostream&, const complex<__float128>&);
  // FIXME!!!
  //template wistream& operator>>(wistream&, complex<__float128>&);
  //template wostream& operator<<(wostream&, const complex<__float128>&);

#if __cplusplus > 201103L

inline namespace literals {
inline namespace complex_literals {

#define __cpp_lib_complex_udls 201309

  inline std::complex<__float128>
  operator""iq(const char* __str)
  { return complex<__float128>(0.0Q, strtoflt128(__str, 0)); }

  constexpr std::complex<__float128>
  operator""iq(unsigned long long __num)
  { return std::complex<__float128>{0.0L, static_cast<__float128>(__num)}; }

} // inline namespace complex_literals
} // inline namespace literals

#endif // C++14

#endif // C++11

} // namespace std

#endif // EMSR_HAVE_FLOAT128

#endif // COMPLEX128_H
