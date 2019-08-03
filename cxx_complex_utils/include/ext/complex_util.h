// Special functions -*- C++ -*-

// Copyright (C) 2016-2019 Free Software Foundation, Inc.
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

/** @file ext/complex_util.h
 * This is an internal header file, included by other library headers.
 * You should not attempt to use it directly.
 */

#ifndef _GLIBCXX_EXT_COMPLEX_UTIL_H
#define _GLIBCXX_EXT_COMPLEX_UTIL_H 1

#pragma GCC system_header

#include <complex>
#include <ratio>
#include <limits>
#include <type_traits>
#include <ext/fp_type_util.h>
#include <ext/math_util.h>

namespace std _GLIBCXX_VISIBILITY(default)
{

  /**
   * We need isnan to be extended to std::complex.
   * Return true if one component of a complex number is NaN.
   */
  template<typename _Tp>
    inline bool
    isnan(const std::complex<_Tp>& __z)
    { return std::isnan(std::real(__z)) || std::isnan(std::imag(__z)); }

  /**
   * Return true if one component of a complex number is inf.
   */
  template<typename _Tp>
    inline bool
    isinf(const std::complex<_Tp>& __z)
    { return isinf(std::real(__z)) || isinf(std::imag(__z)); }

} // namespace std

namespace __gnu_cxx _GLIBCXX_VISIBILITY(default)
{
_GLIBCXX_BEGIN_NAMESPACE_VERSION

//
// The base definitions of __num_traits and fp_promote,
// reside in ext/fp_type_util.h
//

  /**
   * A class to reach into compound numeric types to extract the
   * value or element type - specialized for complex.
   */
  template<typename _Tp>
    struct __num_traits<std::complex<_Tp>>
    {
      using __value_type = typename std::complex<_Tp>::value_type;
    };

  /**
   * Create a complex number NaN.
   */
  template<typename _Tp>
    struct __make_NaN<std::complex<_Tp>>
    {
      constexpr std::complex<_Tp>
      operator()()
      {
	auto __NaN = std::numeric_limits<_Tp>::quiet_NaN();
	return std::complex<_Tp>{__NaN, __NaN};
      }
    };

  /**
   * A function to reliably test a complex number for realness.
   *
   * @param __w  The complex argument
   * @param __mul The multiplier for numeric epsilon for comparison tolerance
   * @return  @c true if @f$ Im(w) @f$ is zero within @f$ mul * epsilon @f$,
   *          @c false otherwize.
   */
  template<typename _Tp>
    bool
    __fp_is_real(const std::complex<_Tp>& __w, const _Tp __mul = _Tp{1})
    { return __gnu_cxx::__fp_is_zero(std::imag(__w), __mul); }

  // Specialize for real numbers.
  template<typename _Tp>
    bool
    __fp_is_real(const _Tp)
    { return true; }

  /**
   * A function to reliably test a complex number for imaginaryness [?].
   *
   * @param __w  The complex argument
   * @param __mul The multiplier for numeric epsilon for comparison tolerance
   * @return  @c true if @f$ Re(w) @f$ is zero within @f$ mul * epsilon @f$,
   *          @c false otherwize.
   */
  template<typename _Tp>
    bool
    __fp_is_imag(const std::complex<_Tp>& __w, const _Tp __mul = _Tp{1})
    { return __gnu_cxx::__fp_is_zero(std::real(__w), __mul); }

  // Specialize for real numbers.
  template<typename _Tp>
    bool
    __fp_is_imag(const _Tp)
    { return false; }

  //  Overloads of integral queries in math_util.h

  /**
   * A function to reliably compare two complex numbers.
   *
   * @param __a The left hand side
   * @param __b The right hand side
   * @param __mul The multiplier for numeric epsilon for comparison
   * @return @c true if a and b are equal to zero
   *         or differ only by @f$ max(a,b) * mul * epsilon @f$
   */
  template<typename _Tp>
    inline bool
    __fp_is_equal(const std::complex<_Tp>& __a, const std::complex<_Tp>& __b,
		  _Tp __mul = _Tp{1})
    {
      const auto _S_eps = std::numeric_limits<_Tp>::epsilon();
      const auto _S_tol = __mul * _S_eps;
      bool __retval = true;
      if (!__fp_is_zero(std::abs(__a), __mul) || !__fp_is_zero(std::abs(__b), __mul))
	// Looks mean, but is necessary that the next line has sense.
	__retval = (std::abs(__a - __b) < __fp_max_abs(__a, __b) * _S_tol);
      return __retval;
    }

  /**
   * A function to reliably compare a complex number and a real number.
   *
   * @param __a The complex left hand side
   * @param __b The real right hand side
   * @param __mul The multiplier for numeric epsilon for comparison
   * @return @c true if a and b are equal to zero
   *         or differ only by @f$ max(a,b) * mul * epsilon @f$
   */
  template<typename _Tp>
    inline bool
    __fp_is_equal(const std::complex<_Tp>& __a, _Tp __b,
		  _Tp __mul = _Tp{1})
    {
      const auto _S_eps = std::numeric_limits<_Tp>::epsilon();
      const auto _S_tol = __mul * _S_eps;
      bool __retval = true;
      if (__fp_is_real(__a, __mul))
	return __fp_is_equal(std::real(__a), __b, __mul);
      else
      	return false;
    }

  /**
   * A function to reliably compare a real number and a complex number.
   *
   * @param __a The real left hand side
   * @param __b The complex right hand side
   * @param __mul The multiplier for numeric epsilon for comparison
   * @return @c true if a and b are equal to zero
   *         or differ only by @f$ max(a,b) * mul * epsilon @f$
   */
  template<typename _Tp>
    inline bool
    __fp_is_equal(const _Tp __a, std::complex<_Tp>& __b,
		  _Tp __mul = _Tp{1})
    {
      const auto _S_eps = std::numeric_limits<_Tp>::epsilon();
      const auto _S_tol = __mul * _S_eps;
      bool __retval = true;
      if (__fp_is_real(__b, __mul))
	return __fp_is_equal(__a, std::real(__b), __mul);
      else
      	return false;
    }

  /**
   * A function to reliably compare a complex number with zero.
   *
   * @param __a The complex point number
   * @param __mul The multiplier for numeric epsilon for comparison
   * @return @c true if a and b are equal to zero
   *         or differ only by @f$ max(a,b) * mul * epsilon @f$
   */
  template<typename _Tp>
    inline bool
    __fp_is_zero(const std::complex<_Tp>& __a, _Tp __mul = _Tp{1})
    { return __fp_is_zero(std::abs(__a), __mul); }

  /**
   * A function to reliably detect if a complex number is a real integer.
   *
   * @param __a The complex number
   * @return @c true if a is an integer within mul * epsilon.
   */
  template<typename _Tp>
    inline __fp_is_integer_t
    __fp_is_integer(const std::complex<_Tp>& __a, _Tp __mul = _Tp{1})
    {
      if (__fp_is_real(__a, __mul))
	return __fp_is_integer(std::real(__a), __mul);
      else
      	return __fp_is_integer_t{false, 0};
    }

  /**
   * A function to reliably detect if a complex number is a real half-integer.
   *
   * @param __a The complex number
   * @return @c true if 2a is an integer within mul * epsilon
   *            and the returned value is half the integer, int(a) / 2.
   */
  template<typename _Tp>
    inline __fp_is_integer_t
    __fp_is_half_integer(const std::complex<_Tp>& __a, _Tp __mul = _Tp{1})
    {
      if (__fp_is_real(__a, __mul))
	return __fp_is_half_integer(std::real(__a), __mul);
      else
      	return __fp_is_integer_t{false, 0};
    }

  /**
   * A function to reliably detect if a complex number is a real
   * half-odd-integer.
   *
   * @param __a The complex number
   * @return @c true if 2a is an odd integer within mul * epsilon
   *            and the returned value is int(a - 1) / 2.
   */
  template<typename _Tp>
    inline __fp_is_integer_t
    __fp_is_half_odd_integer(const std::complex<_Tp>& __a, _Tp __mul = _Tp{1})
    {
      if (__fp_is_real(__a, __mul))
	return __fp_is_half_odd_integer(std::real(__a), __mul);
      else
      	return __fp_is_integer_t{false, 0};
    }

  /**
   * A function to reliably detect if a floating point number is an even integer.
   *
   * @param __a The floating point number
   * @param __mul The multiplier of machine epsilon for the tolerance
   * @return @c true if a is an even integer within mul * epsilon.
   */
  template<typename _Tp>
    inline __fp_is_integer_t
    __fp_is_even_integer(const std::complex<_Tp>& __a, _Tp __mul = _Tp{1})
    {
      if (__fp_is_real(__a, __mul))
	{
	  const auto __integ = __fp_is_integer(std::real(__a), __mul);
	  return __fp_is_integer_t{__integ && !(__integ() & 1), __integ()};
	}
      else
	return __fp_is_integer_t{false, 0};
    }

  /**
   * A function to reliably detect if a floating point number is an odd integer.
   *
   * @param __a The floating point number
   * @param __mul The multiplier of machine epsilon for the tolerance
   * @return @c true if a is an odd integer within mul * epsilon.
   */
  template<typename _Tp>
    inline __fp_is_integer_t
    __fp_is_odd_integer(const std::complex<_Tp>& __a, _Tp __mul = _Tp{1})
    {
      if (__fp_is_real(__a, __mul))
	{
	  const auto __integ = __fp_is_integer(std::real(__a), __mul);
	  return __fp_is_integer_t{__integ && (__integ() & 1), __integ()};
	}
      else
	return __fp_is_integer_t{false, 0};
    }

#if __cplusplus >= 201103L

  /**
   * This is a more modern version of __promote_N in ext/type_traits
   * specialized for complex.
   * This is used for numeric argument promotion of complex and cmath
   */
  template<typename _Tp>
    struct fp_promote_help<std::complex<_Tp>, false>
    {
    private:
      using __vtype = typename std::complex<_Tp>::value_type;
    public:
      using __type = decltype(std::complex<fp_promote_help_t<__vtype>>{});
    };

  /**
   * Type introspection for complex.
   */
  template<typename _Tp>
    struct is_complex
    : public std::false_type
    { };

  template<typename _Tp>
    struct is_complex<const _Tp>
    : public is_complex<_Tp>
    { };

  template<typename _Tp>
    struct is_complex<volatile _Tp>
    : public is_complex<_Tp>
    { };

  template<typename _Tp>
    struct is_complex<const volatile _Tp>
    : public is_complex<_Tp>
    { };

  template<typename _Tp>
    struct is_complex<std::complex<_Tp>>
    : public std::true_type
    { };

  template<typename _Tp>
    using is_complex_t = typename is_complex<_Tp>::type;

  template<typename _Tp>
    constexpr bool is_complex_v = is_complex<_Tp>::value;

#endif // __cplusplus >= 201103L

_GLIBCXX_END_NAMESPACE_VERSION
} // namespace __gnu_cxx

#endif // _GLIBCXX_EXT_COMPLEX_UTIL_H
