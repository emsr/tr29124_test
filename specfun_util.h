// Special functions -*- C++ -*-

// Copyright (C) 2006-2015 Free Software Foundation, Inc.
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

/** @file bits/specfun_util.h
 *  This is an internal header file, included by other library headers.
 *  Do not attempt to use it directly. @headername{cmath}
 */

//
// ISO C++ 14882 TR29124: Mathematical Special Functions
//

// Written by Edward Smith-Rowland based on numerous mathematics books.

#ifndef _GLIBCXX_BITS_SPECFUN_UTIL_H
#define _GLIBCXX_BITS_SPECFUN_UTIL_H 1

#include <complex>

namespace std _GLIBCXX_VISIBILITY(default)
{
namespace __detail
{
_GLIBCXX_BEGIN_NAMESPACE_VERSION

  /// A class to encapsulate type dependent floating point
  /// constants.  Not everything will be able to be expressed as
  /// type logic.
  template<typename _Tp>
    struct __floating_point_constant
    {
      static const _Tp __value;
    };

  template<typename _Tp>
    struct __num_traits
    {
      using __value_type = _Tp;
    };

  template<>
    template<typename _Tp>
      struct __num_traits<std::complex<_Tp>>
      {
	using __value_type = typename std::complex<_Tp>::value_type;
      };

  template<typename _Tp>
    using __num_traits_t = typename __num_traits<_Tp>::__value_type;

  /// A structure for numeric constants.
  template<typename _Tp>
    struct __numeric_constants
    {
      ///  Constant @f$ \pi @f$.
      static constexpr _Tp __pi() noexcept
      { return static_cast<_Tp>(3.1415926535897932384626433832795029L); }
      ///  Constant @f$ \pi / 2 @f$.
      static constexpr _Tp __pi_2() noexcept
      { return static_cast<_Tp>(1.5707963267948966192313216916397514L); }
      ///  Constant @f$ \pi / 3 @f$.
      static constexpr _Tp __pi_3() noexcept
      { return static_cast<_Tp>(1.0471975511965977461542144610931676L); }
      ///  Constant @f$ \pi / 4 @f$.
      static constexpr _Tp __pi_4() noexcept
      { return static_cast<_Tp>(0.7853981633974483096156608458198757L); }
      /// Constant: @f$ 4 \pi / 3 @f$.
      static constexpr _Tp __4pi_3() noexcept
      { return static_cast<_Tp>(4.1887902047863909846168578443726705L); }
      /// Constant: @f$ 2 \pi @f$.
      static constexpr _Tp __2pi() noexcept
      { return static_cast<_Tp>(6.2831853071795864769252867665590057L); }
      /// Constant: @f$ 4 \pi @f$.
      static constexpr _Tp __4pi() noexcept
      { return static_cast<_Tp>(12.566370614359172953850573533118011L); }
      /// Constant: degrees per radian @f$ 180 / \pi @f$.
      static constexpr _Tp __deg_rad() noexcept
      { return static_cast<_Tp>(57.295779513082320876798154814105170L); }
      /// Constant: radians per degree @f$ \pi / 180 @f$.
      static constexpr _Tp __rad_deg() noexcept
      { return static_cast<_Tp>(0.017453292519943295769236907684886127L); }
      /// Constant: @f$ \sqrt(\pi / 2) @f$.
      static constexpr _Tp __sqrt_pi_2()
      { return static_cast<_Tp>(1.2533141373155002512078826424055226L); }
      ///  Constant @f$ 1 / \pi @f$.
      static constexpr _Tp __1_pi() noexcept
      { return static_cast<_Tp>(0.3183098861837906715377675267450287L); }
      ///  Constant @f$ 2 / \sqrt(\pi) @f$.
      static constexpr _Tp __2_sqrtpi() noexcept
      { return static_cast<_Tp>(1.1283791670955125738961589031215452L); }
      ///  Constant @f$ \sqrt(2) @f$.
      static constexpr _Tp __sqrt2() noexcept
      { return static_cast<_Tp>(1.4142135623730950488016887242096981L); }
      ///  Constant @f$ \sqrt(3) @f$.
      static constexpr _Tp __sqrt3() noexcept
      { return static_cast<_Tp>(1.7320508075688772935274463415058723L); }
      /// Constant: @f$ \sqrt(5) @f$.
      static constexpr _Tp __sqrt_5() noexcept
      { return static_cast<_Tp>(2.2360679774997896964091736687312762L); }
      /// Constant: @f$ \sqrt(7) @f$.
      static constexpr _Tp __sqrt_7() noexcept
      { return static_cast<_Tp>(2.6457513110645905905016157536392604L); }
      ///  Constant @f$ \sqrt(\pi/2) @f$.
      static constexpr _Tp __sqrtpio2() noexcept
      { return static_cast<_Tp>(1.2533141373155002512078826424055226L); }
      ///  Constant @f$ 1 / sqrt(2) @f$.
      static constexpr _Tp __sqrt1_2() noexcept
      { return static_cast<_Tp>(0.7071067811865475244008443621048490L); }
      ///  Constant @f$ \log(\pi) @f$.
      static constexpr _Tp __lnpi() noexcept
      { return static_cast<_Tp>(1.1447298858494001741434273513530587L); }
      ///  Constant Euler's constant @f$ \gamma_E @f$.
      static constexpr _Tp __gamma_e() noexcept
      { return static_cast<_Tp>(0.5772156649015328606065120900824024L); }
      /// Constant: Golden Ratio @f$ \phi = (1 + \sqrt{5})/2 @f$.
      static constexpr _Tp __phi() noexcept
      { return static_cast<_Tp>(1.6180339887498948482045868343656381L); }
      ///  Constant Euler-Mascheroni @f$ e @f$
      static constexpr _Tp __euler() noexcept
      { return static_cast<_Tp>(2.7182818284590452353602874713526625L); }
      /// Constant: @f$ \pi^2/6 @f$.
      static constexpr _Tp __pi2_6()
      { return static_cast<_Tp>(1.6449340668482264364724151666460251L); }
      /// Constant: @f$ min(RealType) @f$.
      static constexpr _Tp __min() noexcept
      { return std::numeric_limits<_Tp>::min(); }
      /// Constant: @f$ max(RealType) @f$.
      static constexpr _Tp __max() noexcept
      { return std::numeric_limits<_Tp>::max(); }
      /// Constant: @f$ inf(RealType) @f$.
      static constexpr _Tp __inf() noexcept
      { return std::numeric_limits<_Tp>::infinity(); }
      /// Constant: @f$ NaN(RealType) @f$.
      static constexpr _Tp __NaN() noexcept
      { return std::numeric_limits<_Tp>::quiet_NaN(); }
    };

#if _GLIBCXX_USE_C99_MATH && !_GLIBCXX_USE_C99_FP_MACROS_DYNAMIC

  /// This is a wrapper for the isnan function. Otherwise, for NaN,
  /// all comparisons result in false. If/when we build a std::isnan
  /// out of intrinsics, this will disappear completely in favor of
  /// std::isnan.
  template<typename _Tp>
    inline bool
    __isnan(_Tp __x)
    { return std::isnan(__x); }

#else

  template<typename _Tp>
    inline bool
    __isnan(_Tp __x)
    { return __builtin_isnan(__x); }

  template<>
    inline bool
    __isnan<float>(float __x)
    { return __builtin_isnanf(__x); }

  template<>
    inline bool
    __isnan<long double>(long double __x)
    { return __builtin_isnanl(__x); }

#endif

  template<typename _Tp>
    inline bool
    __isnan(const std::complex<_Tp>& __x)
    { return __isnan(std::real(__x)) || __isnan(std::imag(__x)); }


  template<typename _Tp>
    inline _Tp
    __quiet_NaN(_Tp)
    { return std::numeric_limits<_Tp>::quiet_NaN(); }

  template<typename _Tp>
    inline std::complex<_Tp>
    __quiet_NaN(std::complex<_Tp>)
    {
      auto __nan = std::numeric_limits<_Tp>::quiet_NaN();
      return std::complex<_Tp>(__nan, __nan);
    }


  template<typename _Tp>
    inline _Tp
    __infinity(_Tp)
    { return std::numeric_limits<_Tp>::infinity(); }

  template<typename _Tp>
    inline std::complex<_Tp>
    __infinity(std::complex<_Tp>)
    {
      auto __inf = std::numeric_limits<_Tp>::infinity();
      return std::complex<_Tp>(__inf, __inf);
    }

  template<typename _Tp>
    class _KahanSum
    {
    public:

      _KahanSum() = default;

      explicit _KahanSum(_Tp __sum)
      : _M_sum{__sum}, _M_term{}, _M_temp{}
      { }

      _KahanSum&
      operator+=(_Tp __term)
      {
	this->_M_term = __term - this->_M_temp;
	this->_M_temp = this->_M_sum;
	this->_M_sum += this->_M_term;
	this->_M_temp = this->_M_term - (this->_M_sum - this->_M_temp);
	return *this;
      }

      _KahanSum&
      operator=(_Tp __sum)
      {
	this->_M_term = _Tp{0};
	this->_M_temp = _Tp{0};
	this->_M_sum = __sum;
	return *this;
      }

      _KahanSum&
      operator-=(_Tp __term)
      { return operator+=(-__term); }

      _Tp
      operator()() const
      { return this->_M_sum; }

    private:

      _Tp _M_sum;
      _Tp _M_term;
      _Tp _M_temp;
    };

_GLIBCXX_END_NAMESPACE_VERSION
} // namespace __detail
}

#endif // _GLIBCXX_BITS_SPECFUN_UTIL_H

