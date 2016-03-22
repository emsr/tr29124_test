// Special functions -*- C++ -*-

// Copyright (C) 2006-2016 Free Software Foundation, Inc.
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
 * This is an internal header file, included by other library headers.
 * Do not attempt to use it directly. @headername{cmath}
 */

//
// ISO C++ 14882 TR29124: Mathematical Special Functions
//

#ifndef _GLIBCXX_BITS_SPECFUN_UTIL_H
#define _GLIBCXX_BITS_SPECFUN_UTIL_H 1

#if __cplusplus >= 201103L
#  include <ratio>
#  include <complex>
#endif

#if !defined(__STRICT_ANSI__) && defined(_GLIBCXX_USE_FLOAT128)
#  include <quadmath.h>
#endif // __STRICT_ANSI__ && _GLIBCXX_USE_FLOAT128

namespace std _GLIBCXX_VISIBILITY(default)
{
namespace __detail
{
_GLIBCXX_BEGIN_NAMESPACE_VERSION

  /**
   * A class to reach into compound numeric types to extract the
   * value or element type.  This will be specialized for complex
   * and other types as appropriate.
   */
  template<typename _Tp>
    struct __num_traits
    {
      using __value_type = _Tp;
    };

  template<typename _Tp>
    using __num_traits_t = typename __num_traits<_Tp>::__value_type;


#if __cplusplus >= 201103L
  /**
   * Return a fraction as a real number.
   */
  template<intmax_t _Num, intmax_t _Den = 1, typename _Tp = double>
    inline constexpr _Tp
    __frac()
    {
      using __rat_t = std::ratio<_Num, _Den>;
      return _Tp(__rat_t::num) / _Tp(__rat_t::den);
    }
#endif


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

#if !defined(__STRICT_ANSI__) && defined(_GLIBCXX_USE_FLOAT128)
  template<>
    inline bool
    __isnan<__float128>(__float128 __x)
    { return __builtin_isnanq(__x); }
#endif // __STRICT_ANSI__ && _GLIBCXX_USE_FLOAT128

#endif // _GLIBCXX_USE_C99_MATH && !_GLIBCXX_USE_C99_FP_MACROS_DYNAMIC

_GLIBCXX_END_NAMESPACE_VERSION
} // namespace __detail
} // namespace std

namespace std _GLIBCXX_VISIBILITY(default)
{
namespace __detail
{
_GLIBCXX_BEGIN_NAMESPACE_VERSION

  /**
   * A Sum is default constructable
   * It has operator+=(Tp)
   * It has operator-=(Tp)
   * It has operator() -> _Tp
   * It has operator bool()
   * It has num_terms() -> std::size_t
   * 
   */

  /**
   * This is a basic naive sum.
   */
  template<typename _Tp>
    class _BasicSum
    {
    public:

      using value_type = _Tp;

      ///  Default constructor.
      _BasicSum() = default;

      ///  Constructor taking the first term.
      explicit
      _BasicSum(value_type __first_term)
      : _BasicSum{}
      {
	operator+=(__first_term);
      }

      /// Add a new term to the sum.
      _BasicSum&
      operator+=(value_type __term)
      {
	if (!this->_M_converged)
	  {
	    if (__isnan(__term))
	      throw std::runtime_error("_BasicSum: bad term");
	    if (std::abs(__term) == std::numeric_limits<value_type>::infinity())
	      throw std::runtime_error("_BasicSum: infinite term");
	    ++this->_M_num_terms;
	    this->_M_term = __term;
	    this->_M_sum += this->_M_term;
	  }
	return *this;
      }

      /// Subtract a new term from the sum.
      _BasicSum&
      operator-=(value_type __term)
      { return operator+=(-__term); }

      /// Return true if the sum converged.
      operator
      bool() const
      { return !this->_M_converged; }

      /// Return the current value of the sum.
      value_type
      operator()() const
      { return this->_M_sum; }

      /// Return the current number of terms contributing to the sum.
      std::size_t
      num_terms() const
      { return this->_M_num_terms; }

      ///  Reset the sum to it's initial state.
      _BasicSum&
      reset()
      {
	this->_M_sum = value_type{0};
	this->_M_term = value_type{0};
	this->_M_num_terms = 0;
	this->_M_converged = false;
	return *this;
      }

      ///  Restart the sum with the first new term.
      _BasicSum&
      reset(value_type __first_term)
      {
	this->reset();
	this->operator+=(__first_term);
	return *this;
      }

    private:

      value_type _M_sum;
      value_type _M_term;
      std::size_t _M_num_terms;
      bool _M_converged;
    };

  /**
   * This is a Kahan sum which tries to account for roundoff error.
   */
  template<typename _Tp>
    class _KahanSum
    {
    public:

      using value_type = _Tp;

      ///  Default constructor.
      _KahanSum() = default;

      ///  Constructor taking the first term.
      explicit
      _KahanSum(value_type __first_term)
      : _KahanSum{}
      {
	operator+=(__first_term);
      }

      /// Add a new term to the sum.
      _KahanSum&
      operator+=(value_type __term)
      {
	if (!this->_M_converged)
	  {
	    if (__isnan(__term))
	      throw std::runtime_error("_KahanSum: bad term");
	    if (std::abs(__term) == std::numeric_limits<value_type>::infinity())
	      throw std::runtime_error("_KahanSum: infinite term");
	    ++this->_M_num_terms;
	    this->_M_term = __term - this->_M_temp;
	    this->_M_temp = this->_M_sum;
	    this->_M_sum += this->_M_term;
	    this->_M_temp = this->_M_term - (this->_M_sum - this->_M_temp);
	  }
	return *this;
      }

      /// Subtract a new term from the sum.
      _KahanSum&
      operator-=(value_type __term)
      { return operator+=(-__term); }

      /// Return true if the sum converged.
      operator
      bool() const
      { return !this->_M_converged; }

      /// Return the current value of the sum.
      value_type
      operator()() const
      { return this->_M_sum; }

      /// Return the current number of terms contributing to the sum.
      std::size_t
      num_terms() const
      { return this->_M_num_terms; }

      ///  Reset the sum to it's initial state.
      _KahanSum&
      reset()
      {
	this->_M_sum = value_type{0};
	this->_M_term = value_type{0};
	this->_M_temp = value_type{0};
	this->_M_num_terms = 0;
	this->_M_converged = false;
	return *this;
      }

      ///  Restart the sum with the first new term.
      _KahanSum&
      reset(value_type __first_term)
      {
	this->reset();
	this->operator+=(__first_term);
	return *this;
      }

    private:

      value_type _M_sum;
      value_type _M_term;
      value_type _M_temp;
      std::size_t _M_num_terms;
      bool _M_converged;
    };

_GLIBCXX_END_NAMESPACE_VERSION
} // namespace __detail
} // namespace std

namespace __gnu_cxx _GLIBCXX_VISIBILITY(default)
{

#if __cplusplus >= 201103L

  /**
   * This is a more modern version of __promote_N in ext/type_traits.
   * This is used for numeric argument promotion of complex and cmath.
   */
  template<typename _Tp, bool = std::is_integral<_Tp>::value>
    struct __promote_help
    { using __type = double; };

  // No nested __type member for non-integer non-floating point types,
  // allows this type to be used for SFINAE to constrain overloads in
  // <cmath> and <complex> to only the intended types.
  template<typename _Tp>
    struct __promote_help<_Tp, false>
    { };

  template<>
    struct __promote_help<float>
    { using __type = float; };

  template<>
    struct __promote_help<double>
    { using __type = double; };

  template<>
    struct __promote_help<long double>
    { using __type = long double; };

#if !defined(__STRICT_ANSI__) && defined(_GLIBCXX_USE_FLOAT128)
  template<>
    struct __promote_help<__float128>
    { using __type = __float128; };
#endif

  template<typename... _Tps>
    using __promote_help_t = typename __promote_help<_Tps...>::__type;

  // Decay refs and cv...
  // Alternatively we could decay refs and propagate cv to promoted type.
  template<typename _Tp, typename... _Tps>
    struct __promote_num
    { using __type = decltype(__promote_help_t<std::decay_t<_Tp>>{}
		   + typename __promote_num<_Tps...>::__type{}); };

  template<>
    template<typename _Tp>
      struct __promote_num<_Tp>
      { using __type = decltype(__promote_help_t<std::decay_t<_Tp>>{}); };

  template<typename... _Tps>
    using __promote_num_t = typename __promote_num<_Tps...>::__type;

#endif // __cplusplus >= 201103L

} // namespace __gnu_cxx

#endif // _GLIBCXX_BITS_SPECFUN_UTIL_H

