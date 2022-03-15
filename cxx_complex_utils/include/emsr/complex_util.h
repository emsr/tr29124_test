
// Copyright (C) 2016-2019 Free Software Foundation, Inc.
// Copyright (C) 2020-2022 Edward M. Smith-Rowland
//
// This program is free software; you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation; either version 3 of the License, or (at
// your option) any later version.
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

#ifndef COMPLEX_UTIL_H
#define COMPLEX_UTIL_H 1

#include <ratio>
#include <limits>
#include <type_traits>
#include <complex>

#include <emsr/fp_type_util.h>
#include <emsr/math_util.h>

namespace std
{

  /**
   * We need isnan to be extended to std::complex.
   * Return true if one component of a complex number is NaN.
   */
  template<typename Tp>
    inline bool
    isnan(const std::complex<Tp>& z)
    { return std::isnan(std::real(z)) || std::isnan(std::imag(z)); }

  /**
   * Return true if one component of a complex number is inf.
   */
  template<typename Tp>
    inline bool
    isinf(const std::complex<Tp>& z)
    { return isinf(std::real(z)) || isinf(std::imag(z)); }

} // namespace std

namespace emsr
{

//
// The base definitions of num_traits and fp_promote,
// reside in ext/fp_type_util.h
//

  /**
   * A class to reach into compound numeric types to extract the
   * value or element type - specialized for complex.
   */
  template<typename Tp>
    struct num_traits<std::complex<Tp>>
    {
      using value_type = typename std::complex<Tp>::value_type;
    };

  /**
   * Create a complex number NaN.
   */
  template<typename Tp>
    struct make_NaN<std::complex<Tp>>
    {
      constexpr std::complex<Tp>
      operator()()
      {
	auto NaN = std::numeric_limits<Tp>::quiet_NaN();
	return std::complex<Tp>{NaN, NaN};
      }
    };

  /**
   * A function to reliably test a complex number for realness.
   *
   * @param w  The complex argument
   * @param mul The multiplier for numeric epsilon for comparison tolerance
   * @return  @c true if @f$ Im(w) @f$ is zero within @f$ mul * epsilon @f$,
   *          @c false otherwize.
   */
  template<typename Tp>
    bool
    fp_is_real(const std::complex<Tp>& w, const Tp mul = Tp{1})
    { return fp_is_zero(std::imag(w), mul); }

  // Specialize for real numbers.
  template<typename Tp>
    bool
    fp_is_real(const Tp)
    { return true; }

  /**
   * A function to reliably test a complex number for imaginaryness [?].
   *
   * @param w  The complex argument
   * @param mul The multiplier for numeric epsilon for comparison tolerance
   * @return  @c true if @f$ Re(w) @f$ is zero within @f$ mul * epsilon @f$,
   *          @c false otherwize.
   */
  template<typename Tp>
    bool
    fp_is_imag(const std::complex<Tp>& w, const Tp mul = Tp{1})
    { return fp_is_zero(std::real(w), mul); }

  // Specialize for real numbers.
  template<typename Tp>
    bool
    fp_is_imag(const Tp)
    { return false; }

  //  Overloads of integral queries in math_util.h

  /**
   * A function to reliably compare two complex numbers.
   *
   * @param a The left hand side
   * @param b The right hand side
   * @param mul The multiplier for numeric epsilon for comparison
   * @return @c true if a and b are equal to zero
   *         or differ only by @f$ max(a,b) * mul * epsilon @f$
   */
  template<typename Tp>
    inline bool
    fp_is_equal(const std::complex<Tp>& a, const std::complex<Tp>& b,
		  Tp mul = Tp{1})
    {
      const auto s_eps = std::numeric_limits<Tp>::epsilon();
      const auto s_tol = mul * s_eps;
      bool retval = true;
      if (!fp_is_zero(std::abs(a), mul) || !fp_is_zero(std::abs(b), mul))
	// Looks mean, but is necessary that the next line has sense.
	retval = (std::abs(a - b) < fp_max_abs(a, b) * s_tol);
      return retval;
    }

  /**
   * A function to reliably compare a complex number and a real number.
   *
   * @param a The complex left hand side
   * @param b The real right hand side
   * @param mul The multiplier for numeric epsilon for comparison
   * @return @c true if a and b are equal to zero
   *         or differ only by @f$ max(a,b) * mul * epsilon @f$
   */
  template<typename Tp>
    inline bool
    fp_is_equal(const std::complex<Tp>& a, Tp b,
		  Tp mul = Tp{1})
    {
      if (fp_is_real(a, mul))
	return fp_is_equal(std::real(a), b, mul);
      else
      	return false;
    }

  /**
   * A function to reliably compare a real number and a complex number.
   *
   * @param a The real left hand side
   * @param b The complex right hand side
   * @param mul The multiplier for numeric epsilon for comparison
   * @return @c true if a and b are equal to zero
   *         or differ only by @f$ max(a,b) * mul * epsilon @f$
   */
  template<typename Tp>
    inline bool
    fp_is_equal(const Tp a, std::complex<Tp>& b,
		  Tp mul = Tp{1})
    {
      if (fp_is_real(b, mul))
	return fp_is_equal(a, std::real(b), mul);
      else
      	return false;
    }

  /**
   * A function to reliably compare a complex number with zero.
   *
   * @param a The complex point number
   * @param mul The multiplier for numeric epsilon for comparison
   * @return @c true if a and b are equal to zero
   *         or differ only by @f$ max(a,b) * mul * epsilon @f$
   */
  template<typename Tp>
    inline bool
    fp_is_zero(const std::complex<Tp>& a, Tp mul = Tp{1})
    { return fp_is_zero(std::abs(a), mul); }

  /**
   * A function to reliably detect if a complex number is a real integer.
   *
   * @param a The complex number
   * @return @c true if a is an integer within mul * epsilon.
   */
  template<typename Tp>
    inline fp_is_integer_t
    fp_is_integer(const std::complex<Tp>& a, Tp mul = Tp{1})
    {
      if (fp_is_real(a, mul))
	return fp_is_integer(std::real(a), mul);
      else
      	return fp_is_integer_t{false, 0};
    }

  /**
   * A function to reliably detect if a complex number is a real half-integer.
   *
   * @param a The complex number
   * @return @c true if 2a is an integer within mul * epsilon
   *            and the returned value is half the integer, int(a) / 2.
   */
  template<typename Tp>
    inline fp_is_integer_t
    fp_is_half_integer(const std::complex<Tp>& a, Tp mul = Tp{1})
    {
      if (fp_is_real(a, mul))
	return fp_is_half_integer(std::real(a), mul);
      else
      	return fp_is_integer_t{false, 0};
    }

  /**
   * A function to reliably detect if a complex number is a real
   * half-odd-integer.
   *
   * @param a The complex number
   * @return @c true if 2a is an odd integer within mul * epsilon
   *            and the returned value is int(a - 1) / 2.
   */
  template<typename Tp>
    inline fp_is_integer_t
    fp_is_half_odd_integer(const std::complex<Tp>& a, Tp mul = Tp{1})
    {
      if (fp_is_real(a, mul))
	return fp_is_half_odd_integer(std::real(a), mul);
      else
      	return fp_is_integer_t{false, 0};
    }

  /**
   * A function to reliably detect if a floating point number is an even integer.
   *
   * @param a The floating point number
   * @param mul The multiplier of machine epsilon for the tolerance
   * @return @c true if a is an even integer within mul * epsilon.
   */
  template<typename Tp>
    inline fp_is_integer_t
    fp_is_even_integer(const std::complex<Tp>& a, Tp mul = Tp{1})
    {
      if (fp_is_real(a, mul))
	{
	  const auto integ = fp_is_integer(std::real(a), mul);
	  return fp_is_integer_t{integ && !(integ() & 1), integ()};
	}
      else
	return fp_is_integer_t{false, 0};
    }

  /**
   * A function to reliably detect if a floating point number is an odd integer.
   *
   * @param a The floating point number
   * @param mul The multiplier of machine epsilon for the tolerance
   * @return @c true if a is an odd integer within mul * epsilon.
   */
  template<typename Tp>
    inline fp_is_integer_t
    fp_is_odd_integer(const std::complex<Tp>& a, Tp mul = Tp{1})
    {
      if (fp_is_real(a, mul))
	{
	  const auto integ = fp_is_integer(std::real(a), mul);
	  return fp_is_integer_t{integ && (integ() & 1), integ()};
	}
      else
	return fp_is_integer_t{false, 0};
    }

  /**
   * This is a more modern version of promote_N in ext/type_traits
   * specialized for complex.
   * This is used for numeric argument promotion of complex and cmath
   */
  template<typename Tp>
    struct fp_promote_help<std::complex<Tp>, false>
    {
    private:
      using vtype = typename std::complex<Tp>::value_type;
    public:
      using type = decltype(std::complex<fp_promote_help_t<vtype>>{});
    };

  /**
   * Type introspection for complex.
   */
  template<typename Tp>
    struct is_complex
    : public std::false_type
    { };

  template<typename Tp>
    struct is_complex<const Tp>
    : public is_complex<Tp>
    { };

  template<typename Tp>
    struct is_complex<volatile Tp>
    : public is_complex<Tp>
    { };

  template<typename Tp>
    struct is_complex<const volatile Tp>
    : public is_complex<Tp>
    { };

  template<typename Tp>
    struct is_complex<std::complex<Tp>>
    : public std::true_type
    { };

  template<typename Tp>
    using is_complex_t = typename is_complex<Tp>::type;

  template<typename Tp>
    constexpr bool is_complex_v = is_complex<Tp>::value;

} // namespace emsr

#endif // COMPLEX_UTIL_H
