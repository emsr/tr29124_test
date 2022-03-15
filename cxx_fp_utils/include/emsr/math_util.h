
// Copyright (C) 2016-2019 Free Software Foundation, Inc.
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

#ifndef MATH_UTIL_H
#define MATH_UTIL_H 1

#include <limits>
#include <cmath>
#include <algorithm> // min, max

namespace emsr
{

  /**
   * Return -1 if the integer argument is odd and +1 if it is even.
   */
  template<typename Tp, typename _IntTp>
    inline constexpr Tp
    parity(_IntTp k) noexcept
    { return (k & 1) ? Tp{-1} : Tp{+1}; }

  /**
   * A function to return the maximum of the absolute values of two numbers
   * ... so we won't include everything.
   * @param a The left hand side
   * @param b The right hand side
   */
  template<typename Tp>
    inline constexpr Tp
    fp_max_abs(Tp a, Tp b) noexcept
    {
      if (std::isnan(a) || std::isnan(b))
	return std::numeric_limits<Tp>::quiet_NaN();
      else
	{
	  const auto aa = std::abs(a);
	  const auto bb = std::abs(b);
	  return std::max(aa, bb);
	}
    }

  /**
   * A function to reliably compare two floating point numbers.
   *
   * @param a The left hand side
   * @param b The right hand side
   * @param mul The multiplier for numeric epsilon for comparison
   * @return @c true if a and b are equal to zero
   *         or differ only by @f$ max(a,b) * mul * epsilon @f$
   */
  template<typename Tp>
    inline constexpr bool
    fp_is_equal(Tp a, Tp b, Tp mul = Tp{1}) noexcept
    {
      if (std::isnan(a) || std::isnan(b) || std::isnan(mul))
	return false;
      else
	{
	  const auto _S_tol = mul * std::numeric_limits<Tp>::epsilon();
	  bool retval = true;
	  if ((a != Tp{0}) || (b != Tp{0}))
	    // Looks mean, but is necessary that the next line has sense.
	    retval = (std::abs(a - b) < fp_max_abs(a, b) * _S_tol);
	  return retval;
	}
    }

  /**
   * A function to reliably compare a floating point number with zero.
   *
   * @param a The floating point number
   * @param mul The multiplier for numeric epsilon for comparison
   * @return @c true if a and b are equal to zero
   *         or differ only by @f$ max(a,b) * mul * epsilon @f$
   */
  template<typename Tp>
    inline constexpr bool
    fp_is_zero(Tp a, Tp mul = Tp{1}) noexcept
    {
      if (std::isnan(a) || std::isnan(mul))
	return false;
      else
	{
	  const auto _S_tol = mul * std::numeric_limits<Tp>::epsilon();
	  if (a != Tp{0})
	    return (std::abs(a) < _S_tol);
	  else
            return true;
	}
    }

  /**
   * A struct returned by floating point integer queries.
   */
  struct fp_is_integer_t
  {
    // A flag indicating whether the floating point number is integralish.
    bool fp_is_integral = false;

    // An integer related to the floating point integral value.
    int value = 0;

    // Return is_integral in a boolean context.
    constexpr operator bool() const noexcept
    { return this->fp_is_integral; }

    // Return value.
    constexpr int
    operator()() const noexcept
    { return this->value; }
  };

  /**
   * A function to reliably detect if a floating point number is an integer.
   *
   * @param a The floating point number
   * @param mul The multiplier of machine epsilon for the tolerance
   * @return @c true if a is an integer within mul * epsilon.
   */
  template<typename Tp>
    inline constexpr fp_is_integer_t
    fp_is_integer(Tp a, Tp mul = Tp{1}) noexcept
    {
      if (std::isnan(a) || std::isnan(mul))
	return fp_is_integer_t{false, 0};
      else
	{
	  const auto n = static_cast<int>(std::nearbyint(a));
	  const auto eq = fp_is_equal(a, Tp(n), mul);
	  return fp_is_integer_t{eq, n};
	}
    }

  /**
   * A function to reliably detect if a floating point number is a half-integer.
   *
   * @param a The floating point number
   * @param mul The multiplier of machine epsilon for the tolerance
   * @return @c true if 2a is an integer within mul * epsilon
   *            and the returned value is half the integer, int(a) / 2. 
   */
  template<typename Tp>
    inline constexpr fp_is_integer_t
    fp_is_half_integer(Tp a, Tp mul = Tp{1}) noexcept
    {
      if (std::isnan(a) || std::isnan(mul))
	return fp_is_integer_t{false, 0};
      else
	{
	  const auto n = static_cast<int>(std::nearbyint(Tp{2} * a));
	  const auto eq = fp_is_equal(Tp{2} * a, Tp(n), mul);
	  return fp_is_integer_t{eq, n / 2};
	}
    }

  /**
   * A function to reliably detect if a floating point number is a
   * half-odd-integer.
   *
   * @param a The floating point number
   * @param mul The multiplier of machine epsilon for the tolerance
   * @return @c true if 2a is an odd integer within mul * epsilon
   *            and the returned value is int(a - 1) / 2.
   */
  template<typename Tp>
    inline constexpr fp_is_integer_t
    fp_is_half_odd_integer(Tp a, Tp mul = Tp{1}) noexcept
    {
      if (std::isnan(a) || std::isnan(mul))
	return fp_is_integer_t{false, 0};
      else
	{
	  const auto n = static_cast<int>(std::nearbyint(Tp{2} * a));
	  const bool halfodd = ((n & 1) == 1)
				&& fp_is_equal(Tp{2} * a, Tp(n), mul);
	  return fp_is_integer_t{halfodd, (n - 1) / 2};
	}
    }

  /**
   * A function to reliably detect if a floating point number is an even integer.
   *
   * @param a The floating point number
   * @param mul The multiplier of machine epsilon for the tolerance
   * @return @c true if a is an even integer within mul * epsilon.
   */
  template<typename Tp>
    inline constexpr fp_is_integer_t
    fp_is_even_integer(Tp a, Tp mul = Tp{1}) noexcept
    {
      if (std::isnan(a) || std::isnan(mul))
	return fp_is_integer_t{false, 0};
      else
	{
	  const auto integ = fp_is_integer(a, mul);
	  return fp_is_integer_t{integ && !(integ() & 1), integ()};
	}
    }

  /**
   * A function to reliably detect if a floating point number is an odd integer.
   *
   * @param a The floating point number
   * @param mul The multiplier of machine epsilon for the tolerance
   * @return @c true if a is an odd integer within mul * epsilon.
   */
  template<typename Tp>
    inline constexpr fp_is_integer_t
    fp_is_odd_integer(Tp a, Tp mul = Tp{1}) noexcept
    {
      if (std::isnan(a) || std::isnan(mul))
	return fp_is_integer_t{false, 0};
      else
	{
	  const auto integ = fp_is_integer(a, mul);
	  return fp_is_integer_t{integ && (integ() & 1), integ()};
	}
    }

} // namespace emsr

#endif // MATH_UTIL_H
