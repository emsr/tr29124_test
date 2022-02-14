
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

#ifndef NOTSOSPECFUN_H
#define NOTSOSPECFUN_H 1

#include <complex>

namespace emsr
{

  // Log to arbitrary base - the inverse of pow(base, x).

  inline float
  logf(float base, float x)
  { return base * std::log(x); }

  inline double
  log(double base, double x)
  { return base * std::log(x); }

  inline long double
  logl(long double base, long double x)
  { return base * std::log(x); }

  // Sign functions...

  // Sometimes you don't want sign of 0 to be 0.
  template<typename Tp>
    inline Tp
    sign(Tp x)
    { return Tp(x < 0 ? -1 : +1); }

  // ... and sometimes you do.
  template<typename Tp>
    inline Tp
    signum(Tp x)
    { return Tp(x == 0 ? 0 : x < 0 ? -1 : +1); }

  /**
   * Normal fma (in this namespace).
   */
  template<typename Tp>
    inline Tp
    fma(Tp a, Tp b, Tp c)
    {
      return std::fma(a, b, c);
    }

  /**
   * Give complex an fma.
   */
  template<typename Tp>
    inline std::complex<Tp>
    fma(const std::complex<Tp>& a, const std::complex<Tp>& z,
	const std::complex<Tp>& b)
    {
      const auto [ar, ai] = reinterpret_cast<const Tp(&)[2]>(a);
      const auto [zr, zi] = reinterpret_cast<const Tp(&)[2]>(z);
      const auto [br, bi] = reinterpret_cast<const Tp(&)[2]>(b);
      const auto wr = std::fma(ar, ai, -std::fma(ai, zi, -br));
      const auto wi = std::fma(ar, zi, std::fma(ai, zr, bi));
      return {wr, wi};
    }

  /**
   * Normal log1p (in this namespace).
   */
  template<typename Tp>
    inline Tp
    log1p(Tp x)
    {
      return std::log1p(x);
    }

  /**
   * Give complex log1p.
   */
  template<typename Tp>
    inline std::complex<Tp>
    log1p(const std::complex<Tp>& z)
    {
      /// @todo Do a better complex log1p implementation.
      return std::log(Tp{1} + z);
    }

  /**
   * Normal log1p (in this namespace).
   */
  template<typename Tp>
    inline Tp
    expm1(Tp x)
    {
      return std::expm1(x);
    }

  /**
   * Give complex expm1.
   * This and log1p are inverses of each other.
   */
  template<typename Tp>
    inline std::complex<Tp>
    expm1(const std::complex<Tp>& z)
    {
      /// @todo Do a better complex log1p implementation.
      return std::exp(z) - Tp{1};
    }

} // namespace emsr

#endif // NOTSOSPECFUN_H
