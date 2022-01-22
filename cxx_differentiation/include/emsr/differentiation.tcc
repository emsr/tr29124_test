// Numerical Differentiation for -*- C++ -*-

// Copyright (C) 2019 Free Software Foundation, Inc.
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

#ifndef DIFFERENTIATION_TCC
#define DIFFERENTIATION_TCC 1

/**
 * @file  differentiation.tcc
 *
 * This file contains the outline implementation of the differentiation
 * routines.
 */

#include <cmath>
#include <limits>
#include <array>
#include <complex>

  /**
   * Compute the derivative using the 5-point rule
   * (x-h, x-h/2, x, x+h/2, x+h).
   * Note that the central point is not used.
   * @f[
   *   f'_3(x) = \frac{f(x + h) - f(x - h)}{2h}
   * @f]
   * @f[
   *   f'_5(x) = \frac{4}{3}\frac{f(x + h/2) - f(x - h/2)}{2h}
   *           - \frac{1}{3}f'_3(x)
   * @f]
   *
   * Compute the error using the difference between the 5-point
   * and the 3-point rule (x-h, x, x+h).
   * Again the central point is not used.
   *
   * @tparam Func A function type callable with a numeric type
   *          and returning the same.
   * @tparam Tp  The foating point type of the argument and stepsize.
   *
   * @param func The function to be differentated.
   * @param x The location at which the derivative is required.
   * @param h The stepsize.
   */
  template<typename Func, typename Tp>
    derivative_t<Tp>
    derivative_central(Func func, Tp x, Tp h)
    {
      const auto s_eps = std::numeric_limits<Tp>::epsilon();

      const auto fmh = func(x - h);
      const auto fph = func(x + h);
      const auto r3 = (fph - fmh) / (Tp{2} * h);

      const auto fmhh = func(x - h / Tp{2});
      const auto fphh = func(x + h / Tp{2});
      const auto r5 = (Tp{4} / Tp{3}) * (fphh - fmhh) / h
		      - (Tp{1} / Tp{3}) * r3;

      const auto e3 = (std::abs(fph) + std::abs(fmh)) * s_eps;
      const auto e5 = Tp{2} * (std::abs(fphh)
				+ std::abs(fmhh)) * s_eps + e3;
      // Rounding error:
      const auto dy = std::max(std::abs(r3), std::abs(r5))
      		      * std::abs(x / h) * s_eps;
      const auto err_round = std::abs(e5 / h) + dy;

      return {r5, std::abs(r5 - r3), err_round};
    }

  /**
   * Compute the derivative using the 4-point rule (x + h/4, x + h/2,
   * x + 3h/4, x + h).
   * @f[
   *    f'_2(x) = \frac{f(x + h) - f(x + h/2)}{h/2}
   * @f]
   * @f[
   *    f'_4(x) = \frac{22(f(x + h) - f(x + 3h/4))
   *                  - 62(f(x + 3h/4) - f(x + h/2))
   *                  + 52(f(x + h/2 - f(x + h/4))}{3h}
   * @f]
   *
   * Compute the error using the difference between the 4-point and
   * the 2-point rule (x + h/2, x+h).
   *
   * @tparam Func A function type callable with a numeric type
   *          and returning the same.
   * @tparam Tp  The foating point type of the argument and stepsize.
   *
   * @param func The function to be differentated.
   * @param x The location at which the derivative is required.
   * @param h The stepsize.
   */
  template<typename Func, typename Tp>
    derivative_t<Tp>
    derivative_forward(Func func, Tp x, Tp h)
    {
      const auto s_eps = std::numeric_limits<Tp>::epsilon();

      const auto f1 = func(x + h / Tp{4});
      const auto f2 = func(x + h / Tp{2});
      const auto f3 = func(x + Tp{3} * h / Tp{4});
      const auto f4 = func(x + h);

      const auto r2 = Tp{2} * (f4 - f2) / h;
      const auto r4 = (Tp{22} * (f4 - f3)
		       - Tp{62} * (f3 - f2)
		       + Tp{52} * (f2 - f1)) / (Tp{3} * h);


      const auto e4 = Tp{2 * 20.67}
		      * (std::abs(f4) + std::abs(f3)
		       + std::abs(f2) + std::abs(f1)) * s_eps;

      // The next term is due to finite precision in x+h = O (eps * x).
      const auto dy = std::max(std::abs(r2), std::abs(r4))
		      * std::abs(x / h) * s_eps;
      const auto err_round = std::abs(e4 / h) + dy;

      return {r4, std::abs(r4 - r2), err_round};
    }

  /**
   * Compute the derivative of a function func at a point x by Ridder's method
   * of polynomial extrapolation.
   * @f[
   *
   * @f]
   *
   * The value h is input as an estimated stepsize; it should not be small
   * but rather it should be an interval over which the function changes
   * substantially.
   *
   * @tparam Func A function type callable with a numeric type
   *          and returning the same.
   * @tparam Tp  The foating point type of the argument and stepsize.
   *
   * @param func The function to be differentated.
   * @param x The location at which the derivative is required.
   * @param h The initial stepsize.
   */
  template<typename Func, typename Tp>
    derivative_t<Tp>
    derivative_ridder(Func func, Tp x, Tp h)
    {
      const auto s_eps = std::numeric_limits<Tp>::epsilon();

      constexpr Tp s_scale = Tp{1.4};
      constexpr Tp s_scale2 = s_scale * s_scale;
      constexpr Tp s_big = std::numeric_limits<Tp>::max();
      constexpr Tp s_safe = Tp{2};

      if (h <= Tp{0})
	std::throw_domain_error("derivative_ridder:"
				  " Stepsize must be positive.");

      constexpr int _Num = 10;
      std::array<std::array<Tp, _Num>, _Num> a;
      std::array<Tp, _Num> rerr, dy;
      rerr.fill(Tp{0});
      dy.fill(Tp{0});

      auto hh = h;
      auto fph = func(x + hh);
      auto fmh = func(x - hh);

      a[0][0] = (fph - fmh) / (Tp{2} * hh);
      rerr[0] += (std::abs(fph) + std::abs(fmh)) * s_eps;
      dy[0] += std::abs(a[0][0]);

      auto ans = Tp{0};
      auto err_trunc = s_big;
      auto err_round = s_big;
      // Successive columns of the Neville tableau will go to smaller stepsizes
      // and to higher orders of extrapolation.
      for (int j = 1; j < _Num; ++j)
	{
	  // Try a new, smaller stepsize.
	  hh /= s_scale;
	  fph = func(x + hh);
	  fmh = func(x - hh);
	  a[0][j] = (fph - fmh) / (Tp{2} * hh);
	  rerr[j] += (std::abs(fph) + std::abs(fmh)) * s_eps;
	  dy[j] = std::abs(a[0][j]);//std::max(dy[j], std::abs(a[0][j]));
	  auto fac = Tp{1};
	  for (int i = 1; i <= j; ++i)
	    {
	      // Compute extrapolations of various orders, requiring
	      // no new function evaluations.
	      fac *= s_scale2;
	      a[i][j] = (fac * a[i-1][j] - a[i-1][j-1])
			    / (fac - Tp{1});
	      // Compare each new extrapolation to one order lower,
	      // both at the present stepsize and to the previous one.
	      auto errt = std::max(std::abs(a[i][j] - a[i-1][j]),
				   std::abs(a[i][j] - a[i-1][j-1]));
	      if (errt <= err_trunc)
		{
		  ans = a[i][j];
		  err_trunc = errt;
		  err_round = rerr[j] / h
			      + dy[j] * std::abs(x / h) * s_eps;
		}
	    }

	  // Quit if higher order is worse by a significant factor Safe.
	  if (std::abs(a[j][j])
	    - std::abs(a[j - 1][j - 1]) >= s_safe * err_trunc)
	    break;
	}

      // @todo Figure out a rounding error for the Ridder derivative.
      // e[j] = (std::abs(fp) + std::abs(fm)) * s_eps;
      return {ans, err_trunc, err_round};
    }

  /**
   * Compute the derivative of a function func at a point x by automatic
   * differentiation:
   * @f[
   *    f'(x) = Im[f(x + i\epsilon)] / \epsilon
   * @f]
   * This routine ignores the input stepsize and uses machine epsilon.
   *
   * If your function has complex overloads and works at least very near
   * the real axis this is a very accurate routine.
   *
   * @tparam Func A function type callable with a numeric type
   *          and returning the same.
   * @tparam Tp  The foating point type of the argument and stepsize.
   *
   * @param func The function to be differentated.
   * @param x The (real) location at which the derivative is required.
   * @param h The initial stepsize.
   */
  template<typename Func, typename Tp>
    derivative_t<Tp>
    derivative_automatic(Func func, Tp x, Tp)
    {
      using _Cmplx = std::complex<Tp>;
      constexpr auto s_i = _Cmplx{0, 1};
      constexpr auto s_eps = std::numeric_limits<Tp>::epsilon();
      const auto w = x + s_i * s_eps;
      const auto ans = std::imag(func(w)) / s_eps;
      // The answer involves a function eval and a division - zero roundoff.
      return {ans, s_eps, Tp{0}};
    }

#endif // DIFFERENTIATION_TCC
