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
   * @tparam _Func A function type callable with a numeric type
   *          and returning the same.
   * @tparam _Tp  The foating point type of the argument and stepsize.
   *
   * @param __func The function to be differentated.
   * @param __x The location at which the derivative is required.
   * @param __h The stepsize.
   */
  template<typename _Func, typename _Tp>
    derivative_t<_Tp>
    derivative_central(_Func __func, _Tp __x, _Tp __h)
    {
      const auto _S_eps = std::numeric_limits<_Tp>::epsilon();

      const auto __fmh = __func(__x - __h);
      const auto __fph = __func(__x + __h);
      const auto __r3 = (__fph - __fmh) / (_Tp{2} * __h);

      const auto __fmhh = __func(__x - __h / _Tp{2});
      const auto __fphh = __func(__x + __h / _Tp{2});
      const auto __r5 = (_Tp{4} / _Tp{3}) * (__fphh - __fmhh) / __h
		      - (_Tp{1} / _Tp{3}) * __r3;

      const auto __e3 = (std::abs(__fph) + std::abs(__fmh)) * _S_eps;
      const auto __e5 = _Tp{2} * (std::abs(__fphh)
				+ std::abs(__fmhh)) * _S_eps + __e3;
      // Rounding error:
      const auto __dy = std::max(std::abs(__r3), std::abs(__r5))
      		      * std::abs(__x / __h) * _S_eps;
      const auto __err_round = std::abs(__e5 / __h) + __dy;

      return {__r5, std::abs(__r5 - __r3), __err_round};
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
   * @tparam _Func A function type callable with a numeric type
   *          and returning the same.
   * @tparam _Tp  The foating point type of the argument and stepsize.
   *
   * @param __func The function to be differentated.
   * @param __x The location at which the derivative is required.
   * @param __h The stepsize.
   */
  template<typename _Func, typename _Tp>
    derivative_t<_Tp>
    derivative_forward(_Func __func, _Tp __x, _Tp __h)
    {
      const auto _S_eps = std::numeric_limits<_Tp>::epsilon();

      const auto __f1 = __func(__x + __h / _Tp{4});
      const auto __f2 = __func(__x + __h / _Tp{2});
      const auto __f3 = __func(__x + _Tp{3} * __h / _Tp{4});
      const auto __f4 = __func(__x + __h);

      const auto __r2 = _Tp{2} * (__f4 - __f2) / __h;
      const auto __r4 = (_Tp{22} * (__f4 - __f3)
		       - _Tp{62} * (__f3 - __f2)
		       + _Tp{52} * (__f2 - __f1)) / (_Tp{3} * __h);


      const auto __e4 = _Tp{2 * 20.67}
		      * (std::abs(__f4) + std::abs(__f3)
		       + std::abs(__f2) + std::abs(__f1)) * _S_eps;

      // The next term is due to finite precision in x+h = O (eps * x).
      const auto __dy = std::max(std::abs(__r2), std::abs(__r4))
		      * std::abs(__x / __h) * _S_eps;
      const auto __err_round = std::abs(__e4 / __h) + __dy;

      return {__r4, std::abs(__r4 - __r2), __err_round};
    }

/* optimize stepsize.
IMHO, this should take an initial stab at the stepsize.
  template<typename _Func, typename _Tp>
    derivative_t<_Tp>
    derivative_forward(_Func __func, _Tp __x)
    {
      auto [__result, __trunc, __round] = derivative_forward(__func, __x, __h);
      auto __error = __round + __trunc;

      if (__round < __trunc && (__round > 0 && __trunc > 0))
	{
	  // Compute an optimised stepsize to minimize the total error,
	  // using the scaling of the truncation error (O(h^2)) and
	  // rounding error (O(1/h)).
	  auto __h_opt = __h * std::cbrt(__round / (_Tp{2} * __trunc));
	  auto [__r_opt, __trunc_opt, __round_opt]
		 = derivative_forward(__func, __x, __h_opt);
	  auto __error_opt = __round_opt + __trunc_opt;

	  // Check that the new error is smaller, and that the new derivative 
	  // is consistent with the error bounds of the original estimate.
	  if (__error_opt < __error
	      && std::abs(__r_opt - __result) < _Tp{4} * __error)
	    {
	      __result = __r_opt;
	      __trunc = __trunc_opt;
	      __round = __round_opt;
	    }
	}

      return {__result, __trunc, __round};
    }
*/
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
   * @tparam _Func A function type callable with a numeric type
   *          and returning the same.
   * @tparam _Tp  The foating point type of the argument and stepsize.
   *
   * @param __func The function to be differentated.
   * @param __x The location at which the derivative is required.
   * @param __h The initial stepsize.
   */
  template<typename _Func, typename _Tp>
    derivative_t<_Tp>
    derivative_ridder(_Func __func, _Tp __x, _Tp __h)
    {
      const auto _S_eps = std::numeric_limits<_Tp>::epsilon();

      constexpr _Tp _S_scale = _Tp{1.4};
      constexpr _Tp _S_scale2 = _S_scale * _S_scale;
      constexpr _Tp _S_big = std::numeric_limits<_Tp>::max();
      constexpr _Tp _S_safe = _Tp{2};

      if (__h <= _Tp{0})
	std::__throw_domain_error("derivative_ridder:"
				  " Stepsize must be positive.");

      constexpr int _Num = 10;
      std::array<std::array<_Tp, _Num>, _Num> __a;
      std::array<_Tp, _Num> __rerr, __dy;
      __rerr.fill(_Tp{0});
      __dy.fill(_Tp{0});

      auto __hh = __h;
      auto __fph = __func(__x + __hh);
      auto __fmh = __func(__x - __hh);

      __a[0][0] = (__fph - __fmh) / (_Tp{2} * __hh);
      __rerr[0] += (std::abs(__fph) + std::abs(__fmh)) * _S_eps;
      __dy[0] += std::abs(__a[0][0]);

      auto __ans = _Tp{0};
      auto __err_trunc = _S_big;
      auto __err_round = _S_big;
      // Successive columns of the Neville tableau will go to smaller stepsizes
      // and to higher orders of extrapolation.
      for (int __j = 1; __j < _Num; ++__j)
	{
	  // Try a new, smaller stepsize.
	  __hh /= _S_scale;
	  __fph = __func(__x + __hh);
	  __fmh = __func(__x - __hh);
	  __a[0][__j] = (__fph - __fmh) / (_Tp{2} * __hh);
	  __rerr[__j] += (std::abs(__fph) + std::abs(__fmh)) * _S_eps;
	  __dy[__j] = std::abs(__a[0][__j]);//std::max(__dy[__j], std::abs(__a[0][__j]));
	  auto __fac = _Tp{1};
	  for (int __i = 1; __i <= __j; ++__i)
	    {
	      // Compute extrapolations of various orders, requiring
	      // no new function evaluations.
	      __fac *= _S_scale2;
	      __a[__i][__j] = (__fac * __a[__i-1][__j] - __a[__i-1][__j-1])
			    / (__fac - _Tp{1});
	      // Compare each new extrapolation to one order lower,
	      // both at the present stepsize and to the previous one.
	      auto __errt = std::max(std::abs(__a[__i][__j] - __a[__i-1][__j]),
				   std::abs(__a[__i][__j] - __a[__i-1][__j-1]));
	      if (__errt <= __err_trunc)
		{
		  __ans = __a[__i][__j];
		  __err_trunc = __errt;
		  __err_round = __rerr[__j] / __h
			      + __dy[__j] * std::abs(__x / __h) * _S_eps;
		}
	    }

	  // Quit if higher order is worse by a significant factor Safe.
	  if (std::abs(__a[__j][__j])
	    - std::abs(__a[__j - 1][__j - 1]) >= _S_safe * __err_trunc)
	    break;
	}

      // @todo Figure out a rounding error for the Ridder derivative.
      // __e[__j] = (std::abs(__fp) + std::abs(__fm)) * _S_eps;
      return {__ans, __err_trunc, __err_round};
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
   * @tparam _Func A function type callable with a numeric type
   *          and returning the same.
   * @tparam _Tp  The foating point type of the argument and stepsize.
   *
   * @param __func The function to be differentated.
   * @param __x The (real) location at which the derivative is required.
   * @param __h The initial stepsize.
   */
  template<typename _Func, typename _Tp>
    derivative_t<_Tp>
    derivative_automatic(_Func __func, _Tp __x, _Tp)
    {
      using _Cmplx = std::complex<_Tp>;
      constexpr auto _S_i = _Cmplx{0, 1};
      constexpr auto _S_eps = std::numeric_limits<_Tp>::epsilon();
      const auto __w = __x + _S_i * _S_eps;
      const auto __ans = std::imag(__func(__w)) / _S_eps;
      return {__ans, _S_eps, _S_eps};
    }

#endif // DIFFERENTIATION_TCC
