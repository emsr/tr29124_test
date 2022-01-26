
// Copyright (C) 2019 Free Software Foundation, Inc.
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

#ifndef DIFFERENTIATION_H
#define DIFFERENTIATION_H 1

/**
 * @file  differentiation.h
 *
 * This file contains the declarations and the inline implementations
 * of the differentiation data structures and routines.
 */

/**
 * @mainpage The cxx_differentiation library
 *
 * \section intro_sec Introduction
 *
 * This is the introduction.
 *
 * \section install_sec Installation
 *
 * \subsection step1 Step 1: Opening the box
 *
 * etc...
 */

  /**
   * A struct conaining the results of a numerical differentiation.
   */
  template<typename Tp>
    struct derivative_t
    {
      /// The value of the derivative.
      Tp value;
      /// The absolute value of the truncation error.
      Tp error_trunc;
      /// The absolute value of the rounding error.
      Tp error_round;
    };

  /**
   * Compute the derivative using the 5-point rule
   * (x-h, x-h/2, x, x+h/2, x+h).
   * Note that the central point is not used.
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
    derivative_central(Func func, Tp x, Tp h);

  /**
   * Compute the derivative using the 4-point rule (x + h/4, x + h/2,
   * x + 3h/4, x + h).
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
    derivative_forward(Func func, Tp x, Tp h);


  /**
   * Compute the derivative using the 4-point rule (x - h/4, x - h/2,
   * x - 3h/4, x - h).
   *
   * Compute the error using the difference between the 4-point and
   * the 2-point rule (x - h/2, x - h).
   *
   * i.e. run the forward differentiator with negative stepsize.
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
    inline derivative_t<Tp>
    derivative_backward(Func func, Tp x, Tp h)
    { return derivative_forward(func, x, -h); }

  /**
   * Compute the derivative of a function func at a point x by Ridder's method
   * of polynomial extrapolation.
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
    derivative_ridder(Func func, Tp x, Tp h);

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
    derivative_automatic(Func func, Tp x, Tp);

#include <ext/differentiation.tcc>

#endif // DIFFERENTIATION_H
