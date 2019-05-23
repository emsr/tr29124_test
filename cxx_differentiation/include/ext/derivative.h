#ifndef DERIVATIVE_H
#define DERIVATIVE_H 1

  /**
   * A struct conaining the results of a numerical differentiation.
   */
  template<typename _Tp>
    struct derivative_t
    {
      /// The value of the derivative.
      _Tp value;
      /// The absolute value of the truncation error.
      _Tp error_trunc;
      /// The absolute value of the rounding error.
      _Tp error_round;
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
    derivative_central(_Func __func, _Tp __x, _Tp __h);

  /**
   * Compute the derivative using the 4-point rule (x + h/4, x + h/2,
   * x + 3h/4, x + h).
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
    derivative_forward(_Func __func, _Tp __x, _Tp __h);

  template<typename _Func, typename _Tp>
    inline derivative_t<_Tp>
    derivative_backward(_Func __func, _Tp __x, _Tp __h)
    { return derivative_forward(__func, __x, -__h); }

  /**
   * Compute the derivative of a function func at a point x by Ridder's method
   * of polynomial extrapolation.
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
    derivative_ridder(_Func __func, _Tp __x, _Tp __h);

#include <ext/derivative.tcc>

#endif // DERIVATIVE_H
