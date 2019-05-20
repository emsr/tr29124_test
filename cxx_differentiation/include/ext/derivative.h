#ifndef DERIVATIVE_H
#define DERIVATIVE_H 1

  template<typename _Tp>
    struct derivative_t
    {
      _Tp value;
      _Tp error_trunc;
      _Tp error_round;
    };

  /**
   * Compute the derivative using the 5-point rule
   *   (x-h, x-h/2, x, x+h/2, x+h).
   * Note that the central point is not used.
   *
   * Compute the error using the difference between the 5-point
   * and the 3-point rule (x-h, x, x+h).
   * Again the central point is not used.
   */
  template<typename _Func, typename _Tp>
    derivative_t<_Tp>
    derivative_central(_Func __func, _Tp __x, _Tp __h);

  /**
   * Compute the derivative using the 4-point rule (x+h/4, x+h/2,
   * x+3h/4, x+h).
   *
   * Compute the error using the difference between the 4-point and
   * the 2-point rule (x+h/2, x+h).
   */
  template<typename _Func, typename _Tp>
    derivative_t<_Tp>
    derivative_forward(_Func __func, _Tp __x, _Tp __h);

  template<typename _Func, typename _Tp>
    inline derivative_t<_Tp>
    derivative_backward(_Func __func, _Tp __x, _Tp __h)
    { return derivative_forward(__func, __x, -__h); }

  /*
   * Returns the derivative of a function func at a point x by Ridder's method
   * of polynomial extrapolation.
   * The value h is input as an estimated stepsize; it should not be small
   * but rather it should be an interval over which func changes substantially.
   */
  template<typename _Func, typename _Tp>
    derivative_t<_Tp>
    derivative_ridder(_Func __func, _Tp __x, _Tp __h);

#include <ext/derivative.tcc>

#endif // DERIVATIVE_H
