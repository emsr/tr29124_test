#ifndef DERIVATIVE_TCC
#define DERIVATIVE_TCC 1

#include <cmath>
#include <limits>
#include <array>

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
    derivative_central(_Func __func, _Tp __x, _Tp __h)
    {
      const auto _S_eps = std::numeric_limits<_Tp>::epsilon();

      const auto __fmh = __func(__x - __h);
      const auto __fph = __func(__x + __h);
      const auto __r3 = (__fph - __fmh) / _Tp{2};

      const auto __fmhh = __func(__x - __h / _Tp{2});
      const auto __fphh = __func(__x + __h / _Tp{2});
      const auto __r5 = (_Tp{4} / _Tp{3}) * (__fphh - __fmhh)
		      - (_Tp{1} / _Tp{3}) * __r3;

      const auto __e3 = (std::abs(__fph) + std::abs(__fmh)) * _S_eps;
      const auto __e5 = _Tp{2} * (std::abs(__fphh)
				+ std::abs(__fmhh)) * _S_eps + __e3;
      // Rounding error:
      const auto __dy = std::max(std::abs(__r3 / __h), std::abs(__r5 / __h))
      		      * (std::abs(__x) / __h) * _S_eps;
      const auto __err_round = std::abs(__e5 / __h) + __dy;

      return {__r5 / __h, std::abs((__r5 - __r3) / __h), __err_round};
    }

  /**
   * Compute the derivative using the 4-point rule (x + h/4, x + h/2,
   * x + 3h/4, x + h).
   *
   * Compute the error using the difference between the 4-point and
   * the 2-point rule (x + h/2, x+h).
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

      const auto __r2 = _Tp{2} * (__f4 - __f2);
      const auto __r4 = (_Tp{22} / _Tp{3}) * (__f4 - __f3)
		      - (_Tp{62} / _Tp{3}) * (__f3 - __f2)
		      + (_Tp{52} / _Tp{3}) * (__f2 - __f1);


      const auto __e4 = _Tp{2 * 20.67}
		      * (std::abs(__f4) + std::abs(__f3)
		       + std::abs(__f2) + std::abs(__f1)) * _S_eps;

      // The next term is due to finite precision in x+h = O (eps * x).
      const auto __dy = std::max(std::abs(__r2 / __h), std::abs(__r4 / __h))
		      * std::abs(__x / __h) * _S_eps;
      const auto __err_round = std::abs(__e4 / __h) + __dy;

      return {__r4 / __h, std::abs((__r4 - __r2) / __h), __err_round};
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
   *
   * The value h is input as an estimated stepsize; it should not be small
   * but rather it should be an interval over which the function changes
   * substantially.
   */
  template<typename _Func, typename _Tp>
    derivative_t<_Tp>
    derivative_ridder(_Func __func, _Tp __x, _Tp __h)
    {
      constexpr _Tp _S_scale = _Tp{1.4};
      constexpr _Tp _S_scale2 = _S_scale * _S_scale;
      constexpr _Tp _S_big = _Tp{1.0e30};
      constexpr _Tp _S_safe = _Tp{2};

      if (__h <= _Tp{0})
	std::__throw_domain_error("derivative_ridder:"
				  " Stepsize must be positive.");

      constexpr int _Num = 10;
      std::array<std::array<_Tp, _Num>, _Num> __a;

      auto __hh = __h;
      __a[0][0] = (__func(__x + __hh) - __func(__x - __hh)) / (_Tp{2} * __hh);

      auto __ans = _Tp{0};
      auto __err = _S_big;
      // Successive columns of the Neville tableau will go to smaller stepsizes
      // and to higher orders of extrapolation.
      for (int __j = 1; __j < _Num; ++__j)
	{
	  // Try a new, smaller stepsize.
	  __hh /= _S_scale;
	  __a[0][__j] = (__func(__x + __hh) - __func(__x - __hh))
		      / (_Tp{2} * __hh);
	  auto __fac = _S_scale2;
	  for (int __i = 1; __i <= __j; ++__i)
	    {
	      // Compute extrapolations of various orders, requiring
	      // no new function evaluations.
	      __a[__i][__j] = (__fac * __a[__i-1][__j] - __a[__i-1][__j-1])
			    / (__fac - _Tp{1});
	      __fac *= _S_scale2;
	      // Compare each new extrapolation to one order lower,
	      // both at the present stepsize and to the previous one.
	      auto __errt = std::max(std::abs(__a[__i][__j] - __a[__i-1][__j]),
				   std::abs(__a[__i][__j] - __a[__i-1][__j-1]));
	      if (__errt <= __err)
		{
		  __err = __errt;
		  __ans = __a[__i][__j];
		}
	    }

	  // Quit if higher order is worse by a significant factor Safe.
	  if (std::abs(__a[__j][__j])
	    - std::abs(__a[__j - 1][__j - 1]) >= _S_safe * __err)
	    break;
	}

      // @todo Figure out a rounding error for the Ridder derivative.
      // __e[__j] = (std::abs(__fp) + std::abs(__fm)) * _S_eps;
      return {__ans, __err, _Tp{0}};
    }

#endif // DERIVATIVE_TCC
