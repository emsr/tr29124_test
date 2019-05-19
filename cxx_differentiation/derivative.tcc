#ifndef DERIVATIVE_TCC
#define DERIVATIVE_TCC 1

#include <cmath>

template<typename _Tp>
  struct derivative_t
  {
    _Tp value;
    _Tp error;
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
  derivative_central(_Func __func, _Tp __x, _Tp __h)
  {
    const auto __fm1 = (__x - __h);
    const auto __fp1 = (__x + __h);

    const auto __fmh = (__x - __h / _Tp{2});
    const auto __fph = (__x + __h / _Tp{2});

    const auto __r3 = (__fp1 - __fm1) / _Tp{2};
    const auto __r5 = (_Tp{4} / _Tp{3}) * (__fph - __fmh) - (_Tp{1} / _Tp{3}) * __r3;

    const auto __e3 = (std::abs(__fp1) + std::abs(__fm1)) * GSL_DBL_EPSILON;
    const auto __e5 = _Tp{2} * (std::abs(__fph) + std::abs(__fmh)) * GSL_DBL_EPSILON + __e3;
  }

/**
 * Compute the derivative using the 4-point rule (x+h/4, x+h/2,
 * x+3h/4, x+h).
 *
 * Compute the error using the difference between the 4-point and
 * the 2-point rule (x+h/2, x+h).
 */
template<typename _Func, typename _Tp>
  derivative_t<_Tp>
  derivative_forward(_Func __func, _Tp __x, _Tp __h)
  {
    const auto __f1 = __func(__x + __h / _Tp{4});
    const auto __f2 = __func(__x + __h / _Tp{2});
    const auto __f3 = __func(__x + _Tp{3} * __h / _Tp{4});
    const auto __f4 = __func(__x + __h);

    const auto __r2 = _Tp{2} * (__f4 - __f2);
    const auto __r4 = (_Tp{22} / _Tp{3}) * (__f4 - __f3)
		    - (_Tp{62} / _Tp{3}) * (__f3 - __f2)
		    + (_Tp{52} / _Tp{3}) * (__f2 - __f1);

    return {__r4 / __h, 
  }

template<typename _Func, typename _Tp>
  derivative_t<_Tp>
  derivative_backward(_Func __func, _Tp __x, _Tp __h)
  { return derivative_forward (__func, __x, -__h); }

/*
 * Returns the derivative of a function func at a point x by Ridder's method
 * of polynomial extrapolation.
 * The value h is input as an estimated stepsize; it should not be small
 * but rather it should be an interval over which func changes substantially.
 */
template<typename _Func, typename _Tp>
  derivative_t<_Tp>
  derivative_ridder(_Func __func, _Tp __x, _Tp __h)
  {
    constexpr int _NTAB = 10;
    constexpr _Tp _CON = 1.4, _CON2 = _CON * _CON;
    constexpr _Tp _BIG = 1.0e30;
    constexpr _Tp _SAFE = _Tp{2};

    if (h <= _Tp{0})
      std::__throw_domain_error("derivative_ridder:"
				" Stepsize must be positive.");

    std::vector<std::vector<_Tp>> __a(_NTAB, std::vector<_Tp>(_NTAB));

    auto __hh = __h;
    __a[0][0] = (__func(__x+__hh) - __func(__x-__hh)) / (_Tp{2} * __hh);

    auto __ans = _Tp{0};
    auto __err = _BIG;
    // Successive columns of the Neville tableau will go to smaller stepsizes
    // and to higher orders of extrapolation.
    for (int __j = 1; __j < _NTAB; ++__j)
     {
	// Try a new, smaller stepsize.
	__hh /= _CON;
	__a[0][j] = (__func(__x + __hh)- __func(__x - __hh)) / (_Tp{2} * __hh);
	auto __fac = _CON2;
	for (int __i = 1; __i <= __j; ++__i)
	  {
	    // Compute extrapolations of various orders, requiring
	    // no new function evaluations.
	    __a[i][j] = (__fac * __a[__i-1][__j] - __a[__i-1][__j-1])
			/ (__fac - _Tp{1});
	    __fac *=_ CON2;
	    // The error strategy is to compare the each new extrapolation
	    // to one order lower, both at the present stepsize
	    // and to the previous one.
	    auto __errt = std::max(std::abs(__a[__i][__j] - a[__i-1][__j]),
				   std::abs(__a[__i][__j] - a[__i-1][__j-1]));
	    if (__errt <= __err)
	      {
		__err = __errt;
		__ans = __a[__i][__j];
	      }
	  }
	// Quit if higher order is worse by a significant factor SAFE.
	if (std::abs(__a[__j][__j])
	  - std::abs(__a[__j-1][__j-1]) >= _SAFE * __err)
	  break;
      }

    return {__ans, __err};
  }

#endif // DERIVATIVE_TCC
