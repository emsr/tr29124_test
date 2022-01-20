/**
 *
 */

#include <iostream>
#include <iomanip>
#include <cmath>
#include <limits>

#include <emsr/math_constants.h> // For math constants.

/**
 * 
 */
template<typename _Tp>
  constexpr _Tp
  __sqrt_newton(_Tp __x)
  {
    const auto _S_eps = std::numeric_limits<_Tp>::epsilon();
    const auto _S_digs = std::numeric_limits<_Tp>::digits;
    constexpr auto __rsqrt2 = _Tp{1} / emsr::sqrt2_v<_Tp>;

    if (__x == _Tp{0})
      return _Tp{0};

    int __e;
    auto __xx = std::frexp(__x, &__e);

    auto __f = 0.4 + 0.6 * __xx;

    for (int i = 0; i < _S_digs; ++i)
      {
	__f -= _Tp{0.5L} * (__f - __xx / __f);
	if (std::abs(__f * __f - __xx) < _S_eps * __f)
	  break;
      }
    if (__e & 1)
      return __rsqrt2 * std::ldexp(__f, (__e + 1) >> 1);
    else
      return std::ldexp(__f, __e >> 1);
  }

/**
 * 
 */
template<typename _Tp>
  constexpr _Tp
  __sqrt_recip_newton(_Tp __x)
  {
    const auto _S_eps = std::numeric_limits<_Tp>::epsilon();
    const auto _S_digs = std::numeric_limits<_Tp>::digits;

    if (__x == _Tp{0})
      return _Tp{0};

    // TODO: Argument reduction...

    auto __f = 0.01;
    auto __ff = __f * __f;
    for (int i = 0; i < _S_digs; ++i)
      {
	__f *= (_Tp{3} - __x * __ff) / _Tp{2};
	__ff = __f * __f;
	if (std::abs(__ff - __x) < _S_eps * __f)
	  break;
      }

    return __f;
  }

/**
 * 
 */
template<typename _Tp>
  void
  test_sqrt()
  {
    std::cout.precision(std::numeric_limits<_Tp>::max_digits10);
    const auto w = 6 + std::cout.precision();

    for (int i = 0; i <= 1000; ++i)
      {
	auto x = _Tp{1.0L} * i;
	auto sqrtr = x * __sqrt_recip_newton(x);
	auto sqrtcw = __sqrt_newton(x);
	auto sqrt = std::sqrt(x);
	auto deltar = sqrtr - sqrt;
	auto deltacw = sqrtcw - sqrt;
	std::cout << ' ' << std::setw(w) << x
		  << ' ' << std::setw(w) << sqrtr
		  << ' ' << std::setw(w) << sqrtcw
		  << ' ' << std::setw(w) << sqrt
		  << ' ' << std::setw(w) << deltar
		  << ' ' << std::setw(w) << deltacw
		  << '\n';
      }
  }

int
main()
{
  test_sqrt<double>();
}
