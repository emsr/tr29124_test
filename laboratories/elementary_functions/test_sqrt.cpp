/**
 *
 */

#include <iostream>
#include <iomanip>
#include <cmath>
#include <limits>

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
	auto sqrt = std::sqrt(x);
	auto delta = sqrtr - sqrt;
	std::cout << ' ' << std::setw(w) << x
		  << ' ' << std::setw(w) << sqrtr
		  << ' ' << std::setw(w) << sqrt
		  << ' ' << std::setw(w) << delta
		  << '\n';
      }
  }

int
main()
{
  test_sqrt<double>();
}
