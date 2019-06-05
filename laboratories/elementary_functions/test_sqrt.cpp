/**
 *
 */

#include <iostream>
#include <cmath>
#include <limits>

/**
 * 
 */
template<typename _Tp>
  _Tp
  __sqrt_recip_newton(_Tp __x)
  {
    if (__x == _Tp{0})
      return _Tp{0};

    auto __f = 0.01;
    auto __ff = __f * __f;
    while (std::abs(__ff - __x) > std::numeric_limits<_Tp>::epsilon() * __f)
      {
	__f *= (_Tp{3} - __x * __ff) / _Tp{2};
	__ff = __f * __f;
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
    for (int i = 0; i <= 1000; ++i)
      {
	auto x = _Tp{1.0} * i;
	auto sqrtr = x * __sqrt_recip_newton(x);
	std::cout << ' ' << x << ' ' << sqrtr << '\n';
      }
  }

int
main()
{
  test_sqrt<double>();
}
