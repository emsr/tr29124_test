/**
 *
 */

#include <cmath>

/**
 * naive Taylor-series implementation of sin(x).
 */
template<typename _Tp>
  _Tp
  __sin_taylor(_Tp __x)
  {
    if (__fp_is_zero(__x))
      return _Tp{0};
    else
      {
	const auto __xx = __x * __x;
	auto __term = _Tp{1};
	auto __sum = _Tp{1};
	for (auto __k = 0; __k < 100; ++__k)
	  {
	    __term *= -__xx / (2 * __k + 1) / (2 * __k + 2);
	    __sum += __term;
	  }
	return __x * __sum;
      }
  }

template<typename _Tp>
  _Tp
  __sin(_Tp __x)
  {
    const auto _S_2pi = emsr::tau_v<_Tp>;
    const auto _S_pi = emsr::pi_v<_Tp>;
    const auto _S_pi_2 = emsr::pi_v<_Tp> / _Tp{2};
    const auto _S_pi_4 = emsr::pi_v<_Tp> / _Tp{4};
    int __sgn = 1;
    if (__x < _Tp{0})
      {
	__x = -__x;
	__sgn = -1;
      }
    __x = fmod(__x, _S_2pi);
    if (__x > _S_pi)
      {
	__x -= _S_pi;
	__sgn *= -1;
      }
    if (__x < _S_pi_4)
	return __sin_taylor(__x);
    else if (__x < _S_pi_2)
	return __cos_taylor(__x);
  }

int
main()
{
}

