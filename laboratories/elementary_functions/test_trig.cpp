/**
 *
 */

#include <cmath>

#include <emsr/math_constants.h>

/**
 * naive Taylor-series implementation of sin(x).
 */
template<typename _Tp>
  _Tp
  sin_taylor(_Tp x)
  {
    if (fp_is_zero(x))
      return _Tp{0};
    else
      {
	const auto xx = x * x;
	auto term = _Tp{1};
	auto sum = _Tp{1};
	for (auto k = 0; k < 100; ++k)
	  {
	    term *= -xx / (2 * k + 1) / (2 * k + 2);
	    sum += term;
	  }
	return x * sum;
      }
  }

template<typename _Tp>
  _Tp
  sin(_Tp x)
  {
    const auto _S_2pi = emsr::tau_v<_Tp>;
    const auto _S_pi = emsr::pi_v<_Tp>;
    const auto _S_pi_2 = emsr::pi_v<_Tp> / _Tp{2};
    const auto _S_pi_4 = emsr::pi_v<_Tp> / _Tp{4};
    int sgn = 1;
    if (x < _Tp{0})
      {
	x = -x;
	sgn = -1;
      }
    x = fmod(x, _S_2pi);
    if (x > _S_pi)
      {
	x -= _S_pi;
	sgn *= -1;
      }
    if (x < _S_pi_4)
	return sin_taylor(x);
    else if (x < _S_pi_2)
	return cos_taylor(x);
  }

int
main()
{
}

