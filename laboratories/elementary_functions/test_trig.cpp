/**
 *
 */

#include <cmath>

#include <emsr/math_constants.h>

/**
 * naive Taylor-series implementation of sin(x).
 */
template<typename Tp>
  Tp
  sin_taylor(Tp x)
  {
    if (fp_is_zero(x))
      return Tp{0};
    else
      {
	const auto xx = x * x;
	auto term = Tp{1};
	auto sum = Tp{1};
	for (auto k = 0; k < 100; ++k)
	  {
	    term *= -xx / (2 * k + 1) / (2 * k + 2);
	    sum += term;
	  }
	return x * sum;
      }
  }

template<typename Tp>
  Tp
  sin(Tp x)
  {
    const auto s_2pi = emsr::tau_v<Tp>;
    const auto s_pi = emsr::pi_v<Tp>;
    const auto s_pi_2 = emsr::pi_v<Tp> / Tp{2};
    const auto s_pi_4 = emsr::pi_v<Tp> / Tp{4};
    int sgn = 1;
    if (x < Tp{0})
      {
	x = -x;
	sgn = -1;
      }
    x = fmod(x, s_2pi);
    if (x > s_pi)
      {
	x -= s_pi;
	sgn *= -1;
      }
    if (x < s_pi_4)
	return sin_taylor(x);
    else if (x < s_pi_2)
	return cos_taylor(x);
  }

int
main()
{
}

