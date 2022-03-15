/**
 *
 */

#include <cmath>
#include <complex>
#include <limits>
#include <iostream>
#include <iomanip>
#include <emsr/math_constants.h>

template<typename Tp>
  struct Kelvin
  {
    Tp x;
    Tp ber_value, ber_deriv;
    Tp bei_value, bei_deriv;
    Tp ker_value, ker_deriv;
    Tp kei_value, kei_deriv;
  };

  /**
   */
  template<typename Tp>
    Kelvin<Tp>
    kelvin_asymp(Tp x)
    {
      constexpr auto s_sqrt2 = emsr::sqrt2_v<Tp>;
      constexpr auto s_pi_4 = emsr::pi_v<Tp> / Tp{4};
      constexpr auto s_pi_8 = emsr::pi_v<Tp> / Tp{8};
      constexpr auto s_mix_iter = 100;

      const auto xrt2 = x / s_sqrt2;
      const auto exp = std::exp(xrt2);
      const auto i8x = Tp{1} / Tp{8} / x;
      const auto dph = std::polar(Tp{1}, s_pi_4);
      auto sign = Tp{1};
      auto term = std::complex<Tp>{1};
      auto be_sum = std::complex<Tp>{0};
      auto be_sump = std::complex<Tp>{0};
      auto ke_sum = std::complex<Tp>{0};
      auto ke_sump = std::complex<Tp>{0};
      for (int k = 1; k < s_mix_iter; ++k)
	{
	  const auto __2km1 = Tp(2 * k - 1);
	  term *= (Tp{1} - __2km1) * (Tp{1} + __2km1) * i8x * dph / Tp(k);
	  be_sum += term;
	  be_sump += -Tp(k) * term / x;
	  sign *= Tp{-1};
	  ke_sum += sign * term;
	  ke_sump += -Tp(k) * sign * term / x;

	}
      const auto _P0 = Tp{1} + be_sum;
      const auto _Q0 = be_sum;
      const auto be_arg = xrt2 - s_pi_8;
      const auto ke_arg = xrt2 - s_pi_8;

      return {0, x,
	      ber, berp, bei, beip,
	      ker, kerp, kei, keip};
    }

int
main()
{
  for (int i = 100; i < 400; ++i)
    {
      const auto x = i * double(0.01);
      const auto kvn = kelvin_asymp(x);
    }
}
