/**
 *
 */

#include <cmath>
#include <complex>
#include <limits>
#include <iostream>
#include <iomanip>
#include <emsr/math_constants.h>

template<typename _Tp>
  struct Kelvin
  {
    _Tp x;
    _Tp ber_value, ber_deriv;
    _Tp bei_value, bei_deriv;
    _Tp ker_value, ker_deriv;
    _Tp kei_value, kei_deriv;
  };

  /**
   */
  template<typename _Tp>
    Kelvin<_Tp>
    kelvin_asymp(_Tp x)
    {
      constexpr auto s_sqrt2 = emsr::sqrt2_v<_Tp>;
      constexpr auto s_pi_4 = emsr::pi_v<_Tp> / _Tp{4};
      constexpr auto s_pi_8 = emsr::pi_v<_Tp> / _Tp{8};
      constexpr auto s_mix_iter = 100;

      const auto xrt2 = x / s_sqrt2;
      const auto exp = std::exp(xrt2);
      const auto i8x = _Tp{1} / _Tp{8} / x;
      const auto dph = std::polar(_Tp{1}, s_pi_4);
      auto sign = _Tp{1};
      auto term = std::complex<_Tp>{1};
      auto be_sum = std::complex<_Tp>{0};
      auto be_sump = std::complex<_Tp>{0};
      auto ke_sum = std::complex<_Tp>{0};
      auto ke_sump = std::complex<_Tp>{0};
      for (int k = 1; k < s_mix_iter; ++k)
	{
	  const auto __2km1 = _Tp(2 * k - 1);
	  term *= (_Tp{1} - __2km1) * (_Tp{1} + __2km1) * i8x * dph / _Tp(k);
	  be_sum += term;
	  be_sump += -_Tp(k) * term / x;
	  sign *= _Tp{-1};
	  ke_sum += sign * term;
	  ke_sump += -_Tp(k) * sign * term / x;

	}
      const auto _P0 = _Tp{1} + be_sum;
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
      const auto kvn = __kelvin_asymp(x);
    }
}
