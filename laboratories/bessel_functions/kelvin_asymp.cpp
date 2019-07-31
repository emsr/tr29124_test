/**
 *
 */

#include <cmath>
#include <complex>
#include <limits>
#include <iostream>
#include <iomanip>
#include <ext/math_constants.h>

template<typename _Tp>
  struct Kelvin
  {
    _Tp __x;
    _Tp __ber_value, __ber_deriv;
    _Tp __bei_value, __bei_deriv;
    _Tp __ker_value, __ker_deriv;
    _Tp __kei_value, __kei_deriv;
  };

  /**
   */
  template<typename _Tp>
    Kelvin<_Tp>
    __kelvin_asymp(_Tp __x)
    {
      constexpr auto _S_sqrt2 = __gnu_cxx::numbers::__root_2_v<_Tp>;
      constexpr auto _S_pi_4 = __gnu_cxx::numbers::__pi_v<_Tp> / _Tp{4};
      constexpr auto _S_pi_8 = __gnu_cxx::numbers::__pi_v<_Tp> / _Tp{8};
      constexpr auto _S_mix_iter = 100;

      const auto __xrt2 = __x / _S_sqrt2;
      const auto __exp = std::exp(__xrt2);
      const auto __i8x = _Tp{1} / _Tp{8} / __x;
      const auto __dph = std::polar(_Tp{1}, _S_pi_4);
      auto __sign = _Tp{1};
      auto __term = std::complex<_Tp>{1};
      auto __be_sum = std::complex<_Tp>{0};
      auto __be_sump = std::complex<_Tp>{0};
      auto __ke_sum = std::complex<_Tp>{0};
      auto __ke_sump = std::complex<_Tp>{0};
      for (int __k = 1; __k < _S_mix_iter; ++__k)
	{
	  const auto __2km1 = _Tp(2 * __k - 1);
	  __term *= (_Tp{1} - __2km1) * (_Tp{1} + __2km1) * __i8x * __dph / _Tp(__k);
	  __be_sum += __term;
	  __be_sump += -_Tp(__k) * __term / __x;
	  __sign *= _Tp{-1};
	  __ke_sum += __sign * __term;
	  __ke_sump += -_Tp(__k) * __sign * __term / __x;

	}
      const auto _P0 = _Tp{1} + __be_sum;
      const auto _Q0 = __be_sum;
      const auto __be_arg = __xrt2 - _S_pi_8;
      const auto __ke_arg = __xrt2 - _S_pi_8;

      return {0, __x,
	      __ber, __berp, __bei, __beip,
	      __ker, __kerp, __kei, __keip};
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
