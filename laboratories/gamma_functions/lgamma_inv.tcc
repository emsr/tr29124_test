#include <emsr/fp_type_util.h>

#include <cmath>

namespace __detail
{
  /**
   * Invert @f$ log(\Gamma(x)) @f$ using Newton's rule.
   * To support ilfactorial you could set lsg = logGammaStar to zero.
   */
  template<typename _Tp>
    _Tp
    __lgamma_inv_newton(_Tp __y)
    {
      constexpr auto _S_eps = std::numeric_limits<_Tp>::epsilon();
      constexpr auto _S_ln_2 = emsr::ln2_v<_Tp>;
      constexpr auto _S_ln_pi = emsr::lnpi_v<_Tp>;
      constexpr auto _S_log_sqrt_2pi = (_S_ln_2 + _S_ln_pi) / _Tp{2};
      constexpr auto _S_max_iter = 100;
      auto __a0 = __y - _S_log_sqrt_2pi;
      auto __a = __a0;
      if (__a <= _Tp{0})
	__a = _Tp{2};
      for (int __i = 0; __i < _S_max_iter; ++__i)
	{
	  const auto __a_prev = __a;
	  const auto __lgs = std::lgamma(__a)
		   - (__a - _Tp{0.5L}) * std::log(__a) + __a - _S_log_sqrt_2pi;
	  __a = _Tp{0.5L} + (__a0 + __a - __lgs) / std::log(__a);
	  if (std::abs(__a - __a_prev) < std::abs(__a) * _S_eps)
	    break;
	}
      return __a;
    }
}

  template<typename _Tp>
    inline _Tp
    __lgamma_inv(_Tp __y)
    {
      if (std::isnan(__y))
	return __y;
      else
	return __detail::__lgamma_inv_newton(__y);
    };

  inline float
  lgamma_invf(float __y)
  { return __lgamma_inv<float>(__y); }

  inline long double
  lgamma_invl(long double __y)
  { return __lgamma_inv<long double>(__y); }

  // Do arg promotion
  template<typename _Tp>
    inline emsr::fp_promote_t<_Tp>
    lgamma_inv(_Tp __y)
    {
      using __type = emsr::fp_promote_t<_Tp>;
      return __lgamma_inv<__type>(__y);
    }
