#include <cmath>

#include <emsr/fp_type_util.h>
#include <emsr/math_constants.h>

namespace detail
{
  /**
   * Invert @f$ log(\Gamma(x)) @f$ using Newton's rule.
   * To support ilfactorial you could set lsg = logGammaStar to zero.
   */
  template<typename _Tp>
    _Tp
    lgamma_inv_newton(_Tp y)
    {
      constexpr auto s_eps = std::numeric_limits<_Tp>::epsilon();
      constexpr auto s_ln_2 = emsr::ln2_v<_Tp>;
      constexpr auto s_ln_pi = emsr::lnpi_v<_Tp>;
      constexpr auto s_log_sqrt_2pi = (s_ln_2 + s_ln_pi) / _Tp{2};
      constexpr auto s_max_iter = 100;
      auto a0 = y - s_log_sqrt_2pi;
      auto a = a0;
      if (a <= _Tp{0})
	a = _Tp{2};
      for (int i = 0; i < s_max_iter; ++i)
	{
	  const auto a_prev = a;
	  const auto lgs = std::lgamma(a)
		   - (a - _Tp{0.5L}) * std::log(a) + a - s_log_sqrt_2pi;
	  a = _Tp{0.5L} + (a0 + a - lgs) / std::log(a);
	  if (std::abs(a - a_prev) < std::abs(a) * s_eps)
	    break;
	}
      return a;
    }

  template<typename _Tp>
    inline _Tp
    lgamma_inv(_Tp y)
    {
      if (std::isnan(y))
	return y;
      else
	return detail::lgamma_inv_newton(y);
    };
}

  inline float
  lgamma_invf(float y)
  { return detail::lgamma_inv<float>(y); }

  inline long double
  lgamma_invl(long double y)
  { return detail::lgamma_inv<long double>(y); }

  // Do arg promotion
  template<typename _Tp>
    inline emsr::fp_promote_t<_Tp>
    lgamma_inv(_Tp y)
    {
      using type = emsr::fp_promote_t<_Tp>;
      return detail::lgamma_inv<type>(y);
    }
