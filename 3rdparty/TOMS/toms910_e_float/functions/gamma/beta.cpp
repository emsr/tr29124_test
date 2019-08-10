#include <functions/complex/e_float_complex.h>
#include <functions/constants/constants.h>
#include <functions/elementary/elementary.h>
#include <functions/gamma/gamma.h>
#include <functions/hypergeometric/hypergeometric.h>

namespace Beta_Series
{
  template<typename T> inline T Beta_Template(const T& a, const T& b)
  {
    using ef::gamma;
    using efz::gamma;

    return (gamma(a) * gamma(b)) / gamma(a + b);
  }
}

e_float ef::beta(const e_float& a, const e_float& b)
{
  return Beta_Series::Beta_Template<e_float>(a, b);
}

ef_complex efz::beta(const ef_complex& a, const ef_complex& b)
{
  return Beta_Series::Beta_Template<ef_complex>(a, b);
}

e_float ef::incomplete_beta(const e_float& x, const e_float& a, const e_float& b)
{
  if(ef::iszero(x))
  {
    return ef::zero();
  }
  else if(ef::isone(x))
  {
    return ef::beta(a, b);
  }
  else if((x < ef::zero()) || (x > ef::one()))
  {
    return std::numeric_limits<e_float>::quiet_NaN();
  }
  else
  {
    const e_float x_pow_a = ef::pow(x, a);

    if(x < ef::half())
    {
      return (x_pow_a / a) * ef::hyperg_2f1(a, ef::one() - b, a + ef::one(), x);
    }
    else
    {
      const e_float one_minus_x       = ef::one() - x;
      const e_float one_minus_x_pow_b = ef::pow(one_minus_x, b);

      const e_float term1 = ef::beta(a, b);
      const e_float term2 = ((one_minus_x_pow_b * x_pow_a) / b) * ef::hyperg_2f1(ef::one(), a + b, b + ef::one(), one_minus_x);

      return term1 - term2;
    }
  }
}
