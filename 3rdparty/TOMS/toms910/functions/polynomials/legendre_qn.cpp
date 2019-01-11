
#include <e_float/e_float.h>
#include <functions/constants/constants.h>
#include <functions/elementary/elementary.h>
#include <functions/gamma/gamma.h>
#include <functions/polynomials/polynomials.h>

namespace LegendreQ_Series
{
  static e_float AtIdenticallyZero(const INT32 n);
}

static e_float LegendreQ_Series::AtIdenticallyZero(const INT32 n)
{
  const INT32 n_is_odd = static_cast<INT32>(n % static_cast<INT32>(2)) != static_cast<INT32>(0);

  if(n_is_odd)
  {
    const INT32 n_plus_one          = static_cast<INT32>(n + static_cast<INT32>(1));
    const INT32 n_plus_one_over_two = static_cast<INT32>(n_plus_one / static_cast<INT32>(2));
    const bool  b_negate            = static_cast<INT32>(n_plus_one_over_two % static_cast<INT32>(2)) != static_cast<INT32>(0);

    const UINT32 n_minus_one = static_cast<UINT32>(static_cast<UINT32>(n) - static_cast<UINT32>(1u));

    const e_float Qn = ef::factorial2(n_minus_one) / ef::factorial2(static_cast<UINT32>(n));

    return !b_negate ? Qn : -Qn; 
  }
  else
  {
    return ef::zero();
  }
}

e_float ef::legendre_q(const INT32 n, const e_float& x)
{
  // Perform a range-check and an order-check.
  if(   (ef::isone(ef::abs(x)))
     || (x < ef::one_minus())
     || (x > ef::one())
     || (n < static_cast<INT32>(0))
    )
  {
    return std::numeric_limits<e_float>::quiet_NaN();
  }

  if(ef::iszero(x))
  {
    return LegendreQ_Series::AtIdenticallyZero(n);
  }

  if(ef::isneg(x))
  {
    // Use the reflection formula for negative argument.
    const bool b_negate = static_cast<INT32>(n  % static_cast<INT32>(2)) != static_cast<INT32>(0);

    const e_float Qx = ef::legendre_q(n, -x);

    return !b_negate ? -Qx : Qx;
  }

  // Compute the logarithmic part of LegendreQ.
  const e_float log_term = ef::small_arg(x) ?  ef::log1p1m2(x)
                                            : (ef::log((ef::one() + x) / (ef::one() - x)) / static_cast<INT32>(2));

  if(n == static_cast<INT32>(0))
  {
    // Compute LegendreQ of order 0.
    return log_term;
  }

  if(n == static_cast<INT32>(1))
  {
    // Compute LegendreQ of order 1.
    return (legendre_p(static_cast<INT32>(1), x) * log_term) - ef::one();
  }

  // Use forward recursion to calculate Legendre Polynomials of degree 2 or higher.
  e_float Qn_m2 = ef::legendre_q(static_cast<INT32>(0), x);
  e_float Qn_m1 = ef::legendre_q(static_cast<INT32>(1), x);
  e_float Qn;

  for(INT32 j = static_cast<INT32>(2); j <= n; j++)
  {
    const INT32 n_minus_one     = static_cast<INT32>(j - static_cast<INT32>(1));
    const INT32 two_n_minus_one = static_cast<INT32>(n_minus_one + j);

    Qn = (((two_n_minus_one * x) * Qn_m1) - (n_minus_one * Qn_m2)) / j;

    Qn_m2 = Qn_m1;
    Qn_m1 = Qn;
  }

  return Qn;
}
