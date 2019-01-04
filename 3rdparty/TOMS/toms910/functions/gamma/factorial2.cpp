
#include <functions/complex/e_float_complex.h>
#include <functions/elementary/elementary.h>
#include <functions/gamma/gamma.h>
#include <functions/tables/tables.h>

namespace Factorial2_Series
{
  static e_float AtInfinity(const INT32 n);
}

static e_float Factorial2_Series::AtInfinity(const INT32 n)
{
  const bool n_is_even = static_cast<INT32>(n % static_cast<INT32>(2)) == static_cast<INT32>(0u);

  if(n_is_even)
  {
    const UINT32 n_half = static_cast<UINT32>(static_cast<UINT32>(n) / static_cast<UINT32>(2));

    return ef::pow2(static_cast<INT64>(n_half)) * ef::factorial(n_half);
  }
  else
  {
    const INT32 n_plus_one = static_cast<INT32>(static_cast<INT32>(n) + static_cast<INT32>(1));

    return ef::factorial(static_cast<UINT32>(n_plus_one)) / ef::factorial2(n_plus_one);
  }
}

e_float ef::factorial2(const INT32 n)
{
  const bool n_is_neg  = (n < static_cast<INT32>(0));

  if(!n_is_neg)
  {
    return (static_cast<std::size_t>(n) < Tables::A006882().size()) ? Tables::A006882()[n]()
                                                                    : Factorial2_Series::AtInfinity(n);
  }
  else
  {
    if(n == static_cast<INT32>(-1))
    {
      return ef::one();
    }

    const INT32 nn        = (!n_is_neg ? n : static_cast<INT32>(-n));
    const bool  n_is_even = (static_cast<INT32>(nn % static_cast<INT32>(2)) == static_cast<INT32>(0u));

    if(n_is_even)
    {
      return std::numeric_limits<e_float>::quiet_NaN();
    }
    else
    {
      // Double factorial for negative odd integers.
      const INT32 n_minus_one_over_two = static_cast<INT32>(static_cast<INT32>(nn - static_cast<INT32>(1)) / static_cast<INT32>(2));
      const bool  b_negate             = static_cast<INT32>(n_minus_one_over_two % static_cast<INT32>(2)) != static_cast<INT32>(0);

      const e_float f2 = ef::one() / ef::factorial2(static_cast<INT32>(nn - static_cast<INT32>(2)));

      return (!b_negate ? f2 : -f2);
    }
  }
}
