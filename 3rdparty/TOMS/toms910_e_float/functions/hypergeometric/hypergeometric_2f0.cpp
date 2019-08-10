
#include <e_float/e_float.h>
#include <functions/constants/constants.h>
#include <functions/elementary/elementary.h>
#include <functions/gamma/gamma.h>
#include <functions/hypergeometric/hypergeometric.h>
#include <functions/hypergeometric/laguerre.h>

e_float ef::hyperg_2f0(const e_float& a, const e_float& b, const e_float& x)
{
  const INT32 Am = ef::to_int32(a);
  const INT32 Bn = ef::to_int32(b);

  if(ef::isint(a) && (Am < static_cast<INT32>(1)))
  {
    // A Laguerre polynomial is used for zero or negative integer value of a.
    const INT32 n = static_cast<INT32>(-Am);

    const e_float factor = ef::pown(-x, static_cast<INT64>(n)) * ef::factorial(static_cast<UINT32>(n));

    const e_float result = factor * ef::laguerre(n, -b - n, ef::one_minus() / x);

    const bool b_neg = (static_cast<INT32>(n % static_cast<INT32>(2)) != static_cast<INT32>(0));

    return (!b_neg ? result : -result);
  }

  if(ef::isint(b) && (Bn < static_cast<INT32>(1)))
  {
    // A Laguerre polynomial is used for zero or negative integer value of b.
    // This is done by switching the order of a and b because the function is
    // symmetric with respect to a and b. The subsequent function call returns
    // the result of the Laguerre polynomial.

    return ef::hyperg_2f0(b, a, x);
  }

  // Hypergeometric2F0 is only implemented for |x| < 1.
  return ((ef::fabs(x) < ef::one()) ? ef::hyp2F0(a, b, x)
                                    : std::numeric_limits<e_float>::quiet_NaN());
}
