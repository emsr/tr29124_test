
#include <e_float/e_float.h>
#include <functions/elementary/elementary.h>
#include <functions/hypergeometric/hypergeometric.h>

e_float ef::hyperg_1f0(const e_float& a, const e_float& x)
{
  // Hypergeometric1F0 is only implemented for |x| < 1.
  return (ef::fabs(x) < ef::one() ? ef::hyp1F0(a, x)
                                  : std::numeric_limits<e_float>::quiet_NaN());
}
