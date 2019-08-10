
#include <e_float/e_float.h>
#include <functions/elementary/elementary.h>
#include <functions/hypergeometric/hypergeometric.h>

e_float ef::hyperg_0f0(const e_float& x)
{
  // Hypergeometric0F0 is only implemented for |x| < 1.
  return (ef::fabs(x) < ef::one() ? ef::hyp0F0(x)
                                  : std::numeric_limits<e_float>::quiet_NaN());
}
