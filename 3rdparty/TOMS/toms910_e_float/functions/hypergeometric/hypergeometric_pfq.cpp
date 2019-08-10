
#include <functions/constants/constants.h>
#include <functions/elementary/elementary.h>
#include <functions/hypergeometric/hypergeometric.h>

e_float ef::hyperg_pfq(const std::deque<e_float>& a, const  std::deque<e_float>& b, const e_float& x)
{
  // This is only implemented for |x| < 1.
  return ef::fabs(x) < ef::one() ? ef::hypPFQ(a, b, x)
                                 : std::numeric_limits<e_float>::quiet_NaN();
}
