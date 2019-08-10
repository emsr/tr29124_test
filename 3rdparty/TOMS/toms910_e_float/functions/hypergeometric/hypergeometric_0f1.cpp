
#include <e_float/e_float.h>
#include <functions/elementary/elementary.h>
#include <functions/gamma/gamma.h>
#include <functions/hypergeometric/hypergeometric.h>

e_float ef::hyperg_0f1(const e_float& b, const e_float& x)
{
  // This implementation is not complete.
  return ef::hyp0F1(b, x);
}

e_float ef::hyperg_0f1_reg(const e_float& b, const e_float& x)
{
  // This implementation is not complete.
  return ef::hyp0F1(b, x) / ef::gamma(b);
}
