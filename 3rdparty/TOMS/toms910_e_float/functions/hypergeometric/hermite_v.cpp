
#include <functions/complex/e_float_complex.h>
#include <functions/constants/constants.h>
#include <functions/elementary/elementary.h>
#include <functions/hypergeometric/hermite.h>
#include <functions/hypergeometric/hypergeometric.h>
#include <functions/hypergeometric/parabolic_cylinder.h>
#include <functions/polynomials/polynomials.h>

e_float ef::hermite(const e_float& v, const e_float& x)
{
  if(ef::isint(v) && ef::ispos(v))
  {
    return ef::hermite(static_cast<INT32>(ef::to_int64(v)), x);
  }
  else
  {
    const e_float factor = ef::pow(ef::two(), v / static_cast<INT32>(2)) * ef::exp((x * x) / static_cast<INT32>(2));

    return factor * ef::weber_d(v, ef::sqrt2() * x);
  }
}
