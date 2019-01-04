
#include <e_float/e_float.h>
#include <functions/constants/constants.h>
#include <functions/elementary/elementary.h>
#include <functions/gamma/gamma_util.h>
#include <functions/gamma/gamma.h>
#include <functions/gamma/polygamma.h>

void GammaUtil::GammaOfPlusXMinusX(const e_float& x, e_float& gamma_plus_x, e_float& gamma_minus_x)
{
  // Calculate gamma(x) and gamma(-x) without regard for the actual sign of x.
  // In other words, calculate gamma(abs(x)) and set the results appropriately
  // using the reflection formula.
  const bool    bo_x_isneg    =  ef::isneg(x);
  const e_float abs_x         =  ef::fabs(x);
  const e_float gamma_x_pos   =  ef::gamma(abs_x);
  const e_float gamma_x_neg   = -ef::pi() / ((abs_x * gamma_x_pos) * ef::sin(ef::pi() * abs_x));

  gamma_plus_x  =  bo_x_isneg ? gamma_x_neg : gamma_x_pos;
  gamma_minus_x =  bo_x_isneg ? gamma_x_pos : gamma_x_neg;
}

void GammaUtil::DiGammaOfPlusXMinusX(const e_float& x, e_float& psi_plus_x,   e_float& psi_minus_x)
{
  // Calculate di_gamma(x) and di_gamma(-x) without regard for the actual sign of x.
  // In other words, calculate di_gamma(abs(x)) and set the results appropriately
  // using the reflection formula.
  const bool    bo_x_isneg =  ef::isneg(x);
  const e_float abs_x      =  ef::fabs(x);
  const e_float psi_x_pos  =  ef::poly_gamma(abs_x);
  const e_float psi_x_neg  = (psi_x_pos + (ef::one() / abs_x)) + (ef::pi() * ef::cot(ef::pi() * abs_x));

  psi_plus_x  =  bo_x_isneg ? psi_x_neg : psi_x_pos;
  psi_minus_x =  bo_x_isneg ? psi_x_pos : psi_x_neg;
}
