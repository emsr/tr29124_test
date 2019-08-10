
#include <e_float/e_float.h>
#include <functions/constants/constants.h>
#include <functions/elementary/elementary.h>
#include <functions/integer/integer.h>
#include <functions/hypergeometric/legendre_xvu.h>

e_float Legendre_Series::LegendreXvu::MyLegendre(const e_float& v, const e_float& u, const e_float& x) const
{
  if(ef::isone(ef::fabs(x)) || (ef::fabs(x) > ef::one()))
  {
    return std::numeric_limits<e_float>::quiet_NaN();
  }

  if(ef::isint(v))
  {
    const INT32 n = ef::to_int32(v);

    return Legendre_nu(n, u, x);
  }

  if(ef::isint(u))
  {
    const INT32 m = ef::to_int32(u);

    return Legendre_vm(v, m, x);
  }

  if(v < -ef::one())
  {
    return AtReflectNegativeDegree(v, u, x);
  }

  if(ef::isneg(x))
  {
    return AtReflectNegativeArgument(v, u, x);
  }

  if(v < ef::ten())
  {
    // Use the Taylor series near +1 for small values of v.
    return AtOnePlus(v, u, x);
  }

  if(NeedsReflectNegativeOrder(u))
  {
    return AtReflectNegativeOrder(v, u, x);
  }

  // The Taylor series was not used because it was not expected to converge.
  // Use upward recursion of the degree and upward (rightward) recursion of
  // the order to navigate to the desired Legendre function value.

  const INT32   m          = ef::to_int32(u);
  const bool    b_recur_mu = (m > static_cast<INT32>(1));
  const e_float um         = (!b_recur_mu ? u : (u - m));

  // Perform upward recursion of the degree v.
  const e_float v_frac = ef::decimal_part(v);

  e_float Pvu_m2 = MyLegendre(v_frac, um, x);
  e_float Pvu_m1 = MyLegendre(v_frac + static_cast<INT32>(1), um, x);
  e_float Pvu;

  e_float vr = v_frac + static_cast<INT32>(2);

  const INT32 nv = ef::to_int32(v);

  for(INT32 i = static_cast<INT32>(2); i <= nv; i++)
  {
    const e_float two_v_minus_one    = (vr * static_cast<INT32>(2)) - ef::one();
    const e_float u_plus_v_minus_one = (vr + um) - ef::one();
    const e_float v_minus_u          =  vr - um;

    Pvu = (((two_v_minus_one * x) * Pvu_m1) - (u_plus_v_minus_one * Pvu_m2)) / v_minus_u;

    Pvu_m2 = Pvu_m1;
    Pvu_m1 = Pvu;

    ++vr;
  }

  // Perform upward (rightward) recursion of the order u.
  if(b_recur_mu)
  {
    const e_float x_squared = x * x;

    const e_float sqrt_x_term = ef::sqrt(ef::one() - x_squared);

    e_float ur      = um + ef::two();
    e_float Pvu_um2 = Pvu;
    e_float Pvu_um1 = MyLegendre(v, um + ef::one(), x);
    e_float Pvu_u;

    const e_float v_times_v_plus_one  = v * (v  + ef::one());

    for(INT32 j = static_cast<INT32>(2); j <= m; j++)
    {
      const e_float u_minus_one = ur - ef::one();
      const e_float u_minus_two = ur - ef::two();

      e_float term1 = -(((u_minus_one * static_cast<INT32>(2)) * x) * Pvu_um1) / sqrt_x_term;
      e_float term2 =   ((u_minus_one * u_minus_two) - v_times_v_plus_one) * Pvu_um2;

      Pvu_u = term1 + term2;

      Pvu_um2 = Pvu_um1;
      Pvu_um1 = Pvu_u;

      ++ur;
    }

    return Pvu_u;
  }
  else
  {
    return Pvu;
  }
}
