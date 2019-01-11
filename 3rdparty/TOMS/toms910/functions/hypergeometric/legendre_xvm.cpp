
#include <e_float/e_float.h>
#include <functions/constants/constants.h>
#include <functions/elementary/elementary.h>
#include <functions/hypergeometric/hypergeometric_util.h>
#include <functions/hypergeometric/legendre_xvm.h>
#include <functions/integer/integer.h>

e_float Legendre_Series::LegendreXvm::MyLegendre(const e_float& v, const INT32 m, const e_float& x) const
{
  if(ef::isone(ef::fabs(x)) || (ef::fabs(x) > ef::one()))
  {
    return std::numeric_limits<e_float>::quiet_NaN();
  }

  if(ef::isint(v))
  {
    return Legendre_nm(ef::to_int32(v), m, x);
  }

  if(m == static_cast<INT32>(0))
  {
    return Legendre_v(v, x);
  }

  if(m < static_cast<INT32>(0))
  {
    return AtReflectNegativeOrder(v, m, x);
  }

  if(v < -ef::one())
  {
    return AtReflectNegativeDegree(v, m, x);
  }

  if(ef::iszero(x))
  {
    return AtIdenticallyZero(v, m);
  }

  if(ef::isneg(x))
  {
    return AtReflectNegativeArgument(v, m, x);
  }

  if(v < ef::ten())
  {
    return AtOnePlus(v, m, x);
  }

  // The Taylor series was not used because it was not expected to converge.
  // Use upward recursion of the degree and upward (rightward) recursion of
  // the order to navigate to the desired Legendre function value.

  const bool b_needs_recur_of_m = NeedsRecurOfM(x);

  const INT32 mr = (!b_needs_recur_of_m ? m : static_cast<INT32>(1));

  const INT32 nv = ef::to_int32(v);

  const e_float v_frac = ef::decimal_part(v);

  e_float vr     = v_frac + static_cast<INT32>(2);
  e_float Pvu_m2 = MyLegendre(v_frac,             mr, x);
  e_float Pvu_m1 = MyLegendre(v_frac + ef::one(), mr, x);
  e_float Pvu;

  if(nv == static_cast<INT32>(0))
  {
    Pvu = Pvu_m2;
  }
  else if(nv == static_cast<INT32>(1))
  {
    Pvu = Pvu_m1;
  }
  else
  {
    for(INT32 i = static_cast<INT32>(2); i <= nv; i++)
    {
      // Use upward recursion of the degree.
      const e_float two_v_minus_one    = (vr * static_cast<INT32>(2)) - ef::one();
      const e_float u_plus_v_minus_one =  (mr - 1) + vr;
      const e_float v_minus_u          =  vr - mr;

      Pvu = (((two_v_minus_one * x) * Pvu_m1) - (u_plus_v_minus_one * Pvu_m2)) / v_minus_u;

      Pvu_m2 = Pvu_m1;
      Pvu_m1 = Pvu;

      ++vr;
    }
  }

  if(!b_needs_recur_of_m)
  {
    return Pvu;
  }

  const e_float x_squared = x * x;
  const e_float sqrt_x_term = ef::sqrt(ef::one() - x_squared);

  // Use rightward (upward) recursion of the order.
  e_float Pvu_um2 = Legendre_v(v, x);
  e_float Pvu_um1 = Pvu;
  e_float Pvu_u;

  if(m == static_cast<INT32>(1))
  {
    return Pvu;
  }
  else
  {
    const e_float v_times_v_plus_one  = v * (v  + ef::one());

    for(INT32 j = static_cast<INT32>(2); j <= m; j++)
    {
      const INT32 m_minus_one = static_cast<INT32>(j - static_cast<INT32>(1));
      const INT32 m_minus_two = static_cast<INT32>(j - static_cast<INT32>(2));

      e_float term1 = -((static_cast<INT32>(m_minus_one * static_cast<INT32>(2)) * x) * Pvu_um1) / sqrt_x_term;
      e_float term2 =   (static_cast<INT32>(m_minus_one * m_minus_two) - v_times_v_plus_one) * Pvu_um2;

      Pvu_u = term1 + term2;

      Pvu_um2 = Pvu_um1;
      Pvu_um1 = Pvu_u;
    }

    return Pvu_u;
  }
}
