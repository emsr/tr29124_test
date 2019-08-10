
#include <e_float/e_float.h>
#include <functions/constants/constants.h>
#include <functions/elementary/elementary.h>
#include <functions/hypergeometric/legendre_xv.h>
#include <functions/integer/integer.h>

e_float Legendre_Series::LegendreXv::MyLegendre(const e_float& v, const e_float& x) const
{
  if(ef::isint(v))
  {
    const INT32 n = ef::to_int32(v);

    return Legendre_n(n, x);
  }

  if(v < ef::one_minus())
  {
    return AtReflectNegativeDegree(v, x);
  }

  if(ef::iszero(x))
  {
    return AtIdenticallyZero(v);
  }

  if(v < ef::twenty())
  {
    // Use series expansion for moderate values of degree.
    if(ef::ispos(x))
    {
      return AtOnePlus(v, x);
    }
    else
    {
      return AtOneMinus(v, x);
    }
  }

  // Use forward recursion to calculate Legendre Functions of higher degree.
  const e_float v_frac = ef::decimal_part(v);

  e_float vr     = v_frac + static_cast<INT32>(2);
  e_float Pv_m2 = MyLegendre(v_frac,             x);
  e_float Pv_m1 = MyLegendre(v_frac + ef::one(), x);
  e_float Pv;

  const INT32 nv = ef::to_int32(v);

  if(nv == static_cast<INT32>(0))
  {
    Pv = Pv_m2;
  }
  else if(nv == static_cast<INT32>(1))
  {
    Pv = Pv_m1;
  }
  else
  {
    for(INT32 i = static_cast<INT32>(2); i <= nv; i++)
    {
      // Use upward recursion of the degree.
      const e_float two_v_minus_one = (vr * static_cast<INT32>(2)) - ef::one();
      const e_float v_minus_one     = vr - ef::one();

      Pv = (((two_v_minus_one * x) * Pv_m1) - (v_minus_one * Pv_m2)) / vr;

      Pv_m2 = Pv_m1;
      Pv_m1 = Pv;

      ++vr;
    }
  }

  return Pv;  
}

