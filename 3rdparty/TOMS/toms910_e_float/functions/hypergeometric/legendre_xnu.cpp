
#include <e_float/e_float.h>
#include <functions/constants/constants.h>
#include <functions/elementary/elementary.h>
#include <functions/hypergeometric/legendre_xnu.h>

e_float Legendre_Series::LegendreXnu::MyLegendre(const INT32 n, const e_float& u, const e_float& x) const
{
  if(ef::isint(u))
  {
    const INT32 m = ef::to_int32(u);

    return Legendre_nm(n, m, x);
  }

  if(n < static_cast<INT32>(0))
  {
    return AtReflectNegativeDegree(n, u, x);
  }

  if(ef::isneg(x))
  {
    return AtReflectNegativeArgument(n, u, x);
  }

  if(n < static_cast<INT32>(10))
  {
    return AtOnePlus(n, u, x);
  }

  if(NeedsReflectNegativeOrder(u))
  {
    return AtReflectNegativeOrder(n, u, x);
  }

  // None of the attempted series has converged successfully. Use upward recursion of the
  // degree and upward (rightward) recursion of the order to navigate to the desired value.

  const INT32   m          = ef::to_int32(u);
  const bool    b_recur_mu = (m > static_cast<INT32>(1));
  const e_float um         = (!b_recur_mu ? u : (u - m));

  e_float Pvu_m2 = MyLegendre(static_cast<INT32>(0), um, x);
  e_float Pvu_m1 = MyLegendre(static_cast<INT32>(1), um, x);
  e_float Pvu;

  if(n == static_cast<INT32>(0))
  {
    Pvu = Pvu_m2;
  }
  else if(n == static_cast<INT32>(1))
  {
    Pvu = Pvu_m1;
  }
  else
  {
    for(INT32 i = static_cast<INT32>(2); i <= n; i++)
    {
      // Use upward recursion of the degree.
      const e_float two_v_minus_one    = (i * static_cast<INT32>(2)) - ef::one();
      const e_float u_plus_v_minus_one = (i + um) - ef::one();
      const e_float v_minus_u          =  i - um;

      Pvu = (((two_v_minus_one * x) * Pvu_m1) - (u_plus_v_minus_one * Pvu_m2)) / v_minus_u;

      Pvu_m2 = Pvu_m1;
      Pvu_m1 = Pvu;
    }
  }

  if(b_recur_mu)
  {
    // Use rightward (upward) recursion of the order.
    const e_float x_squared = x * x;

    const e_float sqrt_x_term = ef::sqrt(ef::one() - x_squared);

    e_float ur      = um + ef::two();
    e_float Pvu_um2 = Pvu;
    e_float Pvu_um1 = MyLegendre(n, um + ef::one(), x);
    e_float Pvu_u;

    const INT32 n_times_n_plus_one = static_cast<INT32>(n * static_cast<INT32>(n  + static_cast<INT32>(1)));

    for(INT32 j = static_cast<INT32>(2); j <= m; j++)
    {
      const e_float u_minus_one = ur - ef::one();
      const e_float u_minus_two = ur - ef::two();

      e_float term1 = -(((u_minus_one * static_cast<INT32>(2)) * x) * Pvu_um1) / sqrt_x_term;
      e_float term2 =   ((u_minus_one * u_minus_two) - n_times_n_plus_one) * Pvu_um2;

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
