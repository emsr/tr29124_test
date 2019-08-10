
#include <e_float/e_float.h>
#include <functions/constants/constants.h>
#include <functions/elementary/elementary.h>
#include <functions/gamma/gamma.h>
#include <functions/hypergeometric/legendre_xnm.h>

e_float Legendre_Series::LegendreXnm::MyLegendre(const INT32 n, const INT32 m, const e_float& x) const
{
  // Use specialized subroutines for the 0th and 1st orders and degrees.
  if(n == static_cast<INT32>(0) && m == static_cast<INT32>(0)) { return L_00(x); }
  if(n == static_cast<INT32>(1) && m == static_cast<INT32>(0)) { return L_10(x); }
  if(n == static_cast<INT32>(0) && m == static_cast<INT32>(1)) { return L_01(x); }
  if(n == static_cast<INT32>(1) && m == static_cast<INT32>(1)) { return L_11(x); }

  if(n < static_cast<INT32>(0))
  {
    return AtReflectNegativeDegree(n, m, x);
  }

  if(AtResultForDegreeAndOrderIsZero(n, m))
  {
    return ef::zero();
  }

  if(m < static_cast<INT32>(0))
  {
    return AtReflectNegativeOrder(n, m, x);
  }

  if(ef::isneg(x))
  {
    return AtReflectNegativeArgument(n, m, x);
  }

  if(ef::iszero(x))
  {
    // Special treatment for argument = 0.
    return AtIdenticallyZero(n, m);
  }

  if(ef::isone(x))
  {
    // Special treatment for argument = +1.
    return AtIdenticallyOne(n, m);
  }

  // Use recursion formulas to calculate Xnm.
  return MyRecursion(n, m, x);
}

e_float Legendre_Series::LegendreXnm::AtReflectNegativeOrder(const INT32 n, const INT32 m, const e_float& x) const
{
  // Handle negative degree with n < -1 using Computation of Special Functions,
  // Zhang & Jin, 4.5.6, page 100.
  const INT32 mm = static_cast<INT32>(-m);

  const bool b_negate = (static_cast<INT32>(m % static_cast<INT32>(2)) != static_cast<INT32>(0));

  const e_float n_minus_m_fact = ef::factorial(static_cast<UINT32>(n - mm));
  const e_float n_plus_m_fact  = ef::factorial(static_cast<UINT32>(n + mm));

  const e_float Lnm = (n_minus_m_fact * MyLegendre(n, mm, x)) / n_plus_m_fact;

  return !b_negate ? Lnm : -Lnm;
}
