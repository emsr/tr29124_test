
#include <functions/bessel/bessel.h>
#include <functions/bessel/bessel_convergence_limits.h>
#include <functions/bessel/bessel_recursion.h>
#include <functions/bessel/bessel_recursion_order.h>
#include <functions/gamma/gamma.h>
#include <functions/hypergeometric/hypergeometric.h>
#include <functions/hypergeometric/hypergeometric_util.h>
#include <functions/tables/tables.h>

namespace BesselIv_Series
{
  INT32 OrderAsymptoticLimit(void);

  e_float AtZero         (const e_float& v, const e_float& x);
  e_float AtTransition   (const e_float& v, const e_float& x);
  e_float AtInfinity     (const e_float& v, const e_float& x);
  e_float DebyeAsymptotic(const e_float& v, const e_float& x);
}

namespace BesselIn_Series
{
  static e_float AtZero      (const INT32 n, const e_float& x);
  static e_float AtTransition(const INT32 n, const e_float& x);
  static e_float AtInfinity  (const INT32 n, const e_float& x);
}

namespace BesselJn_Series
{
  INT32 MaximumOrderForBesselJn(void);
}

static e_float BesselIn_Series::AtZero(const INT32 n, const e_float& x)
{
  return BesselIv_Series::AtZero(e_float(n), x);
}

static e_float BesselIn_Series::AtTransition(const INT32 n, const e_float& x)
{
  if(n < BesselIv_Series::OrderAsymptoticLimit())
  {
    const e_float two_over_x = ef::two() / x;

    // Do the downward recursion of In in various forms. If the starting point
    // for downward recursion is below the lower bound of the limit for the
    // uniform asymptotic expansion, then use downward recursion combined with
    // a Neumann sum for normalization. Otherwise use downward recursion starting
    // from two cyl_bessel_i functions with orders lying above the lower limit of
    // convergence for the uniform asymptotic expansion.

    // Downward recursion formula for In:
    //
    //                  In+1
    //   In = [ 2 (n+1) ---- ] + In+2
    //                   x 

    const double xd = ef::to_double(x);

    const INT32 N0 = BesselRecursionOrder::RecursionStartOrderI0(   xd);
    const INT32 Nn = BesselRecursionOrder::RecursionStartOrderIn(n, xd);
    const INT32 N2 = static_cast<INT32>(n + static_cast<INT32>(2));
    const INT32 Nm = static_cast<INT32>(static_cast<INT32>(2) * (std::max(N2, std::max(N0, Nn)) / static_cast<INT32>(2)));

    if(Nm < BesselIv_Series::OrderAsymptoticLimit())
    {
      return BesselRecursion::RecurIn(n, x);
    }
    else
    {
      // Use downward recursion starting from the lower bound of the asymptotic limit.

      e_float In_p2 = ef::cyl_bessel_i(static_cast<INT32>(BesselIv_Series::OrderAsymptoticLimit() + static_cast<INT32>(2)), x);
      e_float In_p1 = ef::cyl_bessel_i(static_cast<INT32>(BesselIv_Series::OrderAsymptoticLimit() + static_cast<INT32>(1)), x);
      e_float In;

      // Do the downward recursion.

      for(INT32 nv = BesselIv_Series::OrderAsymptoticLimit(); nv >= n; nv--)
      {
        const INT32 n_plus_one = static_cast<INT32>(nv + static_cast<INT32>(1));
        
        In    = ((In_p1 * two_over_x) * n_plus_one) + In_p2;
        In_p2 = In_p1;
        In_p1 = In;
      }

      return In;
    }
  }
  else
  {
    return BesselIv_Series::DebyeAsymptotic(e_float(n), x);
  }
}

static e_float BesselIn_Series::AtInfinity(const INT32 n, const e_float& x)
{
  return BesselIv_Series::AtInfinity(e_float(n), x);
}

e_float ef::cyl_bessel_i(const INT32 n, const e_float& x)
{
  if(ef::iszero(x))
  {
    return ((n == static_cast<INT32>(0)) ? ef::one() : ef::zero());
  }

  if(ef::isneg(x))
  {
    // Use In(-x) = (-1)^n I_n(x) for x < 0.
    const e_float bn = cyl_bessel_i(n, x);

    return (((n % static_cast<INT32>(2)) == static_cast<INT32>(0)) ? bn : -bn);
  }

  if(n < static_cast<INT32>(0))
  {
    // Use I_-n = I_n for n < 0.
    return cyl_bessel_i(-n, x);
  }

  if(n > BesselJn_Series::MaximumOrderForBesselJn())
  {
    // No support for order with value |n| > 2*10^9.
    return ef::zero();
  }

  // Check for the possible convergence of the Taylor series.

  const double xd = ef::to_double(x);
  const double vd = static_cast<double>(n);

  const bool bo_converge_small = (xd < 2.0) || (xd < BesselConvergenceLimits::UpperLimitSmallX_Iv(vd));

  if(bo_converge_small)
  {
    return BesselIn_Series::AtZero(n, x);
  }
  else
  {
    // Check for the possible convergence of the asymptotic expansion.
    const bool bo_converge_large = (vd < 30000.0) && (xd > BesselConvergenceLimits::LowerLimitLargeX_Iv(vd));

    return bo_converge_large ? BesselIn_Series::AtInfinity(n, x)
                             : BesselIn_Series::AtTransition(n, x);
  }
}
