
#include <numeric>

#include <functions/bessel/bessel.h>
#include <functions/bessel/bessel_asymp.h>
#include <functions/bessel/bessel_convergence_limits.h>
#include <functions/bessel/bessel_recursion.h>
#include <functions/hypergeometric/hypergeometric_util.h>

namespace BesselJv_Series
{
  e_float AtZero    (const e_float& v, const e_float& x);
  e_float AtInfinity(const e_float& v, const e_float& x);
}

namespace BesselJn_Series
{
  INT32 MaximumOrderForBesselJn(void);

  static e_float AtZero      (const INT32 n, const e_float& x);
  static e_float AtTransition(const INT32 n, const e_float& x);
  static e_float AtInfinity  (const INT32 n, const e_float& x);
}

INT32 BesselJn_Series::MaximumOrderForBesselJn(void)
{
  return static_cast<INT32>(2000000000L);
}

static e_float BesselJn_Series::AtZero(const INT32 n, const e_float& x)
{
  return BesselJv_Series::AtZero(e_float(n), x);
}

static e_float BesselJn_Series::AtTransition(const INT32 n, const e_float& x)
{
  const double xd = ef::to_double(x);

  if(xd < BesselAsymp::Method::xmin())
  {
    return BesselRecursion::RecurJn(n, x);
  }
  else
  {
    // Analyze the order and the argument in order to select the appropriate asymptotic method
    // and a possible recursion scheme.
    static const BesselJ_AsympMethod::UniformAsy_z_gt_one uniform_asy_z_gt_one;
    static const BesselJ_AsympMethod::UniformAsy_z_lt_one uniform_asy_z_lt_one;

    if(uniform_asy_z_gt_one.within(n, xd))
    {
      return uniform_asy_z_gt_one.calc(e_float(n), x);
    }
    else if(uniform_asy_z_lt_one.within(n, xd))
    {
      return uniform_asy_z_lt_one.calc(e_float(n), x);
    }

    // The order was not within a convergence band. Use downward recursion
    // starting from the lower limit of the next higher convergence band.

    // Find the lower limit of the next higher convergence band.
    const INT32 n_lo_z_gt_one = uniform_asy_z_gt_one.lower(xd);
    const INT32 n_lo_z_lt_one = uniform_asy_z_lt_one.lower(xd);

    const INT32 nr = n < n_lo_z_gt_one ? n_lo_z_gt_one : n_lo_z_lt_one;

    e_float Jm_p2 = ef::cyl_bessel_j(static_cast<INT32>(nr + static_cast<INT32>(2)), x);
    e_float Jm_p1 = ef::cyl_bessel_j(static_cast<INT32>(nr + static_cast<INT32>(1)), x);
    e_float Jm;

    const e_float two_over_x = ef::two() / x;

    // Do the downward recursion.

    for(INT32 m = nr; m >= n; m--)
    {
      Jm    = ((Jm_p1 * two_over_x) * static_cast<INT32>(m + static_cast<INT32>(1))) - Jm_p2;
      Jm_p2 = Jm_p1;
      Jm_p1 = Jm;
    }

    return Jm;
  }
}

static e_float BesselJn_Series::AtInfinity(const INT32 n, const e_float& x)
{
  return BesselJv_Series::AtInfinity(e_float(n), x);
}

e_float ef::cyl_bessel_j(const INT32 n, const e_float& x)
{
  if(ef::iszero(x))
  {
    return n == static_cast<INT32>(0) ? ef::one() : ef::zero();
  }

  if(ef::isneg(x))
  {
    // Use J_n(-x) = (-1)^n J_n(x) for x < 0.
    const e_float Jn = cyl_bessel_j(n, ef::fabs(x));

    return ((n % static_cast<INT32>(2)) == static_cast<INT32>(0)) ? Jn : -Jn;
  }

  if(n < static_cast<INT32>(0))
  {
    // Use J_-n = (-1)^n J_n for n < 0.
    const INT32 un = static_cast<INT32>(-n);

    const e_float Jn = cyl_bessel_j(un, x);

    return ((un % static_cast<INT32>(2)) == static_cast<INT32>(0)) ? Jn : -Jn;
  }

  if(n > BesselJn_Series::MaximumOrderForBesselJn())
  {
    // No support for order with value |n| > 2*10^9.
    return ef::zero();
  }

  if(x.has_its_own_cyl_bessel_jn())
  {
    return e_float::my_cyl_bessel_jn(n, x);
  }

  const double xd = ef::to_double(x);
  const double vd = static_cast<double>(n);

  // Check for the possible convergence of the Taylor series.

  const bool bo_converge_small = (xd <2.0) || (xd < BesselConvergenceLimits::UpperLimitSmallX_Jv(vd));

  if(bo_converge_small)
  {
    return BesselJn_Series::AtZero(n, x);
  }
  else
  {

    // Check for the possible convergence of the asymptotic expansion.

    const bool bo_converge_large = (vd < 30000.0) && (xd > BesselConvergenceLimits::LowerLimitLargeX_Jv(vd));

    return bo_converge_large ? BesselJn_Series::AtInfinity(n, x)
                             : BesselJn_Series::AtTransition(n, x);
  }
}
