
#include <functions/bessel/bessel.h>
#include <functions/bessel/bessel_asymp.h>
#include <functions/bessel/bessel_convergence_limits.h>
#include <functions/bessel/bessel_recursion.h>
#include <functions/hypergeometric/hypergeometric.h>
#include <functions/hypergeometric/hypergeometric_util.h>
#include <functions/gamma/gamma.h>

namespace BesselJv_Series
{
         e_float AtZero      (const e_float& v, const e_float& x);
  static e_float AtTransition(const e_float& v, const e_float& x);
         e_float AtInfinity  (const e_float& v, const e_float& x);
}

namespace BesselJn_Series
{
  INT32 MaximumOrderForBesselJn(void);
}

e_float BesselJv_Series::AtZero(const e_float& v, const e_float& x)
{
  // Use the Taylor series expansion for x small (compared to the order).
  // http://functions.wolfram.com/BesselAiryStruveFunctions/BesselJ/06/01/01/

  const e_float x_half = x / static_cast<INT32>(2);

  return ef::pow(x_half, v) * ef::hyperg_0f1_reg(v + ef::one(), -(x_half * x_half));
}

static e_float BesselJv_Series::AtTransition(const e_float& v, const e_float& x)
{
  const double xd = ef::to_double(x);
  const INT32  n  = ef::to_int32(v);

  if(xd < BesselAsymp::Method::xmin())
  {
    return BesselRecursion::RecurJv(v, x);
  }
  else
  {
    // Analyze the order and the argument in order to select the appropriate asymptotic method
    // and a possible recursion scheme.
    static const BesselJ_AsympMethod::UniformAsy_z_gt_one uniform_asy_z_gt_one;
    static const BesselJ_AsympMethod::UniformAsy_z_lt_one uniform_asy_z_lt_one;

    if(uniform_asy_z_gt_one.within(n, xd))
    {
      return uniform_asy_z_gt_one.calc(v, x);
    }
    else if(uniform_asy_z_lt_one.within(n, xd))
    {
      return uniform_asy_z_lt_one.calc(v, x);
    }

    // The order was not within a convergence band. Use downward recursion
    // starting from the lower limit of the next higher convergence band.

    // Find the lower limit of the next higher convergence band.
    const INT32 n_lo_z_gt_one = uniform_asy_z_gt_one.lower(xd);
    const INT32 n_lo_z_lt_one = uniform_asy_z_lt_one.lower(xd);

    const INT32 nr = n < n_lo_z_gt_one ? n_lo_z_gt_one : n_lo_z_lt_one;

    const e_float v_frac   = ef::decimal_part(v);
          e_float n_v_frac = (v_frac + nr) + ef::one();

    e_float Jv_p2 = ef::cyl_bessel_j(v_frac + static_cast<INT32>(nr + static_cast<INT32>(2)), x);
    e_float Jv_p1 = ef::cyl_bessel_j(v_frac + static_cast<INT32>(nr + static_cast<INT32>(1)), x);
    e_float Jv;

    const e_float two_over_x = ef::two() / x;

    // Do the downward recursion.
    for(INT32 nv = nr; nv >= n; nv--)
    {
      Jv    = ((Jv_p1 * two_over_x) * n_v_frac) - Jv_p2;
      Jv_p2 = Jv_p1;
      Jv_p1 = Jv;
      
      --n_v_frac;
    }

    return Jv;
  }
}

e_float BesselJv_Series::AtInfinity(const e_float& v, const e_float& x)
{
  // Use the hypergeometric representation in trigonometric form using
  // Hypergeometric4F1 for large positive values of x.
  // http://functions.wolfram.com/BesselAiryStruveFunctions/BesselJ/06/02/02/

  const e_float two_nu                    = e_float(v * static_cast<INT32>(2));
  const e_float two_nu_plus_one_over_four = (two_nu + ef::one()) / static_cast<INT32>(4);
  const e_float trig_arg                  = x - ((ef::pi() * two_nu_plus_one_over_four));
  
  const std::tr1::array<e_float, 4u> p1_data =
  {{
    (ef::one()   - two_nu) / static_cast<INT32>(4),
    (ef::three() - two_nu) / static_cast<INT32>(4),
    two_nu_plus_one_over_four,
    (ef::three() + two_nu) / static_cast<INT32>(4)
  }};

         const std::deque<e_float> p1(p1_data.begin(), p1_data.end());
  static const std::deque<e_float> q1(static_cast<std::size_t>(1u), ef::half());

  const std::tr1::array<e_float, 4u> p2_data =
  {{
    p1_data[1u],
    (ef::five() - two_nu) / static_cast<INT32>(4),
    p1_data[3u],
    (ef::five() + two_nu) / static_cast<INT32>(4)
  }};

         const std::deque<e_float> p2(p2_data.begin(), p2_data.end());
  static const std::deque<e_float> q2(static_cast<std::size_t>(1u), ef::three_half());

  e_float sin_term;
  e_float cos_term;

  ef::sincos(trig_arg, &sin_term, &cos_term);

  const e_float minus_one_over_xsq = ef::one_minus() / (x * x);

  const e_float sum1 = cos_term * ef::hyperg_pfq(p1, q1, minus_one_over_xsq);
  const e_float sum2 = sin_term * ef::hyperg_pfq(p2, q2, minus_one_over_xsq);

  const e_float sum = sum1 + (sum2 * ((ef::one() - (two_nu * two_nu)) / (x * static_cast<INT32>(8))));

  return sum * ef::sqrt(ef::two() / (ef::pi() * x));
}

e_float ef::cyl_bessel_j(const e_float& v, const e_float& x)
{
  if(ef::isint(v))
  {
    return ef::cyl_bessel_j(ef::to_int32(v), x);
  }

  if(ef::iszero(x))
  {
    // Support for zero argument.
    return ef::zero();
  }

  if(ef::isneg(x))
  {
    // Support for negative argument.
    return std::numeric_limits<e_float>::quiet_NaN();
  }

  if(ef::isneg(v))
  {
    // Support for negative order.
    const e_float vv   = ef::fabs(v);
    const e_float pi_v = ef::pi() * vv;

    e_float sin_pi_v;
    e_float cos_pi_v;

    ef::sincos(pi_v, &sin_pi_v, &cos_pi_v);
    
    return (ef::cyl_bessel_j(vv, x) * cos_pi_v) - (ef::cyl_bessel_y(vv, x) * sin_pi_v);
  }

  static const e_float nu_max = e_float(BesselJn_Series::MaximumOrderForBesselJn());

  if(v > nu_max)
  {
    // No support for order with value |v| > 2*10^9.
    return ef::zero();
  }

  // Check for the possible convergence of the Taylor series.

  const double xd = ef::to_double(x);
  const double vd = ef::to_double(v);

  const bool bo_converge_small = (xd < 2.0) || (xd < BesselConvergenceLimits::UpperLimitSmallX_Jv(vd));

  if(bo_converge_small)
  {
    return BesselJv_Series::AtZero(v, x);
  }
  else
  {
    // Check for the possible convergence of the asymptotic expansion.
    const bool bo_converge_large = (vd < 30000.0) && (xd > BesselConvergenceLimits::LowerLimitLargeX_Jv(vd));

    return bo_converge_large ? BesselJv_Series::AtInfinity(v, x)
                             : BesselJv_Series::AtTransition(v, x);
  }
}
