
#include <numeric>

#include <functions/bessel/bessel.h>
#include <functions/bessel/bessel_asymp.h>
#include <functions/bessel/bessel_convergence_limits.h>
#include <functions/bessel/bessel_recursion.h>
#include <functions/hypergeometric/hypergeometric.h>
#include <functions/hypergeometric/hypergeometric_util.h>
#include <functions/gamma/gamma.h>
#include <functions/gamma/gamma_util.h>

namespace BesselYv_Series
{
  static e_float AtZero      (const e_float& v, const e_float& x);
  static e_float AtTransition(const e_float& v, const e_float& x);
         e_float AtInfinity  (const e_float& v, const e_float& x);
}

namespace BesselJn_Series
{
  INT32 MaximumOrderForBesselJn(void);
}

static e_float BesselYv_Series::AtZero(const e_float& v, const e_float& x)
{
  // Use the Taylor series expansion for Yv(x), with x small (compared to the order).
  // http://functions.wolfram.com/Bessel-TypeFunctions/BesselY/06/01/04/01/01/

  // Compute gamma(v + 1) and gamma(-v + 1). First compute gamma(+v) and gamma(-v).
  // Then recur these up by one.
  e_float gamma_plus_v;
  e_float gamma_minus_v;
  
  GammaUtil::GammaOfPlusXMinusX(v, gamma_plus_v, gamma_minus_v);

  const e_float gamma_plus_v_plus_one  = gamma_plus_v  *  v;
  const e_float gamma_minus_v_plus_one = gamma_minus_v * -v;

  const e_float x_half = x / static_cast<INT32>(2);

  const e_float x_half_pow_v = ef::pow(x_half, v);

  const e_float pi_v = ef::pi() * v;
  
  e_float sin_pi_v;
  e_float cos_pi_v;
  ef::sincos(pi_v, &sin_pi_v, &cos_pi_v);
  
  const e_float x_half_squared = x_half * x_half;

  const e_float term1 = ef::hyperg_0f1( v + ef::one(), -x_half_squared) * ((cos_pi_v * x_half_pow_v) / gamma_plus_v_plus_one);
  const e_float term2 = ef::hyperg_0f1(-v + ef::one(), -x_half_squared) * (ef::one() / (x_half_pow_v * gamma_minus_v_plus_one));

  return (term1 - term2) / sin_pi_v;
}

static e_float BesselYv_Series::AtTransition(const e_float& v, const e_float& x)
{
  const double xd = ef::to_double(x);
  const INT32  n  = ef::to_int32(v);

  if(xd < BesselAsymp::Method::xmin())
  {
    // Use the method as described in "Computation with Recurrence Relations",
    // J. Wimp, Chapter 5.7, Eq. 5.48 - Eq. 5.53, page 79.

    std::deque<e_float> Jv;
    static_cast<void>(BesselRecursion::RecurJv(v, x, &Jv));

    // Prepare the recursion of the symbol 'dn' as shown in Eq. 5.49 - Eq. 5.51.
    const e_float two_over_x           = ef::two() / x;
    const e_float v_frac               = ef::decimal_part(v);
    const e_float two_over_x_pow_two_v = ef::pow(two_over_x, v_frac * 2);
    const e_float gamma_one_plus_v     = ef::gamma(ef::one() + v_frac);
    const e_float gamma_sq_one_plus_v  = gamma_one_plus_v * gamma_one_plus_v;
    const e_float pi_v                 = ef::pi() * v_frac;

    const e_float d0 = ef::cot(pi_v) - ((two_over_x_pow_two_v * gamma_sq_one_plus_v) / pi_v);

    e_float d1 = ((two_over_x_pow_two_v * ((v_frac + ef::two()) / (ef::one() - v_frac))) * gamma_sq_one_plus_v) / ef::pi_half();

    e_float Yv =   (d0 * Jv[static_cast<std::size_t>(0u)])
                 + (d1 * Jv[static_cast<std::size_t>(2u)]);

    // Do the recursion of the symbol 'dn' and the summation of the series for Yv.
    const INT32 maxn = static_cast<INT32>((Jv.size() / static_cast<std::size_t>(2u)) - static_cast<std::size_t>(2u));

    for(INT32 nd = static_cast<INT32>(1); nd < maxn; nd++)
    {
      const INT32   n_plus_one   = static_cast<INT32>(nd + static_cast<INT32>(1));
      const e_float two_n_plus_v = v_frac + static_cast<INT32>(static_cast<INT32>(2) * nd);
      const e_float dn_top       = (((two_n_plus_v + static_cast<INT32>(2)) * ((v_frac * static_cast<INT32>(2)) + nd)) * (v_frac + nd)) * d1;
      const e_float dn_bot       = ((-v_frac + n_plus_one) * two_n_plus_v) * n_plus_one;

      const e_float dn = -(dn_top / dn_bot);

      const e_float term = dn * Jv[static_cast<std::size_t>(static_cast<INT32>(2) * n_plus_one)];

      const INT64 order_check = static_cast<INT64>(term.order() - Yv.order());

      if((n > static_cast<INT32>(20)) && (order_check < -ef::tol()))
      {
        break;
      }

      Yv += term;

      d1 = dn;
    }

    if(n == static_cast<INT32>(0))
    {
      return Yv;
    }
    else
    {
      // Use the Wronskian to compute Yv+1.
      e_float Yv_p1 = ((Yv * Jv[static_cast<std::size_t>(1u)]) - (ef::one() / (x * ef::pi_half()))) / Jv[static_cast<std::size_t>(0u)];

      if(n == static_cast<INT32>(1))
      {
        return Yv_p1;
      }
      else
      {
        // Do the upward recursion of Yv:
        //
        //                  Yv+1
        // Yv+2 = [ 2 (v+1) ---- ] - Yv
        //                   x 

        e_float Yv_p2;

        for(INT32 nv = static_cast<INT32>(1); nv < n; nv++)
        {
          const e_float n_plus_v_frac_plus_one = v_frac + nv;

          // Upward recursion for Yv_frac_p2.
          Yv_p2 = ((n_plus_v_frac_plus_one * Yv_p1) * two_over_x) - Yv;
          Yv    = Yv_p1;
          Yv_p1 = Yv_p2;
        }

        return Yv_p2;
      }
    }
  }
  else
  {
    // Analyze the order and the argument in order to select the appropriate asymptotic method
    // and a possible recursion scheme.
    static const BesselY_AsympMethod::UniformAsy_z_gt_one uniform_asy_z_gt_one;
    static const BesselY_AsympMethod::UniformAsy_z_lt_one uniform_asy_z_lt_one;

    if(uniform_asy_z_gt_one.within(n, xd))
    {
      return uniform_asy_z_gt_one.calc(v, x);
    }
    else if(uniform_asy_z_lt_one.within(n, xd))
    {
      return uniform_asy_z_lt_one.calc(v, x);
    }

    // The order was not within a convergence band. Use upward recursion
    // starting from the upper limit of the previous lower convergence band.
    // This upper limit of the previous lower convergence band is actually
    // equal to zero for the lower convergence band.

    // Find the upper limit of the previous lower convergence band.
    const INT32 n_hi_z_gt_one = static_cast<INT32>(0);
    const INT32 n_hi_z_lt_one = uniform_asy_z_gt_one.upper(xd);

    const INT32 nr = n > n_hi_z_lt_one ? n_hi_z_lt_one : n_hi_z_gt_one;

    const e_float v_frac   = ef::decimal_part(v);
          e_float n_v_frac = v_frac + nr;

    e_float Yv    = ef::cyl_bessel_y(v_frac + static_cast<INT32>(nr - static_cast<INT32>(2)), x);
    e_float Yv_p1 = ef::cyl_bessel_y(v_frac + static_cast<INT32>(nr - static_cast<INT32>(1)), x);
    e_float Yv_p2;

    const e_float two_over_x = ef::two() / x;

    // Do the upward recursion.
    for(INT32 nv = nr; nv <= n; nv++)
    {
      Yv_p2 = (((n_v_frac - static_cast<INT32>(1)) * Yv_p1) * two_over_x) - Yv;
      Yv    = Yv_p1;
      Yv_p1 = Yv_p2;

      ++n_v_frac;
    }

    return Yv_p2;
  }
}

e_float BesselYv_Series::AtInfinity(const e_float& v, const e_float& x)
{
  // Use the hypergeometric representation in trigonometric form using
  // Hypergeometric4F1 for large positive values of x.
  // http://functions.wolfram.com/Bessel-TypeFunctions/BesselY/06/02/01/01/02/

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

  const e_float sum1 = sin_term * ef::hyperg_pfq(p1, q1, minus_one_over_xsq);
  const e_float sum2 = cos_term * ef::hyperg_pfq(p2, q2, minus_one_over_xsq);

  const e_float sum = sum1 - (sum2 * ((ef::one() - (two_nu * two_nu)) / (x * static_cast<INT32>(8))));

  return sum * ef::sqrt(ef::two() / (ef::pi() * x));
}

e_float ef::cyl_bessel_y(const e_float& v, const e_float& x)
{
  if(ef::isint(v))
  {
    return ef::cyl_bessel_y(ef::to_int32(v), x);
  }

  if(ef::iszero(x) || ef::isneg(x))
  {
    // Support for zero argument.
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
    
    return (ef::cyl_bessel_j(vv, x) * sin_pi_v) + (ef::cyl_bessel_y(vv, x) * cos_pi_v);
  }

  static const e_float nu_max = e_float(BesselJn_Series::MaximumOrderForBesselJn());

  if(v > nu_max)
  {
    // No support for order with value |v| > 2*10^9.
    return ef::zero();
  }

  // Check for possible convergence of the Taylor series.

  const double xd = ef::to_double(x);
  const double vd = ef::to_double(v);

  const bool bo_converge_small = (xd < 2.0) || (xd < BesselConvergenceLimits::UpperLimitSmallX_Yv(vd));

  if(bo_converge_small)
  {
    return BesselYv_Series::AtZero(v, x);
  }
  else
  {
    // Check for the possible convergence of the asymptotic expansion.

    const bool bo_converge_large = (vd < 30000.0) && (xd > BesselConvergenceLimits::LowerLimitLargeX_Yv(vd));

    return bo_converge_large ? BesselYv_Series::AtInfinity(v, x)
                             : BesselYv_Series::AtTransition(v, x);
  }
}
