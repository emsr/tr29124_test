
#include <functions/bessel/bessel.h>
#include <functions/bessel/bessel_asymp.h>
#include <functions/bessel/bessel_convergence_limits.h>
#include <functions/bessel/bessel_recursion.h>
#include <functions/bessel/bessel_recursion_order.h>
#include <functions/gamma/gamma.h>
#include <functions/hypergeometric/hypergeometric.h>
#include <functions/hypergeometric/hypergeometric_util.h>
#include <functions/tables/tables.h>
#include <utility/util_digit_scale.h>

namespace BesselIv_Series
{
  INT32 OrderAsymptoticLimit(void)
  {
    return static_cast<INT32>(static_cast<double>(50000.0) * Util::DigitScale());
  }

  e_float AtZero         (const e_float& v, const e_float& x);
  e_float AtTransition   (const e_float& v, const e_float& x);
  e_float AtInfinity     (const e_float& v, const e_float& x);
  e_float DebyeAsymptotic(const e_float& v, const e_float& x);
}

namespace BesselJn_Series
{
  INT32 MaximumOrderForBesselJn(void);
}

e_float BesselIv_Series::AtZero(const e_float& v, const e_float& x)
{
  // http://functions.wolfram.com/Bessel-TypeFunctions/BesselI/06/01/04/01/01/
  const e_float x_half = x / static_cast<INT32>(2);

  return ef::pow(x_half, v) * ef::hyperg_0f1_reg(v + ef::one(), x_half * x_half);
}

e_float BesselIv_Series::AtTransition(const e_float& v, const e_float& x)
{
  static const e_float v_asymp_limit = e_float(BesselIv_Series::OrderAsymptoticLimit());

  if(v < v_asymp_limit)
  {
    const e_float v_frac     = ef::decimal_part(v);
    const e_float two_over_x = ef::two() / x;

    // Do the downward recursion of Iv in various forms. If the starting point
    // for downward recursion is below the lower bound of the limit for the
    // uniform asymptotic expansion, then use downward recursion combined with
    // a Neumann sum for normalization. Otherwise use downward recursion starting
    // from two cyl_bessel_i functions with orders lying above the lower limit of
    // convergence for the uniform asymptotic expansion.

    // Downward recursion formula for Iv:
    //
    //                  Iv+1
    //   Iv = [ 2 (v+1) ---- ] + Iv+2
    //                   x 

    const double xd = ef::to_double(x);
    const INT32  n  = ef::to_int32(v);

    const INT32 N0 = BesselRecursionOrder::RecursionStartOrderI0(   xd);
    const INT32 Nn = BesselRecursionOrder::RecursionStartOrderIn(n, xd);
    const INT32 N2 = static_cast<INT32>(n + static_cast<INT32>(2));
    const INT32 Nm = static_cast<INT32>(static_cast<INT32>(2) * (std::max(N2, std::max(N0, Nn)) / static_cast<INT32>(2)));

    if(Nm < BesselIv_Series::OrderAsymptoticLimit())
    {
      return BesselRecursion::RecurIv(v, x);
    }
    else
    {
      // Use downward recursion starting from the lower bound of the asymptotic limit.

      e_float Iv_p2 = ef::cyl_bessel_i(v_frac + static_cast<INT32>(BesselIv_Series::OrderAsymptoticLimit() + static_cast<INT32>(2)), x);
      e_float Iv_p1 = ef::cyl_bessel_i(v_frac + static_cast<INT32>(BesselIv_Series::OrderAsymptoticLimit() + static_cast<INT32>(1)), x);
      e_float Iv;

      // Do the downward recursion.

      e_float n_plus_one_plus_v_frac = (v_frac + v_asymp_limit) + + static_cast<INT32>(1);

      for(INT32 nv = BesselIv_Series::OrderAsymptoticLimit(); nv >= n; nv--)
      {
        Iv    = ((Iv_p1 * two_over_x) * n_plus_one_plus_v_frac) + Iv_p2;
        Iv_p2 = Iv_p1;
        Iv_p1 = Iv;

        --n_plus_one_plus_v_frac;
      }

      return Iv;
    }
  }
  else
  {
    return BesselIv_Series::DebyeAsymptotic(v, x);
  }
}

e_float BesselIv_Series::AtInfinity(const e_float& v, const e_float& x)
{
  // Use the large-argument exponential series from Abramowitz and Stegun 9.7.1.

  const e_float mu      = (v * v) * static_cast<INT32>(4);
  const e_float eight_x = x * static_cast<INT32>(8);

  e_float eight_x_pow_k_k_fact = eight_x;
  
  e_float mu_term = mu - ef::one();
  e_float sum     = ef::one() - (mu_term / eight_x_pow_k_k_fact);

  bool b_neg_term = false;

  for(INT32 k = static_cast<INT32>(2); k < ef::max_iteration(); k++)
  {
    eight_x_pow_k_k_fact *= (eight_x * k);

    const INT32 two_k_minus_one = (static_cast<INT32>(2) * k) - static_cast<INT32>(1);

    mu_term *= (mu - (e_float(two_k_minus_one) * two_k_minus_one));

    const e_float term = mu_term / eight_x_pow_k_k_fact;

    const INT64 order_check = static_cast<INT64>(term.order() - sum.order());

    if((k > static_cast<INT32>(20)) && (order_check < -ef::tol()))
    {
      break;
    }

    !b_neg_term ? sum += term : sum -= term;

    b_neg_term = !b_neg_term;
  }

  return (ef::exp(x) * sum) / ef::sqrt(ef::two_pi() * x);
}

e_float BesselIv_Series::DebyeAsymptotic(const e_float& v, const e_float& x)
{
  // Abramowitz and Stegun 9.7.7, page 378.

  const e_float z                 = x / v;
  const e_float z_sq              = z * z;
  const e_float sqrt_one_plus_zsq = ef::sqrt(ef::one() + z_sq);
  const e_float t                 = ef::one() / sqrt_one_plus_zsq;

  const e_float one_over_nu = ef::one() / v;

  e_float sum               = ef::one();
  e_float one_over_nu_pow_k = one_over_nu;

  std::deque<e_float> vt;

  for(INT32 k = static_cast<INT32>(1); k < static_cast<INT32>(Tables::A144618().size()); k++)
  {
    const e_float term = BesselAsymp::DebyeU(k, t, vt) * one_over_nu_pow_k;

    if(term.order() < -ef::tol())
    {
      break;
    }

    sum += term;

    one_over_nu_pow_k *= one_over_nu;
  }
  
  const e_float eta    = sqrt_one_plus_zsq + ef::log(z / (ef::one() + sqrt_one_plus_zsq));
  const e_float factor = (ef::sqrt(one_over_nu / ef::two_pi()) * ef::exp(v * eta)) / ef::sqrt(sqrt_one_plus_zsq);

  return factor * sum;
}

e_float ef::cyl_bessel_i(const e_float& v, const e_float& x)
{
  if(ef::isint(v))
  {
    return cyl_bessel_i(ef::to_int32(v), x);
  }

  if(ef::isneg(x))
  {
    return std::numeric_limits<e_float>::quiet_NaN();
  }

  if(ef::isneg(v))
  {
    const e_float vv       = ef::fabs(v);
    const e_float sin_v_pi = ef::sin(vv * ef::pi());
    const e_float kv       = ef::cyl_bessel_k(vv, x);
    const e_float iv       = ef::cyl_bessel_i(vv, x);

    return iv + ((sin_v_pi * kv) / ef::pi_half());
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

  const bool bo_converge_small = (xd < 2.0) || (xd < BesselConvergenceLimits::UpperLimitSmallX_Iv(vd));

  if(bo_converge_small)
  {
    return BesselIv_Series::AtZero(v, x);
  }
  else
  {
    // Check for the possible convergence of the asymptotic expansion.
    const bool bo_converge_large = (vd < 30000.0) && (xd > BesselConvergenceLimits::LowerLimitLargeX_Iv(vd));

    return bo_converge_large ? BesselIv_Series::AtInfinity(v, x)
                             : BesselIv_Series::AtTransition(v, x);
  }
}
