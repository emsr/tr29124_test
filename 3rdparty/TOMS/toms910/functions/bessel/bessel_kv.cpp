
#include <algorithm>

#include <functions/bessel/bessel.h>
#include <functions/bessel/bessel_asymp.h>
#include <functions/bessel/bessel_convergence_limits.h>
#include <functions/bessel/bessel_recursion_order.h>
#include <functions/gamma/gamma.h>
#include <functions/gamma/gamma_util.h>
#include <functions/hypergeometric/hypergeometric.h>
#include <functions/hypergeometric/hypergeometric_util.h>
#include <functions/tables/tables.h>

namespace BesselKv_Series
{
  static e_float AtZero           (const e_float& v, const e_float& x);
         e_float AtTransition     (const e_float& v, const e_float& x);
  static e_float AtTransitionTemme(const e_float& v, const e_float& x, e_float& Kv_p1);
         e_float AtInfinity       (const e_float& v, const e_float& x);
  static e_float DebyeAsymptotic  (const e_float& v, const e_float& x);
}

namespace BesselJn_Series
{
  INT32 MaximumOrderForBesselJn(void);
}

namespace BesselIv_Series
{
  INT32 OrderAsymptoticLimit(void);
}

e_float BesselKv_Series::AtZero(const e_float& v, const e_float& x)
{
  // http://functions.wolfram.com/BesselAiryStruveFunctions/BesselK/06/01/01/01/

  // Compute gamma(v + 1) and gamma(-v + 1). First compute gamma(+v) and gamma(-v).
  // Then recur these up by one.
  e_float gamma_plus_v;
  e_float gamma_minus_v;

  GammaUtil::GammaOfPlusXMinusX(v, gamma_plus_v, gamma_minus_v);

  const e_float gamma_plus_v_plus_one  = gamma_plus_v  *  v;
  const e_float gamma_minus_v_plus_one = gamma_minus_v * -v;

  const e_float x_half         = x / static_cast<INT32>(2);
  const e_float x_half_squared = x_half * x_half;

  const e_float x_half_pow_v = ef::pow(x_half, v);

  const e_float term1 =  ef::hyperg_0f1(ef::one() - v, x_half_squared) / (x_half_pow_v * gamma_minus_v_plus_one);
  const e_float term2 = (ef::hyperg_0f1(ef::one() + v, x_half_squared) * x_half_pow_v) / gamma_plus_v_plus_one;

  return ef::pi_half() * (ef::csc(ef::pi() * v) * (term1 - term2));
}

e_float BesselKv_Series::AtTransition(const e_float& v, const e_float& x)
{
  static const e_float v_asymp_limit = e_float(BesselIv_Series::OrderAsymptoticLimit());

  if(v < v_asymp_limit)
  {
    // Use the series expansion from Temme in transition region for orders
    // less than the asymptotic limit.
    INT32 nr = static_cast<INT32>(0);
    
    e_float vv = v;

    if(vv > ef::one())
    {
      // Calculate the integer number of recursions which are necessary.
      static const INT64 n32_max = static_cast<INT64>(std::numeric_limits<INT32>::max());

      const INT64 vn = ef::to_int64(vv);

      nr = static_cast<INT32>(std::min(vn, n32_max));

      vv = ef::decimal_part(vv);
    }
    
    // Limit the range of the order to within +/- 1/2.
    // For  1/2 <  v <= 1   scale the order to vv = -(1 - |v|).
    if((vv > ef::half() && vv <=  ef::one()))
    {
      vv = -(ef::one() - vv);
      ++nr;
    }

    e_float Kv_p1;
    e_float Kv;

    if(vv == ef::half() || vv == -ef::half())
    {
      // There is special handling for order 1/2.
      const e_float one_over_x = ef::one() / x;

      Kv    = ef::sqrt(ef::pi_half() * one_over_x) / ef::exp(x);
      Kv_p1 = Kv * (ef::one() + one_over_x);
    }
    else
    {
      // Carry out the series expansion from Temme.
      Kv = BesselKv_Series::AtTransitionTemme(vv, x, Kv_p1);
    }

    // Handle any recursions, if these are necessary.
    if(nr == static_cast<INT32>(0))
    {
      return Kv;
    }
    else if(nr == static_cast<INT32>(1))
    {
      return Kv_p1;
    }
    else
    {
      // Do the recursions.
      e_float Kv_m1;

      const e_float two_over_x = ef::two() / x;

      for(INT32 n = static_cast<INT32>(1); n < nr; n++)
      {
        Kv_m1 = Kv;
        Kv    = Kv_p1;

        ++vv;

        Kv_p1 = ((vv * two_over_x) * Kv) + Kv_m1;
      }

      return Kv_p1;
    }
  }
  else
  {
    // Use the uniform asymptotic expansion for orders greater than the asymptotic limit.
    return BesselKv_Series::DebyeAsymptotic(v, x);
  }
}

static e_float BesselKv_Series::AtTransitionTemme(const e_float& v, const e_float& x, e_float& Kv_p1)
{
  // Set up the recursion in Temme 3.5, 3.6 and the sum in 3.13.
  
  const e_float v_sq  = v * v;
  const e_float two_x = x * static_cast<INT32>(2);

  const e_float half_minus_v_fabs = ef::half() - ef::fabs(v);
  const bool    b_v_near_half     = ef::small_arg(half_minus_v_fabs);

  const INT32 N = (!b_v_near_half ? BesselRecursionOrder::RecursionStartOrderKv            (ef::to_double(v),                 ef::to_double(x))
                                  : BesselRecursionOrder::RecursionStartOrderKv_v_near_half(ef::to_double(half_minus_v_fabs), ef::to_double(x)));

  e_float n_minus_half = N - ef::half();
  e_float an           = (((n_minus_half * n_minus_half) - v_sq) / static_cast<INT32>(N + static_cast<INT32>(1))) / N;
  e_float bn           = (static_cast<INT32>(static_cast<INT32>(2) * N) + two_x) / static_cast<INT32>(N + static_cast<INT32>(1));

  e_float kn_p1;
  e_float kn    = ef::one();
  e_float kn_m1 = bn / an;

  e_float sum_kn = kn + kn_m1;

  for(INT32 n = static_cast<INT32>(N - static_cast<INT32>(1)); n >= static_cast<INT32>(1); n--)
  {
    --n_minus_half;

    an = (((n_minus_half * n_minus_half) - v_sq) / static_cast<INT32>(n + static_cast<INT32>(1))) / n;
    bn = (static_cast<INT32>(static_cast<INT32>(2) * n) + two_x) / static_cast<INT32>(n + 1);

    kn_p1 = kn;
    kn    = kn_m1;

    kn_m1 = ((bn * kn) - kn_p1) / an;
    
    sum_kn += kn_m1;
  }
  
  Kv_p1 = kn;

  const e_float Kv = (ef::sqrt_pi() * ef::exp(-x)) * (kn_m1 / (ef::sqrt(two_x) * sum_kn));
  
  Kv_p1 = Kv * (((v + x) + (ef::half() - (kn / kn_m1))) / x);
  
  return Kv;
}

e_float BesselKv_Series::AtInfinity(const e_float& v, const e_float& x)
{
  const e_float factor = ef::sqrt(ef::pi_half() / x) * ef::exp(-x);

  return factor * ef::hyperg_2f0( v + ef::half(),
                                 -v + ef::half(),
                                  ef::one_minus() / (x * static_cast<INT32>(2)));
}

static e_float BesselKv_Series::DebyeAsymptotic(const e_float& v, const e_float& x)
{
  // Abramowitz and Stegun 9.7.8, page 378.

  const e_float z                 = x / v;
  const e_float z_sq              = z * z;
  const e_float sqrt_one_plus_zsq = ef::sqrt(ef::one() + z_sq);
  const e_float t                 = ef::one() / sqrt_one_plus_zsq;

  const e_float one_over_nu = ef::one() / v;

  e_float sum               = ef::one();
  e_float one_over_nu_pow_k = one_over_nu;
  bool    neg_term          = true;

  std::deque<e_float> vt;

  for(INT32 k = static_cast<INT32>(1); k < static_cast<INT32>(Tables::A144618().size()); k++)
  {
    const e_float term = BesselAsymp::DebyeU(k, t, vt) * one_over_nu_pow_k;

    if(term.order() < -ef::tol())
    {
      break;
    }

    !neg_term ? sum += term : sum -= term;

    neg_term = !neg_term;

    one_over_nu_pow_k *= one_over_nu;
  }
  
  const e_float eta    = sqrt_one_plus_zsq + ef::log(z / (ef::one() + sqrt_one_plus_zsq));
  const e_float factor = (ef::sqrt(ef::pi_half() * one_over_nu) * ef::exp(-v * eta)) / ef::sqrt(sqrt_one_plus_zsq);

  return factor * sum;
}

e_float ef::cyl_bessel_k(const e_float& v, const e_float& x)
{
  if(ef::isneg(x))
  {
    return std::numeric_limits<e_float>::quiet_NaN();
  }

  const e_float vv = ef::fabs(v);

  static const e_float nu_max = e_float(BesselJn_Series::MaximumOrderForBesselJn());

  if(vv > nu_max)
  {
    // No support for order with value |v| > 2*10^9.
    return ef::zero();
  }

  // Check for the possible convergence of the Taylor series.

  const double xd = ef::to_double(x);
  const double vd = ef::to_double(v);

  const bool bo_converge_small = (xd < 2.0) || (xd < BesselConvergenceLimits::UpperLimitSmallX_Kv(vd));

  if(bo_converge_small)
  {
    return BesselKv_Series::AtZero(vv, x);
  }
  else
  {
    // Check for the possible convergence of the asymptotic expansion.

    const bool bo_converge_large = (vd < 30000.0) && (xd > BesselConvergenceLimits::LowerLimitLargeX_Kv(vd));

    return bo_converge_large ? BesselKv_Series::AtInfinity(vv, x)
                             : BesselKv_Series::AtTransition(vv, x);

  }
}
