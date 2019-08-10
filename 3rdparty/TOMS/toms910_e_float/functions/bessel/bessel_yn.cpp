
#include <functions/bessel/bessel.h>
#include <functions/bessel/bessel_asymp.h>
#include <functions/bessel/bessel_convergence_limits.h>
#include <functions/bessel/bessel_recursion.h>
#include <functions/gamma/gamma.h>
#include <functions/hypergeometric/hypergeometric_util.h>
#include <functions/gamma/polygamma.h>

namespace BesselYv_Series
{
  e_float AtInfinity(const e_float& v, const e_float& x);
}

namespace BesselYn_Series
{
  static e_float AtZero      (const INT32 n, const e_float& x);
  static e_float AtTransition(const INT32 n, const e_float& x);
  static e_float AtInfinity  (const INT32 n, const e_float& x);
}

namespace BesselJn_Series
{
  INT32 MaximumOrderForBesselJn(void);
}

e_float BesselYn_Series::AtZero(const INT32 n, const e_float& x)
{
  // Use the Taylor series expansion for Yn(x), with x small (compared to the order).
  // http://functions.wolfram.com/Bessel-TypeFunctions/BesselY/06/01/04/01/02/

  const e_float x_half         = x / static_cast<INT32>(2);
  const e_float x_half_squared = x_half * x_half;
  const e_float n_fact         = ef::factorial(static_cast<UINT32>(n));

  e_float n_minus_k_minus_one_fact     = n_fact / n;
  e_float x_half_pow_two_k_over_k_fact = ef::one();

  e_float sum0 = n == static_cast<INT32>(0) ? ef::zero() : n_minus_k_minus_one_fact;
  
  for(INT32 k = static_cast<INT32>(1); k < n; k++)
  {
    n_minus_k_minus_one_fact     /= static_cast<INT32>(n - k);
    x_half_pow_two_k_over_k_fact *= x_half_squared;
    x_half_pow_two_k_over_k_fact /= k;
    
    const e_float term0 = n_minus_k_minus_one_fact * x_half_pow_two_k_over_k_fact;
    
    const INT64 order_check0 = static_cast<INT64>(term0.order() - sum0.order());

    if((k > static_cast<INT32>(20)) && (order_check0 < -ef::tol()))
    {
      break;
    }

    sum0 += term0;
  }

  e_float one_over_k_fact_k_plus_n_fact = ef::one() / n_fact;
  
  e_float psi_k_plus_one        = -ef::euler_gamma();
  e_float psi_k_plus_n_plus_one =  ef::poly_gamma(e_float(n + 1));

  e_float x_half_pow_two_k = ef::one();

  e_float sum1 = one_over_k_fact_k_plus_n_fact;
  e_float sum2 = (psi_k_plus_one + psi_k_plus_n_plus_one) * one_over_k_fact_k_plus_n_fact;

  bool b_neg_terms = true;

  for(INT32 k = static_cast<INT32>(1); k < ef::max_iteration(); k++)
  {
    x_half_pow_two_k *= x_half_squared;

    one_over_k_fact_k_plus_n_fact /= k;
    one_over_k_fact_k_plus_n_fact /= (static_cast<UINT32>(k + n));

    psi_k_plus_one        += ef::one() / k;
    psi_k_plus_n_plus_one += ef::one() / static_cast<INT32>(k + n);

    const e_float term1 = x_half_pow_two_k * one_over_k_fact_k_plus_n_fact;
    const e_float term2 = (psi_k_plus_one + psi_k_plus_n_plus_one) * term1;

    const INT64 order_check1 = static_cast<INT64>(term1.order() - sum1.order());
    const INT64 order_check2 = static_cast<INT64>(term2.order() - sum2.order());

    if((k > static_cast<INT32>(20)) && (order_check1 < -ef::tol()) && (order_check2 < -ef::tol()))
    {
      break;
    }

    if(!b_neg_terms)
    {
      sum1 += term1;
      sum2 += term2;
    }
    else
    {
      sum1 -= term1;
      sum2 -= term2;
    }

    b_neg_terms = !b_neg_terms;
  }

  const e_float x_half_pow_n = ef::pown(x_half, static_cast<INT64>(n));

  const e_float part0 = sum0 / x_half_pow_n;
  const e_float part1 = (ef::two() * ef::log(x_half)) * sum1;

  return ((x_half_pow_n * (part1 - sum2)) - part0) / ef::pi();
}

static e_float BesselYn_Series::AtTransition(const INT32 n, const e_float& x)
{
  const double xd = ef::to_double(x);

  if(xd < BesselAsymp::Method::xmin())
  {
    // Compute Y0, Y1 and subsequently use upward recursion if necessary. See Zhang & Jin,
    // Eq. 5.3.8 and Eq. 5.3.9 to expand Y0 and Y1 using a series of Jn.

    // Compute the series of Jn.
    std::deque<e_float> Jn;
    static_cast<void>(BesselRecursion::RecurJn(n, x, &Jn));

    // Compute Y0 and Y1 using a series of Jn.
    e_float sum_j2k        = ef::zero();
    e_float sum_j2k_plus_1 = ef::zero();

    bool is_neg_term = true;    

    for(INT32 k = static_cast<INT32>(1); k < static_cast<INT32>(Jn.size() / 2u); k++)
    {
      const INT32 two_k          = static_cast<INT32>(static_cast<INT32>(2) * k);
      const INT32 two_k_plus_one = static_cast<INT32>(static_cast<INT32>(1) + two_k);

      const e_float term_j2k        =   Jn[static_cast<std::size_t>(two_k)] / k;
      const e_float term_j2k_plus_1 = ((Jn[static_cast<std::size_t>(two_k_plus_one)] * two_k_plus_one) / k) / static_cast<INT32>(k + static_cast<INT32>(1));

      if(!is_neg_term)
      {
        sum_j2k        += term_j2k;
        sum_j2k_plus_1 += term_j2k_plus_1;
      }
      else
      {
        sum_j2k        -= term_j2k;
        sum_j2k_plus_1 -= term_j2k_plus_1;
      }

      is_neg_term = !is_neg_term;
    }

    // Compute and set Y0.
    const e_float x_half = x / static_cast<INT32>(2);
    const e_float lg_x_half_plus_euler_gamma = ef::log(x_half) + ef::euler_gamma();

    e_float Yn = ((lg_x_half_plus_euler_gamma * Jn[0u]) / ef::pi_half()) - (sum_j2k / ef::pi_quarter());

    if(n == static_cast<INT32>(0))
    {
      return Yn;
    }
    else
    {
      // Compute and set Y1.
      e_float Yn_p1 = -(Jn[0u] / (x * ef::pi_half())) + (((lg_x_half_plus_euler_gamma - ef::one()) * Jn[1u]) / ef::pi_half()) - (sum_j2k_plus_1 / ef::pi_half());

      if(n == static_cast<INT32>(1))
      {
        return Yn_p1;
      }
      else
      {
        // Do the upward recursion of Yn:
        //
        //                  Yn+1
        // Yn+2 = [ 2 (n+1) ---- ] - Yn
        //                   x 

        e_float Yn_p2;

        const e_float two_over_x = ef::two() / x;

        for(INT32 nv = static_cast<INT32>(0); nv <= static_cast<INT32>(n - static_cast<INT32>(2)); nv++)
        {
          const INT32 n_plus_one = static_cast<INT32>(nv + static_cast<INT32>(1));

          // Upward recursion for Yn_p2.
          Yn_p2 = ((n_plus_one * Yn_p1) * two_over_x) - Yn;
          Yn    = Yn_p1;
          Yn_p1 = Yn_p2;
        }

        return Yn_p2;
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
      return uniform_asy_z_gt_one.calc(e_float(n), x);
    }
    else if(uniform_asy_z_lt_one.within(n, xd))
    {
      return uniform_asy_z_lt_one.calc(e_float(n), x);
    }

    // The order was not within a convergence band. Use upward recursion
    // starting from the upper limit of the previous lower convergence band.
    // This upper limit of the previous lower convergence band is actually
    // equal to zero for the lower convergence band.

    // Find the upper limit of the previous lower convergence band.
    const INT32 n_hi_z_gt_one = static_cast<INT32>(0);
    const INT32 n_hi_z_lt_one = uniform_asy_z_gt_one.upper(xd);

    const INT32 nr = n > n_hi_z_lt_one ? n_hi_z_lt_one : n_hi_z_gt_one;

    e_float Yn    = ef::cyl_bessel_y(static_cast<INT32>(nr - static_cast<INT32>(2)), x);
    e_float Yn_p1 = ef::cyl_bessel_y(static_cast<INT32>(nr - static_cast<INT32>(1)), x);
    e_float Yn_p2;

    const e_float two_over_x = ef::two() / x;

    // Do the upward recursion.
    for(INT32 nv = nr; nv <= n; nv++)
    {
      Yn_p2 = (((nv - static_cast<INT32>(1)) * Yn_p1) * two_over_x) - Yn;
      Yn    = Yn_p1;
      Yn_p1 = Yn_p2;
    }

    return Yn_p2;
  }
}

e_float BesselYn_Series::AtInfinity(const INT32 n, const e_float& x)
{
  return BesselYv_Series::AtInfinity(e_float(n), x);
}

e_float ef::cyl_bessel_y(const INT32 n, const e_float& x)
{
  if(ef::iszero(x) || ef::isneg(x))
  {
    // Support for zero argument.
    // Support for negative argument.
    return std::numeric_limits<e_float>::quiet_NaN();
  }

  if(n < static_cast<INT32>(0))
  {
    // Use Y_n = (-1)^n Yn for n < 0.
    const INT32 un = static_cast<INT32>(-n);

    const e_float Yn = cyl_bessel_y(un, x);

    return ((un % static_cast<INT32>(2)) == static_cast<INT32>(0)) ? Yn : -Yn;
  }

  if(n > BesselJn_Series::MaximumOrderForBesselJn())
  {
    // No support for order with value |n| > 2*10^9.
    return ef::zero();
  }

  if(x.has_its_own_cyl_bessel_yn())
  {
    return e_float::my_cyl_bessel_yn(n, x);
  }

  const double xd = ef::to_double(x);
  const double vd = static_cast<double>(n);

  // Check for the possible convergence of the Taylor series.

  const bool bo_converge_small = (xd < 2.0) || (xd < BesselConvergenceLimits::UpperLimitSmallX_Yv(vd));

  if(bo_converge_small)
  {
    return BesselYn_Series::AtZero(n, x);
  }
  else
  {
    // Check for the possible convergence of the asymptotic expansion.

    const bool bo_converge_large = (vd < 30000.0) && (xd > BesselConvergenceLimits::LowerLimitLargeX_Yv(vd));

    return bo_converge_large ? BesselYn_Series::AtInfinity(n, x)
                             : BesselYn_Series::AtTransition(n, x);
  }
}
