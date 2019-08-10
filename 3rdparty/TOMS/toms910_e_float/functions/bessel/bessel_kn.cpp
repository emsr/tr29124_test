
#include <algorithm>

#include <functions/bessel/bessel.h>
#include <functions/bessel/bessel_convergence_limits.h>
#include <functions/bessel/bessel_recursion.h>
#include <functions/hypergeometric/hypergeometric.h>
#include <functions/hypergeometric/hypergeometric_util.h>
#include <functions/tables/tables.h>

namespace BesselKv_Series
{
  e_float AtTransition(const e_float& v, const e_float& x);
  e_float AtInfinity  (const e_float& v, const e_float& x);
}

namespace BesselK0_Series
{
  static e_float AtZero(const e_float& x);
}

namespace BesselK1_Series
{
  static e_float AtZero(const e_float& x);
}

namespace BesselJn_Series
{
  INT32 MaximumOrderForBesselJn(void);
}

static e_float BesselK0_Series::AtZero(const e_float& x)
{
  // Use the Taylor series expansion for K0(x) with argument x <= 2.
  // See Computation of Special Functions, Zhang & Jin, 6.1.11, page 204.
  // See also http://functions.wolfram.com/BesselAiryStruveFunctions/BesselI/06/01/01/
  // for the expansion of I0(x).

  const e_float x_half    = x / static_cast<INT32>(2);
  const e_float x_half_sq = x_half * x_half;

  e_float one_over_k_fact = ef::one();
  e_float sum_one_over_k  = ef::one();
  e_float x_half_sq_pow_k = x_half_sq;

  e_float sum = x_half_sq_pow_k;
  
  for(INT32 k = static_cast<INT32>(2); k < ef::max_iteration(); k++)
  {
    x_half_sq_pow_k *= x_half_sq;
    one_over_k_fact /= k;
    sum_one_over_k  += ef::one() / k;

    const e_float term = sum_one_over_k * (x_half_sq_pow_k * (one_over_k_fact * one_over_k_fact));
    
    const INT64 order_check = static_cast<INT64>(term.order() - sum.order());

    if((k > static_cast<INT32>(20)) && (order_check < -ef::tol()))
    {
      break;
    }
    
    sum += term;
  }

  const e_float I0 = ef::hyperg_0f1_reg(ef::one(), x_half_sq);

  return -((ef::log(x_half) + ef::euler_gamma()) * I0) + sum;
}

static e_float BesselK1_Series::AtZero(const e_float& x)
{
  // Use the Taylor series expansion for K1(x) with argument x <= 2.
  // See Computation of Special Functions, Zhang & Jin, 6.1.12, page 204.
  // See also http://functions.wolfram.com/BesselAiryStruveFunctions/BesselI/06/01/01/
  // for the expansion of I1(x).

  const e_float x_half    = x / static_cast<INT32>(2);
  const e_float x_half_sq = x_half * x_half;

  e_float one_over_k_fact = ef::one();
  e_float sum_one_over_k  = ef::zero();
  e_float x_half_sq_pow_k = ef::one();

  e_float sum = ef::half();
  
  for(INT32 k = static_cast<INT32>(1); k < ef::max_iteration(); k++)
  {
    x_half_sq_pow_k *= x_half_sq;
    one_over_k_fact /= k;
    sum_one_over_k  += ef::one() / k;

    const INT32 k_plus_one = static_cast<INT32>(k + static_cast<INT32>(1));

    const e_float term = (sum_one_over_k + (ef::half() / k_plus_one)) * (x_half_sq_pow_k * (one_over_k_fact * (one_over_k_fact / k_plus_one)));
    
    const INT64 order_check = static_cast<INT64>(term.order() - sum.order());

    if((k > static_cast<INT32>(20)) && (order_check < -ef::tol()))
    {
      break;
    }
    
    sum += term;
  }

  const e_float I1 = x_half * ef::hyperg_0f1_reg(ef::two(), x_half_sq);

  return ((ef::one() / x) + ((ef::log(x_half) + ef::euler_gamma()) * I1)) - (x_half * sum);
}

e_float ef::cyl_bessel_k(const INT32 n, const e_float& x)
{
  if(ef::isneg(x))
  {
    return std::numeric_limits<e_float>::quiet_NaN();
  }

  const INT32 nn = n < 0 ? -n : n;

  if(nn > BesselJn_Series::MaximumOrderForBesselJn())
  {
    // No support for order with value |n| > 2*10^9.
    return ef::zero();
  }

  if(x <= ef::two())
  {
    // Compute K0, (and K1 if necessary) and use recursion if necessary.

    if(nn == static_cast<INT32>(0))
    {
      return BesselK0_Series::AtZero(x);
    }
    else if(nn == static_cast<INT32>(1))
    {
      return BesselK1_Series::AtZero(x);
    }
    else
    {
      // Do the recursions.
      e_float Kn_m1 = BesselK0_Series::AtZero(x);
      e_float Kn    = BesselK1_Series::AtZero(x);

      e_float Kn_p1;
      
      for(INT32 nk = static_cast<INT32>(1); nk < nn; nk++)
      {
        Kn_p1 = (((static_cast<INT32>(2) * nk) / x) * Kn) + Kn_m1;

        Kn_m1 = Kn;
        Kn    = Kn_p1;
      }

      return Kn_p1;
    }
  }
  else
  {
    // Check for the possible convergence of the large-argument hypergeomentric expansion.

    const double xd = ef::to_double(x);
    const double vd = static_cast<double>(nn);

    const bool bo_converge_large = (vd < 30000.0) && xd > BesselConvergenceLimits::LowerLimitLargeX_Kv(vd);

    return bo_converge_large ? BesselKv_Series::AtInfinity(e_float(nn), x)
                             : BesselKv_Series::AtTransition(e_float(nn), x);
  }
}

