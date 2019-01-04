
#include <functions/complex/e_float_complex.h>
#include <functions/constants/constants.h>
#include <functions/elementary/elementary.h>
#include <functions/gamma/gamma.h>
#include <functions/gamma/polygamma.h>
#include <functions/integer/integer.h>
#include <functions/tables/tables.h>
#include <utility/util_digit_scale.h>
#include <functions/zeta/zeta.h>

namespace PolyGamma_Series
{
  static e_float Reflection      (const INT32 n, const e_float& x);
  static e_float AtZero          (const INT32 n, const e_float& x);
  static e_float AtTransitionPlus(const INT32 n, const e_float& x);
  static e_float AtInfinityPlus  (const INT32 n, const e_float& x);
}

static e_float PolyGamma_Series::Reflection(const INT32 n, const e_float& x)
{
  // Use reflection or recursion as shown at Wolfram's Function Site.
  // http://functions.wolfram.com/GammaBetaErf/PolyGamma2/16/01/01/

  const INT64 n_plus_one = static_cast<INT64>(n + static_cast<INT32>(1));
  const bool  n_is_odd   = static_cast<INT32>(n % static_cast<INT32>(2)) != static_cast<INT32>(0);

  if(n > static_cast<INT32>(20))
  {
    // Use downward recursion.
    const e_float xx = ef::one() + ef::decimal_part(x);
    const INT32   m  = static_cast<INT32>(static_cast<INT64>(1) - ef::to_int64(x));

    e_float sum = ef::zero();

    for(INT32 k = static_cast<INT32>(1); k <= m; k++)
    {
      sum += ef::one() / ef::pown(xx - k, n_plus_one);
    }

    const e_float sum_term = ef::factorial(static_cast<UINT32>(n)) * sum;
    const e_float PG       = ef::poly_gamma(n, xx); 

    return !n_is_odd ? PG - sum_term : PG + sum_term;
  }
  else
  {
    // Use reflection which involves Stirling's numbers of the second kind and the
    // use of complex numbers.
    const e_float xx = -x;

    bool b_neg_term = true;

    ef_complex sum = ef::zero();
    
    e_float k_fact    = ef::one();
    e_float two_pow_k = ef::one();

    const ef_complex i_cot_pi_x                 = ef_complex(ef::zero(), ef::cot(ef::pi() * xx));
    const ef_complex i_cot_pi_x_plpus_one       = i_cot_pi_x + ef::one();
          ef_complex i_cot_pi_x_plpus_one_pow_k = ef::one();

    for(INT32 k = static_cast<INT32>(1); k <= n; k++)
    {
      k_fact                     *= k;
      i_cot_pi_x_plpus_one_pow_k *= i_cot_pi_x_plpus_one;
      two_pow_k                  *= static_cast<INT32>(2);

      const ef_complex term = ((k_fact * ef::stirling2(static_cast<UINT32>(n), static_cast<UINT32>(k))) * i_cot_pi_x_plpus_one_pow_k) / two_pow_k;

      !b_neg_term ? sum += term : sum -= term;

      b_neg_term = !b_neg_term;
    }

    const e_float    psi_n_term  = !n_is_odd ? ef::poly_gamma(n, xx) : -ef::poly_gamma(n, xx);
    const e_float    n_z_term    = ef::factorial(static_cast<UINT32>(n)) / ef::pown(xx, n_plus_one);
    const ef_complex two_pi_i    = ef_complex(ef::zero(), ef::two_pi());
    const ef_complex sum_term    = -(((efz::pown(two_pi_i, n_plus_one) / static_cast<INT32>(2)) * (i_cot_pi_x - ef::one())) * sum);

    return ((psi_n_term + n_z_term) + sum_term).real();
  }
}

static e_float PolyGamma_Series::AtZero(const INT32 n, const e_float& x)
{
  // Use a series expansion for x near zero which uses poly_gamma(m, 1) which,
  // in turn, uses the Riemann zeta function for integer arguments.
  // http://functions.wolfram.com/GammaBetaErf/PolyGamma2/06/01/03/01/02/

  const bool b_negate = (static_cast<INT32>(n % static_cast<INT32>(2)) == static_cast<INT32>(0u));

  const e_float n_fact               = ef::factorial(static_cast<UINT32>(n));
  const e_float z_pow_n_plus_one     = ef::pown(x, static_cast<INT64>(n + static_cast<INT32>(1)));
  const e_float n_fact_over_pow_term = n_fact / z_pow_n_plus_one;
  const e_float term0                = !b_negate ? n_fact_over_pow_term : -n_fact_over_pow_term;

        e_float one_over_k_fact   = ef::one();
        e_float z_pow_k           = ef::one();
        e_float k_plus_n_fact     = ef::factorial(static_cast<UINT32>(n));
        INT32   k_plus_n_plus_one = static_cast<INT32>(n + static_cast<INT32>(1));
  const e_float pg_kn             = k_plus_n_fact * ef::riemann_zeta(k_plus_n_plus_one);
        bool    b_neg_term        = (static_cast<INT32>(n % static_cast<INT32>(2)) == static_cast<INT32>(0u));
        e_float sum               = !b_neg_term ? pg_kn : -pg_kn;

  for(INT32 k = 1; k < ef::max_iteration(); k++)
  {
    k_plus_n_fact   *= k_plus_n_plus_one++;
    one_over_k_fact /= k;
    z_pow_k         *= x;

    const e_float pg = k_plus_n_fact * ef::riemann_zeta(k_plus_n_plus_one);

    const e_float term = (pg * z_pow_k) * one_over_k_fact;

    const INT64 order_check = static_cast<INT64>(term.order() - sum.order());

    if((k > static_cast<INT32>(20)) && (order_check < -ef::tol()))
    {
      break;
    }

    b_neg_term = !b_neg_term;

    !b_neg_term ? sum += term : sum -= term;
  }

  return term0 + sum;
}

static e_float PolyGamma_Series::AtTransitionPlus(const INT32 n, const e_float& x)
{
  // Use Euler-Maclaurin summation.

  // Use N = (0.4 * digits) + (4 * n)
  static const INT64 d4d  = static_cast<INT64>(static_cast<double>(0.4) * static_cast<double>(static_cast<INT32>(ef::tol())));
         const INT64 N4dn = static_cast<INT64>(d4d + static_cast<INT64>(static_cast<INT64>(4) * static_cast<INT64>(n)));
  static const INT64 n32m = static_cast<INT64>(std::numeric_limits<INT32>::max());
         const INT32 N    = static_cast<INT32>(std::min(N4dn, n32m));
         const INT32 m    = n;

  const INT64 minus_m_minus_one = static_cast<INT64>(static_cast<INT32>(-m) - static_cast<INT32>(1));

  e_float z                              = x;
  e_float sum0                           = ef::zero();
  e_float z_plus_k_pow_minus_m_minus_one = ef::zero();

  for(INT32 k = static_cast<INT32>(1); k <= N; k++)
  {
    z_plus_k_pow_minus_m_minus_one = ef::pown(z, minus_m_minus_one);

    sum0 += z_plus_k_pow_minus_m_minus_one;

    ++z;
  }

  const e_float one_over_z_plus_N_pow_minus_m           = ef::pown(z, static_cast<INT64>(-m));
  const e_float one_over_z_plus_N_pow_minus_m_minus_one = one_over_z_plus_N_pow_minus_m / z;

  const e_float term0 = one_over_z_plus_N_pow_minus_m_minus_one / static_cast<INT32>(2);
  const e_float term1 = one_over_z_plus_N_pow_minus_m           / m;

        e_float sum1                                      = ef::zero();
        e_float one_over_two_k_fact                       = ef::half();
        INT32   mk                                        = static_cast<INT32>(m + static_cast<INT32>(1));
        e_float am                                        = e_float(mk);
  const e_float one_over_z_plus_N_squared                 = ef::one() / (z * z);
        e_float one_over_z_plus_N_pow_minus_m_minus_two_k = one_over_z_plus_N_pow_minus_m * one_over_z_plus_N_squared;

  for(INT32 k = static_cast<INT32>(1); k < ef::max_iteration(); k++)
  {
    const INT32 two_k = static_cast<INT32>(static_cast<INT32>(2) * k);

    const e_float term = ((ef::bernoulli(static_cast<UINT32>(two_k)) * am) * one_over_two_k_fact) * one_over_z_plus_N_pow_minus_m_minus_two_k;

    const INT64 order_check = term.order() - sum1.order();

    if((k > static_cast<INT32>(20)) && (order_check < -ef::tol()))
    {
      break;
    }

    sum1 += term;

    one_over_two_k_fact /= static_cast<INT32>(two_k + static_cast<INT32>(1));
    one_over_two_k_fact /= static_cast<INT32>(two_k + static_cast<INT32>(2));

    am *= static_cast<INT32>(++mk);
    am *= static_cast<INT32>(++mk);

    one_over_z_plus_N_pow_minus_m_minus_two_k *= one_over_z_plus_N_squared;
  }

  const e_float pg = (((sum0 + term0) + term1) + sum1) * ef::factorial(static_cast<UINT32>(m));

  const bool b_negate = (static_cast<INT32>(m % static_cast<INT32>(2)) == static_cast<INT32>(0u));

  return (!b_negate ? pg : -pg);
}

static e_float PolyGamma_Series::AtInfinityPlus(const INT32 n, const e_float& x)
{
  // Use an asymptotic series expansion.
  // http://functions.wolfram.com/GammaBetaErf/PolyGamma2/06/02/

  const bool b_negate = (static_cast<INT32>(n % static_cast<INT32>(2)) == static_cast<INT32>(0u));

  const e_float n_minus_one_fact            = ef::factorial(static_cast<UINT32>(n - static_cast<INT32>(1)));
  const e_float nn                          = e_float(n);
  const e_float n_fact                      = n_minus_one_fact * nn;
  const e_float one_over_z                  = ef::one() / x;
  const e_float one_over_z2                 = one_over_z * one_over_z;
  const e_float one_over_z_pow_n            = ef::one() / ef::pown(x, static_cast<INT64>(n));
        e_float one_over_x_pow_two_k_plus_n = one_over_z_pow_n * one_over_z2;
        e_float two_k_plus_n_minus_one      = nn + ef::one();
        e_float two_k_plus_n_minus_one_fact = n_fact * e_float(static_cast<UINT32>(n + static_cast<INT32>(1)));
        e_float one_over_two_k_fact         = ef::half();
        e_float sum                         = (  (ef::bernoulli(static_cast<UINT32>(2u)) * two_k_plus_n_minus_one_fact)
                                               * (one_over_two_k_fact * one_over_x_pow_two_k_plus_n));

  // Perform the Bernoulli series expansion.
  for(INT32 two_k = static_cast<INT32>(4); two_k < ef::max_iteration(); two_k += static_cast<INT32>(2))
  {
    one_over_x_pow_two_k_plus_n *= one_over_z2;
    two_k_plus_n_minus_one_fact *= ++two_k_plus_n_minus_one;
    two_k_plus_n_minus_one_fact *= ++two_k_plus_n_minus_one;
    one_over_two_k_fact         /= static_cast<INT32>(two_k * static_cast<INT32>(two_k - static_cast<INT32>(1)));

    const e_float term = (  (ef::bernoulli(two_k) * two_k_plus_n_minus_one_fact)
                          * (one_over_two_k_fact * one_over_x_pow_two_k_plus_n));

    const INT64 order_check = term.order() - sum.order();

    if((two_k > static_cast<INT32>(20)) && (order_check < -ef::tol()))
    {
      break;
    }

    sum += term;
  }

  sum += ((((n_minus_one_fact * (nn + (x * static_cast<INT32>(2)))) * one_over_z_pow_n) * one_over_z) / static_cast<INT32>(2));

  return !b_negate ? sum : -sum;
}

e_float ef::poly_gamma(const e_float& x)
{
  if(!ef::isfinite(x))
  {
    return std::numeric_limits<e_float>::quiet_NaN();
  }

  const bool b_neg = ef::isneg(x);

  // Make a local, unsigned copy of the input argument.
  e_float xx(!b_neg ? x : -x);

  // Make a local copy of the unscaled, unsigned argument.
  const e_float xx_unscaled = xx;

  e_float psi;

  if(ef::isint(xx) && (xx < ef::thousand()))
  {
    // Special handling for small, pure-integer arguments (0, 1, 2, 3...).
    const INT32 nx = ef::to_int32(xx);

    if(nx == static_cast<UINT32>(0u))
    {
      return std::numeric_limits<e_float>::quiet_NaN();
    }

    psi = -ef::euler_gamma();

    for(INT32 k = static_cast<INT32>(1); k < nx; k++)
    {
      psi += ef::one() / k;
    }
  }
  else if(isint(xx * static_cast<INT32>(2)) && (xx < ef::thousand()))
  {
    // Special handling for small half-integer arguments (1/2, 1 + 1/2, 2 + 1/2, 3 + 1/2...).
    const INT32 nx = ef::to_int32(xx);

    psi = -ef::euler_gamma() - (ef::ln2() * static_cast<INT32>(2));
    
    for(INT32 k = static_cast<INT32>(1); k <= nx; k++)
    {
      psi += ef::two() / static_cast<INT32>((k * static_cast<INT32>(2)) - static_cast<INT32>(1));
    }
  }
  else
  {
    // Check if the argument should be scaled up for the Bernoulli series expansion.
    static const INT32   min_arg_n = static_cast<INT32>(static_cast<double>(240.0) * Util::DigitScale());
    static const e_float min_arg_x = e_float(min_arg_n);

    const INT32 n_recur = ((xx < min_arg_x) ? static_cast<INT32>((min_arg_n - ef::to_int32(xx)) + 1u)
                                            : static_cast<INT32>(0));

    // Scale the argument up and use downward recursion later for the final result.
    if(n_recur != static_cast<INT32>(0))
    {
      xx += n_recur;
    }

    const e_float one_over_x        = ef::one() / xx;
    const e_float one_over_x2       = one_over_x * one_over_x;
          e_float one_over_x_pow_2k = one_over_x2;
          e_float sum               = one_over_x2 / static_cast<INT32>(12);

    // Perform the Bernoulli series expansion.
    for(INT32 k2 = static_cast<INT32>(4); k2 < ef::max_iteration(); k2 += static_cast<INT32>(2))
    {
      one_over_x_pow_2k *= one_over_x2;

      const e_float term = (ef::bernoulli(k2) * one_over_x_pow_2k) / k2;

      if(term.order() < -ef::tol())
      {
        break;
      }

      sum += term;
    }

    psi = (ef::log(xx) - (one_over_x / static_cast<INT32>(2))) - sum;

    // Rescale the result using downward recursion if necessary.
    if(n_recur != static_cast<UINT32>(0u))
    {
      e_float rescale_sum = ef::one() / xx_unscaled;
    
      for(INT32 k = static_cast<INT32>(1); k < n_recur; k++)
      {
        rescale_sum += ef::one() / (xx_unscaled + k);
      }

      psi -= rescale_sum;
    }
  }

  // Return the result, accounting for possible negative arguments.
  return !b_neg ?  psi
                : (psi + (ef::one() / xx_unscaled)) + (ef::pi() * ef::cot(ef::pi() * xx_unscaled));
}

e_float ef::poly_gamma(const INT32 n, const e_float& x)
{
  if(!ef::isfinite(x))
  {
    return std::numeric_limits<e_float>::quiet_NaN();
  }

  if(n == static_cast<INT32>(0u))
  {
    return ef::poly_gamma(x);
  }

  const bool b_neg = isneg(x);

  if(n < static_cast<INT32>(0u))
  {
    return std::numeric_limits<e_float>::quiet_NaN();
  }

  if(b_neg)
  {
    // Use reflection for poly_gamma for negative argument.
    return PolyGamma_Series::Reflection(n, x);
  }

  if(isint(x))
  {
    if(b_neg || iszero(x))
    {
      return std::numeric_limits<e_float>::quiet_NaN();
    }
    else
    {
      const INT64 nx = ef::to_int64(x);

      // Possible special handling for the pure integer one.
      if(nx == static_cast<INT64>(1))
      {
        const bool    b_negate = (static_cast<INT32>(n % static_cast<INT32>(2)) == static_cast<INT32>(0u));
        const e_float pg       = ef::factorial(static_cast<UINT32>(n)) * ef::riemann_zeta(static_cast<INT32>(n + static_cast<INT32>(1)));

        return (!b_negate ? pg : -pg);
      }
    }
  }

  // Make a local, unsigned copy of the input argument.
  e_float xx(!b_neg ? x : -x);

  // Determine which series or algorithm should be used based on the size of the argument.

  static const INT32 scaled_limit_n = static_cast<INT32>(static_cast<double>(1000.0) * Util::DigitScale());
  static const INT32 scaled_limit_x = static_cast<INT32>(static_cast<double>( 300.0) * Util::DigitScale());

  const INT32   min_arg_n = (n < scaled_limit_n ? scaled_limit_x : static_cast<INT32>(scaled_limit_x + (n - scaled_limit_n)));
  const e_float min_arg_x = e_float(min_arg_n);

  if(xx < ef::one())
  {
    return PolyGamma_Series::AtZero(n, xx);
  }
  else
  {
    return ((xx > min_arg_x) ? PolyGamma_Series::AtInfinityPlus(n, xx)
                             : PolyGamma_Series::AtTransitionPlus(n, xx));
  }
}
