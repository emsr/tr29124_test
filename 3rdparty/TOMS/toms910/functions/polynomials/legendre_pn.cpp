
#include <numeric>

#include <e_float/e_float.h>
#include <functions/constants/constants.h>
#include <functions/elementary/elementary.h>
#include <functions/gamma/gamma.h>
#include <functions/hypergeometric/hypergeometric.h>
#include <functions/hypergeometric/hypergeometric_util.h>
#include <functions/polynomials/polynomials.h>
#include <functions/tables/tables.h>
#include <utility/util_coefficient_expansion.h>
#include <utility/util_interpolate.h>
#include <utility/util_trapezoid.h>

namespace LegendreP_Series
{
  static bool    AsymptoticConvergence     (const INT32 n, const e_float& x);
  static bool    NearOneConvergence        (const INT32 n, const e_float& x);
  static bool    OlverNonUniformConvergence(const INT32 n, const e_float& x);
  static e_float AtIdenticallyZero         (const INT32 n);
  static e_float AtOne                     (const INT32 n, const e_float& x);
  static e_float AtInfinity                (const INT32 n, const e_float& x);
  static e_float AtOlverNonUniform         (const INT32 n, const e_float& x);
}

static bool LegendreP_Series::AsymptoticConvergence(const INT32 n, const e_float& x)
{
  const double minus_n_over_two = static_cast<double>(static_cast<double>(-n) / static_cast<double>(2.0));

  const std::tr1::array<double, static_cast<std::size_t>(2u)> n_data =
  {{
    static_cast<double>(0.5) + minus_n_over_two,
    minus_n_over_two
  }};

  const double xd           = ef::to_double(x);
  const double half_minus_n = static_cast<double>(static_cast<double>(0.5) - static_cast<double>(n));

  return HypergeometricUtil::AsympConverge(std::deque<double>(n_data.begin(), n_data.end()),
                                           std::deque<double>(static_cast<std::size_t>(1u), half_minus_n),
                                           static_cast<double>(-1.0) / (xd * xd));
}

static bool LegendreP_Series::NearOneConvergence(const INT32 n, const e_float& x)
{
  static const double d1 = static_cast<double>(1.0);

  const std::tr1::array<double, static_cast<std::size_t>(2u)> n_data =
  {{
    static_cast<double>(-n),
    static_cast<double>(n) + d1
  }};

  const double xd = ef::to_double(x);

  return HypergeometricUtil::AsympConverge(std::deque<double>(n_data.begin(), n_data.end()),
                                           std::deque<double>(static_cast<std::size_t>(1u), d1),
                                           (static_cast<double>(1.0) - xd) / static_cast<double>(2.0));
}

static bool LegendreP_Series::OlverNonUniformConvergence(const INT32 n, const e_float& x)
{
  if(x > ef::one())
  {
    return false;
  }

  const e_float alpha           = ef::acos(x);
  const e_float two_sin_alpha   = ef::sin(alpha) * static_cast<INT32>(2);
  const INT32   two_sin_alpha_n = ef::to_int32(two_sin_alpha * n);

  return two_sin_alpha_n > static_cast<INT32>(19999);
}

static e_float LegendreP_Series::AtIdenticallyZero(const INT32 n)
{
  const bool n_is_even = static_cast<INT32>(n % static_cast<INT32>(2)) == static_cast<INT32>(0);

  if(n_is_even)
  {
    const INT32 n_over_two = static_cast<INT32>(n / static_cast<INT32>(2));
    const bool  b_negate   = static_cast<INT32>(n_over_two % static_cast<INT32>(2)) != static_cast<INT32>(0);

    const UINT32 n_minus_one = static_cast<UINT32>(n - static_cast<UINT32>(1u));

    const e_float Pn = ef::factorial2(n_minus_one) / ef::factorial2(static_cast<UINT32>(n));

    return !b_negate ? Pn : -Pn; 
  }
  else
  {
    return ef::zero();
  }
}

static e_float LegendreP_Series::AtOne(const INT32 n, const e_float& x)
{
  const e_float nf = e_float(n);

  return ef::hyperg_2f1(-nf, nf + ef::one(), ef::one(), (ef::one() - x) / static_cast<INT32>(2));
}

static e_float LegendreP_Series::AtInfinity(const INT32 n, const e_float& x)
{
  const e_float x_minus_one          = x - ef::one();
  const e_float two_over_one_minus_x = ef::two() / -x_minus_one;

  e_float two_over_one_minus_x_powk_over_kfact = two_over_one_minus_x;

  e_float pochham_minus_n   = e_float(-n);
  e_float minus_n_p         = pochham_minus_n;
  e_float pochham_minus_2n  = e_float(-(n * static_cast<INT32>(2)));
  e_float minus_2n_p        = pochham_minus_2n;

  e_float sum = ef::one() + (((pochham_minus_n * pochham_minus_n) / pochham_minus_2n) * two_over_one_minus_x_powk_over_kfact);

  // Use the series expansion of LegendreP(n, x) for integer n and large x.
  // http://functions.wolfram.com/HypergeometricFunctions/LegendrePGeneral/06/01/04/

  for(INT32 k = static_cast<INT32>(2); k <= n; k++)
  {
    two_over_one_minus_x_powk_over_kfact *= two_over_one_minus_x;
    two_over_one_minus_x_powk_over_kfact /= k;

    pochham_minus_n  *= ++minus_n_p;
    pochham_minus_2n *= ++minus_2n_p;

    const e_float term = (((pochham_minus_n * pochham_minus_n) / pochham_minus_2n) * two_over_one_minus_x_powk_over_kfact);

    if(term.order() < -ef::tol())
    {
      break;
    }

    sum += term;
  }
  
  const e_float x_minus_one_powk_over_nfact = ef::pown(x_minus_one, static_cast<INT32>(n)) / ef::factorial(n);

  const e_float factor = (ef::pow2(static_cast<INT64>(n)) * ef::gamma(e_float(n) + ef::half())) * x_minus_one_powk_over_nfact;

  return (factor * sum) / ef::sqrt_pi();
}

static e_float LegendreP_Series::AtOlverNonUniform(const INT32 n, const e_float& x)
{
  const e_float alpha              = ef::acos(x);
  const e_float N                  = e_float(n);
  const e_float sin_alpha          = ef::sin(alpha);
  const e_float two_sin_alpha      = sin_alpha * static_cast<INT32>(2);
  const e_float alpha_plus_pi_half = alpha + ef::pi_half();

  e_float sum      = ef::zero();
  e_float alpha_ns = ((N + ef::half()) * alpha) + ((N - ef::quarter()) * ef::pi());

  e_float binomial_minus_half_s   = ef::binomial(-ef::half(), ef::zero());
  e_float binomial_s_minus_half_n = ef::binomial(-ef::half(), N);
  e_float two_sin_alpha_pow_s     = ef::one();

  for(INT32 s = static_cast<INT32>(0); s < ef::max_iteration(); s++)
  {
    const e_float term = ((binomial_minus_half_s * binomial_s_minus_half_n) * ef::cos(alpha_ns)) / two_sin_alpha_pow_s;

    const INT64 order_check = static_cast<INT64>(term.order() - sum.order());

    if((s > static_cast<INT32>(20)) && (order_check < -ef::tol()))
    {
      break;
    }

    two_sin_alpha_pow_s *= two_sin_alpha;

    sum      += term;
    alpha_ns -= alpha_plus_pi_half;

    // Use the properties of gamma function recursion to increment the binomials.
    const e_float s_plus_half = e_float(s) + ef::half();

    binomial_minus_half_s   *= -s_plus_half / static_cast<INT32>(s + static_cast<INT32>(1));
    binomial_s_minus_half_n *= (s_plus_half / (s_plus_half - N));
  }

  return sum * ef::sqrt(ef::two() / sin_alpha);
}

e_float ef::legendre_p(const INT32 n, const e_float& x)
{
  if(n < static_cast<INT32>(0))
  {
    // Handle negative degree with v < -1 using
    // Abramowitz & Stegun 8.2.1, P_[-v-1](z) = P_[v](z).
    return legendre_p(static_cast<INT32>(static_cast<INT32>(-n) - static_cast<INT32>(1)), x);
  }

  if(ef::iszero(x))
  {
    return LegendreP_Series::AtIdenticallyZero(n);
  }

  if(ef::isneg(x))
  {
    const bool is_odd_degree = static_cast<INT32>(n % static_cast<INT32>(2)) != static_cast<INT32>(0);

    const e_float Pn = legendre_p(n, -x);

    return !is_odd_degree ? Pn : -Pn;
  }

  if(ef::isone(x)) { return ef::one(); }

  if(n == static_cast<INT32>(0)) { return ef::one(); }
  if(n == static_cast<INT32>(1)) { return x; }

  // Check for convergence near one.
  if(LegendreP_Series::NearOneConvergence(n, x))
  {
    return LegendreP_Series::AtOne(n, x);
  }

  // Check for non-uniform asymptotic convergence.
  if((x < ef::one()) && LegendreP_Series::OlverNonUniformConvergence(n, x))
  {
    return LegendreP_Series::AtOlverNonUniform(n, x);
  }

  // Check for asymptotic convergence.
  if(LegendreP_Series::AsymptoticConvergence(n, x))
  {
    return LegendreP_Series::AtInfinity(n, x);
  }

  // Use forward recursion to calculate Legendre Polynomials of degree 2 or higher.
  e_float Pn_m2 = ef::legendre_p(static_cast<INT32>(0), x);
  e_float Pn_m1 = ef::legendre_p(static_cast<INT32>(1), x);
  e_float Pn;

  for(INT32 j = static_cast<INT32>(2); j <= n; j++)
  {
    const INT32 n_minus_one     = static_cast<INT32>(j - static_cast<INT32>(1));
    const INT32 two_n_minus_one = static_cast<INT32>(n_minus_one + j);

    Pn = (((two_n_minus_one * x) * Pn_m1) - (n_minus_one * Pn_m2)) / j;

    Pn_m2 = Pn_m1;
    Pn_m1 = Pn;
  }

  return Pn;
}
