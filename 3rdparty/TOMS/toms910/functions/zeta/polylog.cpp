
#include <functions/complex/e_float_complex.h>
#include <functions/constants/constants.h>
#include <functions/elementary/elementary.h>
#include <functions/gamma/gamma.h>
#include <functions/gamma/polygamma.h>
#include <functions/hypergeometric/hypergeometric.h>
#include <functions/hypergeometric/hypergeometric_util.h>
#include <functions/integer/integer.h>
#include <functions/tables/tables.h>
#include <functions/zeta/polylog.h>
#include <functions/zeta/zeta.h>
#include <utility/util_power_j_pow_x.h>

namespace PolyLog_Series
{
  static e_float AtZeroForPositiveN       (const INT32 n, const e_float& x);
  static e_float AtZeroForNegativeN       (const INT32 N, const e_float& x);
  static e_float AtOneForNegativeN        (const INT32 N, const e_float& x);
  static e_float AtInfinityMinus          (const INT32 n, const e_float& x);
  static e_float AtInfinityForNegativeN   (const INT32 N, const e_float& x);
  static e_float AtStirlingS2ForNegativeN (const INT32 N, const e_float& x);
  static e_float AtReflectNegativeArgument(const INT32 n, const e_float& x);


  template<typename T> inline T AtOneForPositiveN(const INT32 n, const e_float& x)
  {
    // Use a series expansion for x near one which uses the Riemann zeta function for
    // integer arguments.
    // http://functions.wolfram.com/ZetaFunctionsandPolylogarithms/PolyLog/06/01/02/01/02/
    // The series calculation involves multiple evaluations of the Riemann zeta function
    // for integer values. The implementation of the Riemann zeta function for integer
    // arguments can be found in zeta/zeta.cpp and it is optimized for speed.

    if(n < static_cast<INT32>(2)) { return ef::zero(); }

    const e_float logz        = ef::log(x);
    const INT32   n_minus_one = static_cast<INT32>(n - static_cast<INT32>(1));

    e_float sum1 = ef::zero();

    e_float one_over_j_fact = ef::one();
    e_float logz_pow_j      = ef::one();

    for(INT32 j = static_cast<INT32>(1); j < n_minus_one; j++)
    {
      one_over_j_fact /= j;
      logz_pow_j      *= logz;

      sum1 += (ef::riemann_zeta(static_cast<INT32>(n - j)) * logz_pow_j) * one_over_j_fact;
    }

    e_float sum2 = ef::zero();

    one_over_j_fact /= n_minus_one;

    logz_pow_j      *= logz;

    for(INT32 j = n; j < ef::max_iteration(); j++)
    {
      one_over_j_fact /= j;
      logz_pow_j      *= logz;

      const e_float term = (ef::riemann_zeta(static_cast<INT32>(n - j)) * logz_pow_j) * one_over_j_fact;

      if(term.order() < -ef::tol())
      {
        break;
      }

      sum2 += term;
    }

    // The value of term0 defined below can be either real or complex.
    {
      using ef::log;
      using efz::log;

      const e_float psi_n            = ef::poly_gamma(e_float(n));
      const e_float n_minus_one_fact = ef::factorial(static_cast<UINT32>(n_minus_one));
      const T       log_log_term     = log(T(-logz));

      const T term0 = (((-log_log_term + psi_n) + ef::euler_gamma()) / n_minus_one_fact) * ef::pown(logz, static_cast<INT64>(n_minus_one));

      return (term0 + (sum1 + sum2)) + ef::riemann_zeta(n);
    }
  }
}

static e_float PolyLog_Series::AtZeroForPositiveN(const INT32 n, const e_float& x)
{
  // Use the Taylor series for |x| near zero.
  // http://functions.wolfram.com/ZetaFunctionsandPolylogarithms/PolyLog/06/01/01/01/02/

  return x * ef::hyperg_pfq(std::deque<e_float>(static_cast<std::size_t>(n + 1), ef::one()),
                            std::deque<e_float>(static_cast<std::size_t>(n),     ef::two()),
                            x);
}

static e_float PolyLog_Series::AtZeroForNegativeN(const INT32 N, const e_float& x)
{
  // Use a Taylor series from Wolfram.
  // http://functions.wolfram.com/ZetaFunctionsandPolylogarithms/PolyLog/06/01/01/01/01/

  std::map<UINT32, e_float> k_pow_n_prime_factor_map;

  const e_float nf = e_float(N);

  e_float sum     = x;
  e_float x_pow_k = x;

  for(INT32 k = static_cast<INT32>(2); k < ef::max_iteration(); k++)
  {
    x_pow_k *= x;

    const e_float term = Util::j_pow_x(static_cast<UINT32>(k), nf, k_pow_n_prime_factor_map) * x_pow_k;

    const INT64 order_check = term.order() - sum.order();

    if((k > static_cast<INT32>(20)) && (order_check < -ef::tol()))
    {
      break;
    }

    sum += term;
  }

  return sum;
}

static e_float PolyLog_Series::AtOneForNegativeN(const INT32 N, const e_float& x)
{
  // Use a series called "other series" from Wolfram's Function Site.
  // http://functions.wolfram.com/ZetaFunctionsandPolylogarithms/PolyLog/06/03/
  const INT32 n_plus_one = static_cast<INT32>(N + static_cast<INT32>(1));

        e_float one_over_k_fact = ef::one();
  const e_float log_x           = ef::log(x);
        e_float log_x_pow_k     = ef::one();

  e_float sum = ef::bernoulli(static_cast<UINT32>(n_plus_one)) / n_plus_one;
  
  for(INT32 k = static_cast<INT32>(1); k < ef::max_iteration(); k++)
  {
    one_over_k_fact /= k;
    log_x_pow_k     *= log_x;

    const e_float zeta_term = ef::riemann_zeta(-N - k);

    if(!ef::iszero(zeta_term))
    {
      const e_float term = (zeta_term * log_x_pow_k) * one_over_k_fact;

      const INT64 order_check = term.order() - sum.order();

      if((k > static_cast<INT32>(20)) && (order_check < -ef::tol()))
      {
        break;
      }

      sum += term;
    }
  }

  const e_float term0 = ef::factorial(static_cast<UINT32>(N)) * ef::pown(-log_x, static_cast<INT64>(-n_plus_one));

  return term0 - sum;
}

static e_float PolyLog_Series::AtInfinityMinus(const INT32 n, const e_float& x)
{
  // Use an exponential series from Wolfram's Function Site.
  // http://functions.wolfram.com/ZetaFunctionsandPolylogarithms/PolyLog/06/01/03/01/02/

  const e_float z   = ef::log(-x);
  const e_float zsq = z * z;

  e_float two_pow_one_minus_2k = ef::two();
  e_float z_pow_n_minus_2k     = ef::pown(z, static_cast<INT64>(n));
  e_float n_minus_two_k_fact   = ef::factorial(static_cast<UINT32>(n));
  INT32   n_minus_two_k        = n;

  e_float sum2 = ef::zero();

  const INT32 n_over_two = static_cast<INT32>(n / static_cast<INT32>(2));

  for(INT32 k = static_cast<INT32>(0); k <= n_over_two; k++)
  {
    const INT32 two_k = static_cast<INT32>(static_cast<INT32>(2) * k);

    sum2 += (((ef::one() - two_pow_one_minus_2k) * z_pow_n_minus_2k) * ef::riemann_zeta(two_k)) / n_minus_two_k_fact;

    if(k < n_over_two)
    {
      two_pow_one_minus_2k /= static_cast<INT32>(4);

      z_pow_n_minus_2k /= zsq;
    
      n_minus_two_k_fact /= static_cast<INT32>(n_minus_two_k--);
      n_minus_two_k_fact /= static_cast<INT32>(n_minus_two_k--);
    }
  }

  const e_float exp_minus_z = ef::exp(-z);

  e_float exp_minus_kz = exp_minus_z;

  e_float sum1 = -exp_minus_kz;

  bool is_neg_term = true;

  for(INT32 k = static_cast<INT32>(2); k < ef::max_iteration(); k++)
  {
    is_neg_term = !is_neg_term;

    exp_minus_kz *= exp_minus_z;

    const e_float term = exp_minus_kz / ef::pown(e_float(k), static_cast<INT64>(n));

    const INT64 order_check = term.order() - sum1.order();

    if((k > static_cast<INT32>(20)) && (order_check < -ef::tol()))
    {
      break;
    }

    !is_neg_term ? sum1 += term : sum1 -= term;
  }

  if(static_cast<INT32>(n % static_cast<INT32>(2)) == static_cast<INT32>(0)) { sum1 = -sum1; }

  return sum1 - (sum2 * static_cast<INT32>(2));
}

static e_float PolyLog_Series::AtInfinityForNegativeN(const INT32 N, const e_float& x)
{
  // Use a Taylor series from Wolfram's Function Site.
  // http://functions.wolfram.com/ZetaFunctionsandPolylogarithms/PolyLog/06/01/03/01/02/

  std::map<UINT32, e_float> k_pow_n_prime_factor_map;

  const e_float nf = e_float(N);

  e_float sum     = x;
  e_float x_pow_k = x;

  for(INT32 k = static_cast<INT32>(2); k < ef::max_iteration(); k++)
  {
    x_pow_k *= x;

    const e_float term = Util::j_pow_x(static_cast<UINT32>(k), nf, k_pow_n_prime_factor_map) / x_pow_k;

    const INT64 order_check = term.order() - sum.order();

    if((k > static_cast<INT32>(20)) && (order_check < -ef::tol()))
    {
      break;
    }

    sum += term;
  }

  const INT32 n_minus_one = static_cast<INT32>(N - static_cast<INT32>(1));

  const bool b_negate = static_cast<INT32>(n_minus_one % static_cast<INT32>(2)) != static_cast<INT32>(0);

  return !b_negate ? sum : -sum;
}

static e_float PolyLog_Series::AtStirlingS2ForNegativeN(const INT32 N, const e_float& x)
{
  // Use a series involving Stirling's numbers of the second kind (Ref. Wiki: Wood 1992).
  const INT32 n_plus_one = static_cast<INT32>(N + static_cast<INT32>(1));

        e_float k_minus_one_fact  = ef::one();
  const e_float one_minus_x       = ef::one() - x;
        e_float one_minus_x_pow_k = ef::one();

  e_float sum = ef::zero();

  bool b_neg_term = static_cast<INT32>(N % static_cast<INT32>(2)) != static_cast<INT32>(0);

  for(INT32 k = static_cast<INT32>(1); k <= n_plus_one; k++)
  {
    const INT32 k_minus_one = static_cast<INT32>(k - static_cast<INT32>(1));

    if(k_minus_one > static_cast<INT32>(0))
    {
      k_minus_one_fact *= k_minus_one;
    }

    one_minus_x_pow_k *= one_minus_x;

    const e_float term = (k_minus_one_fact * ef::stirling2(static_cast<UINT32>(n_plus_one), static_cast<UINT32>(k))) / one_minus_x_pow_k;

    !b_neg_term ? sum += term : sum -= term;

    b_neg_term = !b_neg_term;
  }

  return sum;
}

static e_float PolyLog_Series::AtReflectNegativeArgument(const INT32 n, const e_float& x)
{
  const e_float xx = -x;

  const e_float two_pow_one_minus_n = ef::pow2(static_cast<INT64>(static_cast<INT32>(1) - n));

  return (two_pow_one_minus_n * ef::poly_logarithm(n, xx * xx)) - ef::poly_logarithm(n, xx);
}

e_float ef::poly_logarithm(const INT32 n, const e_float& x)
{
  if(n == static_cast<INT32>(1))
  {
    return -ef::log(ef::one() - x);
  }
  else if(n == static_cast<INT32>(0))
  {
    return x / (ef::one() - x);
  }

  if(n < static_cast<INT32>(0))
  {
    if(ef::isone(x))
    {
      return std::numeric_limits<e_float>::quiet_NaN();
    }

    const INT32 N  = static_cast<INT32>(-n);

    if(N == static_cast<INT32>(1))
    {
      return x / ef::pown(ef::one() - x, static_cast<INT64>(2));
    }
    else if(N == static_cast<INT32>(2))
    {
      const e_float one_plus_x = ef::one() + x;
      return (x * one_plus_x) / ef::pown(ef::one() - x, static_cast<INT64>(3));
    }
    else if(N == static_cast<INT32>(3))
    {
      const e_float one_plus_4x_plus_xsq = ef::one() + ((x * static_cast<INT32>(4)) + (x * x));
      return (x * one_plus_4x_plus_xsq) / ef::pown(ef::one() - x, static_cast<INT64>(4));
    }
    else if(N == static_cast<INT32>(4))
    {
      const e_float one_plus_x            = ef::one() + x;
      const e_float one_plus_10x_plus_xsq = ef::one() + ((x * static_cast<INT32>(10)) + (x * x));

      return ((x * one_plus_x) * one_plus_10x_plus_xsq) / ef::pown(ef::one() - x, static_cast<INT64>(5));
    }
    else
    {
      if(N < static_cast<INT32>(65))
      {
        return PolyLog_Series::AtStirlingS2ForNegativeN(N, x);
      }

      if(ef::isneg(x))
      {
        return PolyLog_Series::AtReflectNegativeArgument(n, x);
      }

      static const e_float hundredth = ef::one() / static_cast<INT32>(100);

      const e_float delta_one = N * ef::fabs(ef::one() - x);

      if(delta_one < hundredth)
      {
        return PolyLog_Series::AtOneForNegativeN(N, x);
      }

      if(x < ef::half())
      {
        return PolyLog_Series::AtZeroForNegativeN(N, x);
      }

      if((N > static_cast<INT32>(500)) && (x > e_float(N)))
      {
        return PolyLog_Series::AtInfinityForNegativeN(N, x);
      }

      // Use reflection for positive values of x with negative values of n in the
      // transition region. Note that |x| > 1 here.
      if(x > ef::one())
      {
        const bool b_negate = static_cast<INT32>(N % static_cast<INT32>(2)) != static_cast<INT32>(1);

        const e_float PolyLogOfOneOverX = ef::poly_logarithm(n, ef::one() / x);

        return (!b_negate ? PolyLogOfOneOverX : -PolyLogOfOneOverX);
      }

      // TBD: The following expansion loses precision for large negative n.
      // TBD: Need to find expansions for large negative n in the transition region.
      return PolyLog_Series::AtStirlingS2ForNegativeN(N, x);
    }
  }

  if(ef::isneg(x))
  {
    const e_float two_pow_one_minus_n = ef::pow2(static_cast<INT64>(static_cast<INT32>(1) - n));

    if(ef::isone(-x))
    {
      // The argument is equal to minus one.
      // Use Li_n(-1) = -eta(n) with eta(n) = zeta(n) * [1 - 2^(1 - n)].
      return -ef::riemann_zeta(n) * (ef::one() - two_pow_one_minus_n);
    }
    else if(x < -ef::exp1())
    {
      // The argument is in the range -inf < x < -(e^1).
      return PolyLog_Series::AtInfinityMinus(n, x);
    }
    else if(x < -ef::one())
    {
      // The argument is in the range is -(e^1) <= x < -1.
      // The square relationship is used on the negated argument. This negated
      // argument is a positive argument 1 < |x| < e, which is larger than one
      // but less than e. Thus the values of the polylogarithms of the negated
      // argument as well as the square of the negated argument are complex numbers.
      // The imaginary parts cancel each other and the result is purely real.
      return efz::real( (PolyLog_Series::AtOneForPositiveN<ef_complex>(n,  x * x) * two_pow_one_minus_n)
                       - PolyLog_Series::AtOneForPositiveN<ef_complex>(n, -x));
    }
    else
    {
      static const double d1 = static_cast<double>(1.0);
      static const double d2 = static_cast<double>(2.0);

      // The argument is in the range -1 < x < 0.
      // Use the Taylor series or use the square relationship depending on the convergence.
      const double xd = ef::to_double(x);

      const bool bo_converge_small = HypergeometricUtil::AsympConverge(std::deque<double>(static_cast<std::size_t>(n + 1), d1),
                                                                       std::deque<double>(static_cast<std::size_t>(n),     d2),
                                                                       xd);

      return bo_converge_small ? PolyLog_Series::AtZeroForPositiveN(n, x)
                               : PolyLog_Series::AtReflectNegativeArgument(n, x);
    }
  }
  else
  {
    if(x > ef::one())
    {
      return std::numeric_limits<e_float>::quiet_NaN();
    }
    else
    {
      if(ef::isone(x))
      {
        return ef::riemann_zeta(n);
      }
      else
      {
        static const double d1 = static_cast<double>(1.0);
        static const double d2 = static_cast<double>(2.0);

        const double xd = ef::to_double(x);

        const bool bo_converge_small = HypergeometricUtil::AsympConverge(std::deque<double>(static_cast<std::size_t>(n + 1), d1),
                                                                         std::deque<double>(static_cast<std::size_t>(n),     d2),
                                                                         xd);

        return bo_converge_small ? PolyLog_Series::AtZeroForPositiveN(n, x)
                                 : PolyLog_Series::AtOneForPositiveN<e_float>(n, x);
      }
    }
  }
}
