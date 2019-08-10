
#include <vector>
#include <numeric>
#include <map>

#include <e_float/e_float.h>
#include <functions/constants/constants.h>
#include <functions/elementary/elementary.h>
#include <utility/util_lexical_cast.h>
#include <utility/util_power_x_pow_n.h>

namespace e_floatUtil
{
  static e_float rootn_inv(const e_float& x, const INT32 p);
}

e_float ef::pown(const e_float& x, const INT64 p)
{
  return Util::x_pow_n_template<e_float>(x, p);
}

e_float ef::pow2(const INT64 p)
{
  // Compute two raised to the power of p, in other words 2^p.
  switch(p)
  {
    case static_cast<INT64>(-3):
      return ef::eighth();
    case static_cast<INT64>(-2):
      return ef::quarter();
    case static_cast<INT64>(-1):
      return ef::half();
    case static_cast<INT64>(0):
      return ef::one();
    case static_cast<INT64>(1):
      return ef::two();
    case static_cast<INT64>(2):
      return ef::four();
    case static_cast<INT64>(3):
      return ef::eight();
    default:
      break;
  }

  if(p < static_cast<INT64>(0))
  {
    return ef::pow2(static_cast<INT64>(-p)).calculate_inv();
  }
  else if(p < static_cast<INT64>(std::numeric_limits<UINT64>::digits))
  {
    const UINT64 p2 = static_cast<UINT64>(static_cast<UINT64>(1uLL) << p);

    return e_float(p2);
  }
  else
  {
    return Util::x_pow_n_template(ef::two(), p);
  }
}

static e_float e_floatUtil::rootn_inv(const e_float& x, const INT32 p)
{
  // Compute the value of [1 / (rootn of x)] with n = p.

  // Generate the initial estimate using 1 / rootn.
  // Extract the mantissa and exponent for a "manual"
  // computation of the estimate.
  double dd;
  INT64  ne;
  ef::to_parts(x, dd, ne);

  // Adjust exponent and mantissa such that ne is an even power of p.
  while(ne % static_cast<INT64>(p))
  {
    ++ne;
    dd /= static_cast<double>(10.0);
  }
  
  // Estimate the one over the root using simple manipulations.
  const double one_over_rtn_d = ::pow(dd, static_cast<double>(-1.0) / static_cast<double>(p));

  // Set the result equal to the initial guess.
  e_float result(one_over_rtn_d, static_cast<INT64>(-ne / p));

  static const INT32 double_digits10_minus_one = static_cast<INT32>(static_cast<INT32>(std::numeric_limits<double>::digits10) - static_cast<INT32>(1));

  for(INT32 digits = double_digits10_minus_one; digits <= static_cast<INT32>(ef::tol()); digits *= static_cast<INT32>(2))
  {
    // Adjust precision of the terms.
    result.precision(static_cast<INT32>(digits * static_cast<INT32>(2)));

    // Next iteration
    e_float term = (((-ef::pown(result, p) * x) + ef::one()) / p) + ef::one();

    term.precision(static_cast<INT32>(digits * static_cast<INT32>(2)));

    result *= term;
  
  }

  result.precision(static_cast<INT32>(ef::tol()));

  return result;
}

e_float ef::inv (const e_float& x) { return e_float(x).calculate_inv(); }
e_float ef::sqrt(const e_float& x) { return e_float(x).calculate_sqrt(); }

e_float ef::cbrt(const e_float& x)
{
  return ef::rootn(x, static_cast<INT32>(3));
}

e_float ef::rootn(const e_float& x, const INT32 p)
{
  if(!ef::isfinite(x))
  {
    return std::numeric_limits<e_float>::quiet_NaN();
  }

  if(p < static_cast<INT32>(0))
  {
    return ef::rootn(ef::one() / x, static_cast<INT32>(-p));
  }

  const bool b_x_is_neg  = ef::isneg(x);

  if((p == static_cast<INT32>(0)) || b_x_is_neg)
  {
    return std::numeric_limits<e_float>::quiet_NaN();
  }
  else if(p == static_cast<INT32>(1))
  {
    return x;
  }
  else if(p == static_cast<INT32>(2))
  {
    return ef::sqrt(x);
  }
  else if((p == static_cast<INT32>(3)) && x.has_its_own_cbrt())
  {
    return e_float::my_cbrt(x);
  }

  const e_float xx = ef::fabs(x);

  const e_float rtn =  (x.has_its_own_rootn() ? e_float::my_rootn(xx, p)
                                             : e_floatUtil::rootn_inv(xx, p).calculate_inv());

  return (!b_x_is_neg ? rtn : -rtn);
}

e_float ef::exp(const e_float& x)
{
  if(x.has_its_own_exp())
  {
    return e_float::my_exp(x);
  }

  // Handle special arguments.
  if(ef::isnan(x))
  {
    return std::numeric_limits<e_float>::quiet_NaN();
  }

  if(ef::isinf(x))
  {
    return (!ef::isneg(x) ? std::numeric_limits<e_float>::infinity() : ef::zero());
  }

  if(ef::iszero(x) || x.order() < -ef::tol())
  {
    return ef::one();
  }

  // Get local copy of argument and force it to be positive.
  const bool bo_x_is_neg = ef::isneg(x);
  
  const e_float xx = !bo_x_is_neg ? x : -x;

  // Check the range of the argument.
  static const e_float maximum_arg_for_exp = e_float(std::numeric_limits<e_float>::max_exponent);

  if(xx > maximum_arg_for_exp)
  {
    // Overflow / underflow
    return !bo_x_is_neg ? std::numeric_limits<e_float>::infinity() : ef::zero();
  }

  // Check for pure-integer arguments which can be either signed or unsigned.
  if(ef::isint(x))
  {
    return ef::pown(ef::exp1(), ef::to_int64(x));
  }

  // The algorithm for exp has been taken from MPFUN.
  // exp(t) = [ (1 + r + r^2/2! + r^3/3! + r^4/4! ...)^p2 ] * 2^n
  // where p2 is a power of 2 such as 512, r = t_prime / p2, and
  // t_prime = t - n*ln2, with n chosen to minimize the absolute
  // value of t_prime. In the resulting Taylor series, which is
  // implemented as a hypergeometric function, |r| is bounded by
  // ln2 / p2. For small arguments, no scaling is done.

  static const e_float one_over_ln2 = ef::one() / ef::ln2();

  e_float nf = ef::integer_part(xx * one_over_ln2);
  
  // Scaling.
  static const INT32 p2 = static_cast<INT32>(UINT32(1u) << 11);

  const bool b_scale = xx.order() > static_cast<INT64>(-4);

  // Compute the exponential series of the (possibly) scaled argument.
  e_float exp_series = ef::hyp0F0(b_scale ? (xx - nf * ef::ln2()) / p2 : xx);

  if(b_scale)
  {
    exp_series  = ef::pown(exp_series, p2);
    exp_series *= ef::pow2(ef::to_int64(nf));
  }
  
  return !bo_x_is_neg ? exp_series : exp_series.calculate_inv();
}

namespace Log_Series
{
  static e_float AtOne(const e_float& x)
  {
    // This subroutine computes the series representation of Log[1 + x]
    // for small x without losing precision.

    // http://functions.wolfram.com/ElementaryFunctions/Log/26/01/01/

    return x * ef::hyp2F1( ef::one(), ef::one(), ef::two(), -x);
  }
}

e_float ef::log(const e_float& x)
{
  // Handle special arguments.
  if(ef::isnan(x) || ef::isneg(x) || ef::isinf(x))
  {
    return std::numeric_limits<e_float>::quiet_NaN();
  }

  if(x.has_its_own_log())
  {
    return e_float::my_log(x);
  }

  if(ef::iszero(x))
  {
    return -std::numeric_limits<e_float>::infinity();
  }

  if(ef::isone(x))
  {
    return ef::zero();
  }

  // Make a local copy
  e_float xx = x;

  // Compute the delta of the argument compared to one.
  const e_float x_minus_one  = xx - ef::one();

  if(ef::near_one(xx))
  {
    return Log_Series::AtOne(x_minus_one);
  }

  // For large arguments, the value will be broken into two parts
  // in order to facilitate the convergence of the Newton iteration.
  const bool b_correction = (   (xx.order() > static_cast<INT64>(+1000))
                             || (xx.order() < static_cast<INT64>(-1000)));

  e_float correction;

  if(b_correction)
  {
    // The argument xx is of the form a * 10^b.
    // It will be broken into two parts: log(a) + b * log(10).
    const bool b_neg_exp = xx.order() < static_cast<INT64>(0);

    // Remove a large power of ten from the argument. But be sure to leave the argument
    // large enough (or small enough) to avoid entering the near-one range.
    const INT64 n_order   = xx.order();
    const INT64 n_exp     = !b_neg_exp ? n_order : -n_order;
    const INT64 delta_exp = static_cast<INT64>(n_exp - static_cast<INT64>(8));

    // Convert the scaling power of ten to a string and subsequently to an e_float.
    const e_float ef_delta_exp("1E" + Util::lexical_cast(delta_exp));

    !b_neg_exp ? xx /= ef_delta_exp : xx *= ef_delta_exp;

    correction  = ef::ln10() * e_float(delta_exp);

    if(b_neg_exp)
    {
      correction = -correction;
    }
  }

  // Generate the initial estimate using double precision log combined with
  // the exponent for a "manual" computation of the initial iteration estimate.

  static const double lg10_d = ::log(static_cast<double>(10.0));

  static const INT64 n32_min = static_cast<INT64>(std::numeric_limits<INT32>::min());
  static const INT64 n32_max = static_cast<INT64>(std::numeric_limits<INT32>::max());

  // computation of the estimate.
  double dd;
  INT64  ne;
  ef::to_parts(xx, dd, ne);

  const double nd = static_cast<double>(ne < static_cast<INT64>(0) ? static_cast<INT32>(std::max(ne, n32_min))
                                                                   : static_cast<INT32>(std::min(ne, n32_max)));

  const double dlog = ::log(dd) + (nd * lg10_d);

  const double d10  = !ef::iszero(dlog) ? (::log10(::fabs(dlog)) + static_cast<double>(0.5))
                                        : static_cast<double>(0.0);

  const INT64  p10  =  ef::ispos(dlog)  ? static_cast<INT64>(d10) : static_cast<INT64>(-d10);

  e_float log_val   = !ef::iszero(dlog) ? e_float(dlog / ::pow(static_cast<double>(10.0), static_cast<double>(static_cast<INT32>(p10))), p10)
                      : x_minus_one;

  // Newton-Raphson iteration
  static const INT32 double_digits10_minus_one = static_cast<INT32>(static_cast<INT32>(std::numeric_limits<double>::digits10) - static_cast<INT32>(1));

  for(INT32 digits = double_digits10_minus_one; digits <= static_cast<INT32>(ef::tol()); digits *= static_cast<INT32>(2))
  {
    // Adjust precision of the terms.
    log_val.precision(static_cast<INT32>(digits * static_cast<INT32>(2)));
         xx.precision(static_cast<INT32>(digits * static_cast<INT32>(2)));

    const e_float exp_minus_log = ef::exp(-log_val);

    log_val += (xx * exp_minus_log) - ef::one();
  }

  return !b_correction ? log_val : log_val + correction;
}

e_float ef::log10(const e_float& x)                   { return ef::log(x) / ef::ln10(); }
e_float ef::loga (const e_float& a, const e_float& x) { return ef::log(x) / ef::log(a); }
e_float ef::log1p(const e_float& x)                   { return Log_Series::AtOne(x); }

e_float ef::log1p1m2(const e_float& x)
{
  // This subroutine calculates the series representation of (1/2) Log[(1 + x) / (1 - x)]
  // for small x without losing precision.

  if(!ef::isfinite(x))
  {
    return x;
  }

  if((x <= ef::one_minus()) || (x >= ef::one()))
  {
    return std::numeric_limits<e_float>::quiet_NaN();
  }

  // for values of x near one.
  const e_float x2 = x * x;
        e_float xn = x;

  e_float sum = xn;

  // Series representation of (1/2) Log[(1 + x) / (1 - x)] as given in
  // Schaum's Outlines: Mathematical Handbook of Formulas and Tables,
  // Second Edition, equation 22.8, page 136.
  for(INT32 n = static_cast<INT32>(3); n < ef::max_iteration(); n += static_cast<INT32>(2))
  {
    xn *= x2;

    const e_float term = xn / n;

    const INT64 order_check = static_cast<INT64>(term.order() - sum.order());

    if((n > static_cast<INT32>(20)) && (order_check < -ef::tol()))
    {
      break;
    }

    sum += term;
  }

  return sum;
}

e_float ef::pow(const e_float& x, const e_float& a)
{
  if(!ef::isfinite(x) || ef::isone(a))
  {
    return x;
  }

  if(ef::iszero(x))
  {
    return ef::one();
  }
  
  if(ef::iszero(a))
  {
    return ef::one();
  }
  
  const bool bo_a_isint = ef::isint(a);

  if(ef::isneg(x) && !bo_a_isint)
  {
    return std::numeric_limits<e_float>::quiet_NaN();
  }

  if(a <= ef::one_minus())
  {
    return ef::one() / ef::pow(x, -a);
  }

  const e_float a_int = ef::integer_part(a);
  const INT64   an    = ef::to_int64(a_int);
  const e_float da    = a - a_int;

  if(bo_a_isint)
  {
    return ef::pown(x, an);
  }

  static const e_float nine_tenths = ef::nine() / static_cast<INT32>(10);

  if(ef::ispos(x) && (x > ef::tenth()) && (x < nine_tenths))
  {
    if(ef::small_arg(a))
    {
      // Series expansion for small a.
      return ef::hyp0F0(a * ef::log(x));
    }
    else
    {
      // Series expansion for moderately sized x. Note that for large power of a,
      // the power of the integer part of a is calculated using the pown function.
      return an != static_cast<INT64>(0) ? ef::hyp1F0(-da, ef::one() - x) * ef::pown(x, an)
                                         : ef::hyp1F0( -a, ef::one() - x);
    }
  }
  else
  {
    // Series expansion for pow(x, a). Note that for large power of a, the power
    // of the integer part of a is calculated using the pown function.
    return an ? ef::exp(da * ef::log(x)) * ef::pown(x, an)
              : ef::exp( a * ef::log(x));
  }
}

e_float ef::sinh(const e_float& x)
{
  if(x.has_its_own_sinh())
  {
    return e_float::my_sinh(x);
  }

  if(!ef::isfinite(x))
  {
    return x;
  }

  e_float s;
  ef::sinhcosh(x, &s, static_cast<e_float*>(0u));
  return s;
}

e_float ef::cosh(const e_float& x)
{
  if(x.has_its_own_cosh())
  {
    return e_float::my_cosh(x);
  }

  if(!ef::isfinite(x))
  {
    return x;
  }

  e_float c;
  ef::sinhcosh(x, static_cast<e_float*>(0u), &c);
  return c;
}

void ef::sinhcosh(const e_float& x, e_float* const p_sinh, e_float* const p_cosh)
{
  if(!ef::isfinite(x) || (!p_sinh && !p_cosh))
  {
    return;
  }
  
  if(ef::iszero(x))
  {
    if(p_sinh)
    {
      *p_sinh = ef::zero();
    }

    if(p_cosh)
    {
      *p_cosh = ef::one();
    }
    
    return;
  }

  const e_float e_px = ef::exp(x);
  const e_float e_mx = ef::one() / e_px;

  if(p_sinh)
  {
    *p_sinh  = e_px;
    *p_sinh -= e_mx;
    p_sinh->div_by_int(2);
  }

  if(p_cosh)
  {
    *p_cosh  = e_px;
    *p_cosh += e_mx;
    p_cosh->div_by_int(2);
  }
}

e_float ef::tanh(const e_float& x)
{
  if(x.has_its_own_tanh())
  {
    return e_float::my_tanh(x);
  }

  e_float c, s;
  ef::sinhcosh(x, &s, &c);
  return s * c.calculate_inv();
}

e_float ef::asinh(const e_float& x)
{
  if(x.has_its_own_asinh())
  {
    return e_float::my_asinh(x);
  }

  if(!ef::isfinite(x))
  {
    return std::numeric_limits<e_float>::quiet_NaN();
  }

  if(ef::iszero(x))
  {
    return ef::zero();
  }
  else
  {
    const e_float value = ef::log(ef::fabs(x) + ef::sqrt((x * x) + ef::one()));

    return !ef::isneg(x) ? value : -value;
  }
}

e_float ef::acosh(const e_float& x)
{
  if(x.has_its_own_acosh())
  {
    return e_float::my_acosh(x);
  }

  if(!ef::isfinite(x))
  {
    return std::numeric_limits<e_float>::quiet_NaN();
  }

  if(ef::isneg(x) || x < ef::one())
  {
    return std::numeric_limits<e_float>::quiet_NaN();
  }
  
  if(ef::isone(x))
  {
    return ef::one();
  }

  const e_float x_minus_one = x - ef::one();

  if(ef::small_arg(x_minus_one))
  {
    return   (ef::sqrt2() * ef::sqrt(x_minus_one))
           *  ef::hyp2F1( ef::half(),
                          ef::half(),
                          ef::three_half(),
                         -x_minus_one / static_cast<INT32>(2));
  }
  else
  {
    return ef::log(x + ef::sqrt((x * x) - ef::one()));
  }
}

e_float ef::atanh(const e_float& x)
{
  if(x.has_its_own_atanh())
  {
    return e_float::my_atanh(x);
  }

  if(!ef::isfinite(x))
  {
    return x;
  }

  const e_float xx = ef::fabs(x);
  
  if(xx >= ef::one())
  {
    return std::numeric_limits<e_float>::quiet_NaN();
  }

  const e_float value = (ef::small_arg(x) ?  ef::log1p1m2(x)
                                          : (ef::log((ef::one() + x) / (ef::one() - x)) / static_cast<INT32>(2)));

  return (!ef::isneg(xx) ? value : -value);
}
