
// *****************************************************************************
// Filename    : e_float_math.cpp
// 
// Project     : Multiple precision mathematics
// 
// Date        : 28.02.2004
// 
// Description : Real variable mathematics for e_float.
// 
// *****************************************************************************

#include <e_float/e_float.h>
#include <functions/constants/constants.h>
#include <functions/elementary/elementary.h>

e_float ef::sin(const e_float& x)
{
  if(!ef::isfinite(x))
  {
    return std::numeric_limits<e_float>::quiet_NaN();
  }

  if(x.has_its_own_sin())
  {
    return e_float::my_sin(x);
  }

  // Local copy of the argument
  e_float xx = x;

  // Analyze and prepare the phase of the argument.
  // Make a local, positive copy of the argument, xx.
  // The argument xx will be reduced to 0 <= xx <= pi/2.
  bool b_negate_sin = false;

  if(ef::isneg(x))
  {
    xx = -xx;
    b_negate_sin = !b_negate_sin;
  }

  // Remove even multiples of pi.
  if(xx > ef::pi())
  {
    e_float n_pi = ef::integer_part(xx / ef::pi());
    xx -= n_pi * ef::pi();

    // Adjust signs if the multiple of pi is not even.
    const bool b_n_pi_is_even = ef::iszero(ef::decimal_part(n_pi / static_cast<INT32>(2)));

    if(!b_n_pi_is_even)
    {
      b_negate_sin = !b_negate_sin;
    }
  }

  // Reduce the argument to 0 <= xx <= pi/2.
  if(xx > ef::pi_half())
  {
    xx = ef::pi() - xx;
  }
  
  const bool b_zero    = ef::iszero(xx);
  const bool b_pi_half = ef::iszero(xx - ef::pi_half());

  // Check if the reduced argument is very close to 0 or pi/2.
  const bool    b_near_zero    = ef::small_arg(xx);
  const e_float delta_pi_half  = ef::pi_half() - xx;
  const bool    b_near_pi_half = ef::small_arg(delta_pi_half);
  
  e_float sin_val;

  if(b_zero)
  {
    sin_val = ef::zero();
  }
  else if(b_pi_half)
  {
    sin_val = ef::one();
  }
  else if(b_near_zero)
  {
    const e_float x_squared = xx * xx;

    sin_val = xx * ef::hyp0F1(ef::three_half(), -x_squared / static_cast<INT32>(4));
  }
  else if(b_near_pi_half)
  {
    sin_val = ef::hyp0F1(ef::half(), -(delta_pi_half * delta_pi_half) / static_cast<INT32>(4));
  }
  else
  {
    // Scale to a small argument for an efficient Taylor series,
    // implemented as a hypergeometric function. Use a standard
    // divide by three identity a certain number of times.
    // Here we use division by 3^9 --> (19683 = 3^9).

    const bool b_scale = xx.order() > static_cast<INT64>(-4);

    static const INT32 n_scale           = static_cast<INT32>(9);
    static const INT32 n_three_pow_scale = static_cast<INT32>(static_cast<INT64>(::pow(static_cast<double>(3.0), static_cast<double>(n_scale)) + static_cast<double>(0.5)));
    
    if(b_scale)
    {
      xx /= n_three_pow_scale;
    }

    // Now with small arguments, we are ready for a series expansion.
    sin_val = xx * ef::hyp0F1(ef::three_half(), -(xx * xx) / static_cast<INT32>(4));

    // Convert back using multiple angle identity.
    if(b_scale)
    {
      for(INT32 k = static_cast<INT32>(0); k < n_scale; k++)
      {
        // Rescale the cosine value using the multiple angle identity.
        sin_val  = - ((sin_val * (sin_val * sin_val)) * static_cast<INT32>(4))
                   +  (sin_val * static_cast<INT32>(3));
      }
    }
  }

  return !b_negate_sin ? sin_val : -sin_val;
}

e_float ef::cos(const e_float& x)
{
  if(!ef::isfinite(x))
  {
    return std::numeric_limits<e_float>::quiet_NaN();
  }

  if(x.has_its_own_cos())
  {
    return e_float::my_cos(x);
  }

  // Local copy of the argument
  e_float xx = x;

  // Analyze and prepare the phase of the argument.
  // Make a local, positive copy of the argument, xx.
  // The argument xx will be reduced to 0 <= xx <= pi/2.
  bool b_negate_cos = false;

  if(ef::isneg(x))
  {
    xx = -xx;
  }

  // Remove even multiples of pi.
  if(xx > ef::pi())
  {
    e_float n_pi = ef::integer_part(xx / ef::pi());
    xx -= n_pi * ef::pi();

    // Adjust signs if the multiple of pi is not even.
    const bool b_n_pi_is_even = ef::iszero(ef::decimal_part(n_pi / static_cast<INT32>(2)));

    if(!b_n_pi_is_even)
    {
      b_negate_cos = !b_negate_cos;
    }
  }

  // Reduce the argument to 0 <= xx <= pi/2.
  if(xx > ef::pi_half())
  {
    xx = ef::pi() - xx;
    b_negate_cos = !b_negate_cos;
  }
  
  const bool b_zero    = ef::iszero(xx);
  const bool b_pi_half = ef::iszero(xx - ef::pi_half());

  // Check if the reduced argument is very close to 0 or pi/2.
  const bool    b_near_zero    = ef::small_arg(xx);
  const e_float delta_pi_half  = ef::pi_half() - xx;
  const bool    b_near_pi_half = ef::small_arg(delta_pi_half);
  
  e_float cos_val;

  if(b_zero)
  {
    cos_val = ef::one();
  }
  else if(b_pi_half)
  {
    cos_val = ef::zero();
  }
  else if(b_near_zero)
  {
    const e_float x_squared = xx * xx;

    cos_val = ef::hyp0F1(ef::half(), -x_squared / static_cast<INT32>(4));
  }
  else if(b_near_pi_half)
  {
    cos_val = delta_pi_half * ef::hyp0F1(ef::three_half(), -(delta_pi_half * delta_pi_half) / static_cast<INT32>(4));
  }
  else
  {
    // Scale to a small argument for an efficient Taylor series,
    // implemented as a hypergeometric function. Use a standard
    // divide by three identity a certain number of times.
    // Here we use division by 3^9 --> (19683 = 3^9).

    const bool b_scale = xx.order() > static_cast<INT64>(-4);

    static const INT32 n_scale           = static_cast<INT32>(9);
    static const INT32 n_three_pow_scale = static_cast<INT32>(static_cast<INT64>(::pow(static_cast<double>(3.0), static_cast<double>(n_scale)) + static_cast<double>(0.5)));
    
    if(b_scale)
    {
      xx /= n_three_pow_scale;
    }

    // Now with small arguments, we are ready for a series expansion.
    cos_val = ef::hyp0F1(ef::half(), -(xx * xx) / static_cast<INT32>(4));

    // Convert back using multiple angle identity.
    if(b_scale)
    {
      for(INT32 k = static_cast<INT32>(0); k < n_scale; k++)
      {
        // Rescale the cosine value using the multiple angle identity.
        cos_val  =   ((cos_val * (cos_val * cos_val)) * static_cast<INT32>(4))
                   -  (cos_val * static_cast<INT32>(3));
      }
    }
  }

  return !b_negate_cos ? cos_val : -cos_val;
}

void ef::sincos(const e_float& x, e_float* const p_sin, e_float* const p_cos)
{
  if(p_sin != static_cast<e_float*>(0u)) { *p_sin = ef::sin(x); }
  if(p_cos != static_cast<e_float*>(0u)) { *p_cos = ef::cos(x); }
}

e_float ef::tan(const e_float& x)
{
  if(x.has_its_own_tan())
  {
    return e_float::my_tan(x);
  }
  else
  {
    return ef::sin(x) / ef::cos(x);
  }
}

e_float ef::csc(const e_float& x) { return ef::one()  / ef::sin(x); }
e_float ef::sec(const e_float& x) { return ef::one()  / ef::cos(x); }
e_float ef::cot(const e_float& x) { return ef::cos(x) / ef::sin(x); }

e_float ef::asin(const e_float& x)
{
  if(!ef::isfinite(x))
  {
    return std::numeric_limits<e_float>::quiet_NaN();
  }

  if(x.has_its_own_asin())
  {
    return e_float::my_asin(x);
  }

  const bool b_neg = ef::isneg(x);

  const e_float xx = !b_neg ? x : -x;

  if(xx > ef::one())
  {
    return std::numeric_limits<e_float>::quiet_NaN();
  }

  if(ef::iszero(x))
  {
    return ef::zero();
  }
  
  if(ef::isone(xx))
  {
    return !b_neg ? ef::pi_half() : -ef::pi_half();
  }

  if(ef::small_arg(xx))
  {
    // http://functions.wolfram.com/ElementaryFunctions/ArcSin/26/01/01/
    const e_float asin_value = x * ef::hyp2F1(ef::half(),
                                              ef::half(),
                                              ef::three_half(),
                                              (x * x));

    return !b_neg ? asin_value : -asin_value;
  }
  else if(ef::near_one(xx))
  {
    const e_float dx1 = ef::one() - xx;

    const e_float asin_value =    ef::pi_half()
                               - (  ef::sqrt(dx1 * static_cast<INT32>(2))
                                  * ef::hyp2F1(ef::half(),
                                               ef::half(),
                                               ef::three_half(),
                                               dx1 / static_cast<INT32>(2)));

    return !b_neg ? asin_value : -asin_value;
  }

  // Get initial estimate using standard math function asin.
  double dd;
  INT64  ne;
  ef::to_parts(xx, dd, ne);

  static const INT64 p10_min = static_cast<INT64>(std::numeric_limits<double>::min_exponent10);
  static const INT64 p10_max = static_cast<INT64>(std::numeric_limits<double>::max_exponent10);

  const double de = static_cast<double>(ne < static_cast<INT64>(0) ? static_cast<INT32>(std::max(ne, p10_min))
                                                                   : static_cast<INT32>(std::min(ne, p10_max)));

  e_float value = e_float(::asin(dd * ::pow(static_cast<double>(10.0), de)));

  // Newton-Raphson iteration

  static const INT32 double_digits10_minus_one = static_cast<INT32>(static_cast<INT32>(std::numeric_limits<double>::digits10) - static_cast<INT32>(1));

  for(INT32 digits = double_digits10_minus_one; digits <= static_cast<INT32>(ef::tol()); digits *= static_cast<INT32>(2))
  {
    e_float s, c;
    ef::sincos(value, &s, &c);
    value -= (s - xx) / c;
  }

  return !b_neg ? value : -value;
}

e_float ef::acos(const e_float& x)
{
  if(!ef::isfinite(x))
  {
    return std::numeric_limits<e_float>::quiet_NaN();
  }

  if(x.has_its_own_acos())
  {
    return e_float::my_acos(x);
  }

  if(ef::fabs(x) > ef::one()) { return std::numeric_limits<e_float>::quiet_NaN(); }

  return ef::iszero(x) ? ef::pi_half() : ef::pi_half() - ef::asin(x);
}

namespace Atan_Series
{
  static e_float AtZero    (const e_float& x);
  static e_float AtInfinity(const e_float& x);
}

static e_float Atan_Series::AtZero(const e_float& x)
{
  // http://functions.wolfram.com/ElementaryFunctions/ArcTan/26/01/01/
  return x * ef::hyp2F1( ef::one(),
                         ef::half(),
                         ef::three_half(),
                        -(x * x));
}

static e_float Atan_Series::AtInfinity(const e_float& x)
{
  // http://functions.wolfram.com/ElementaryFunctions/ArcTan/26/01/01/
  return ef::pi_half() - ef::hyp2F1( ef::half(),
                                     ef::one(),
                                     ef::three_half(),
                                    -ef::one() / (x * x)) / x;
}

e_float ef::atan(const e_float& x)
{
  if(!ef::isfinite(x))
  {
    return std::numeric_limits<e_float>::quiet_NaN();
  }

  if(x.has_its_own_atan())
  {
    return e_float::my_atan(x);
  }

  const INT64 order = x.order();

  if(x.isinf() || order > ef::tol())
  {
    return ef::ispos(x) ? ef::pi_half() : -ef::pi_half();
  }
  else if(ef::iszero(x))
  {
    return ef::zero();
  }
  
  if(ef::small_arg(x))
  {
    return Atan_Series::AtZero(x);
  }

  if(ef::large_arg(x))
  {
    return Atan_Series::AtInfinity(x);
  }

  const bool b_neg = ef::isneg(x);
  
  const e_float xx = !b_neg ? x : -x;
  
  // Get initial estimate using standard math function atan or a series
  // expansion for rather large arguments having order 3 or larger.
  double dd;
  INT64  ne;
  ef::to_parts(xx, dd, ne);

  static const INT64 p10_min = static_cast<INT64>(std::numeric_limits<double>::min_exponent10);
  static const INT64 p10_max = static_cast<INT64>(std::numeric_limits<double>::max_exponent10);

  const double de = static_cast<double>(ne < static_cast<INT64>(0) ? static_cast<INT32>(std::max(ne, p10_min))
                                                                   : static_cast<INT32>(std::min(ne, p10_max)));

  e_float value = order < static_cast<INT64>(2) ? e_float(::atan(dd * ::pow(static_cast<double>(10.0), de)))
                                                : Atan_Series::AtInfinity(xx);

  // Newton-Raphson iteration
  static const INT32 double_digits10_minus_one = static_cast<INT32>(static_cast<INT32>(std::numeric_limits<double>::digits10) - static_cast<INT32>(1));

  for(INT32 digits = double_digits10_minus_one; digits <= static_cast<INT32>(ef::tol()); digits *= static_cast<INT32>(2))
  {
    e_float s, c;
    ef::sincos(value, &s, &c);
    value += c * ((xx * c) - s);
  }

  return !b_neg ? value : -value;
}

e_float ef::atan2(const e_float& y, const e_float& x)
{
  if(!ef::isfinite(x) || !ef::isfinite(y))
  {
    return x;
  }

  // y is zero
  if(ef::iszero(y))
  {
    return ef::isneg(x) ? ef::pi() : ef::zero();
  }

  // x is zero
  if(ef::iszero(x))
  {
    return ef::sgn(y) * ef::pi_half();
  }
  
  // y is infinite
  if(y.isinf())
  {
    return (ef::ispos(x) && ef::ispos(y)) || (ef::isneg(x) && ef::isneg(y))
             ?  std::numeric_limits<e_float>::infinity()
             : -std::numeric_limits<e_float>::infinity();
  }

  // Compute atan(y / x), ignoring sign.
  e_float atan_term(ef::atan(y / x));

  // Determine quadrant (sign) based on signs of x, y
  const bool y_neg = ef::isneg(y);
  const bool x_neg = ef::isneg(x);
  
  if(y_neg == x_neg)
  {
    // Both negative or both positive
    return x_neg ? atan_term - ef::pi() : atan_term;
  }
  else
  {
    // Different signs of x, y
    return x_neg ? atan_term + ef::pi() : atan_term;
  }
}
