
#include <functions/constants/constants.h>
#include <functions/elementary/elementary.h>
#include <functions/gamma/gamma.h>
#include <functions/hypergeometric/hypergeometric.h>
#include <functions/hypergeometric/hypergeometric_util.h>
#include <functions/hypergeometric/parabolic_cylinder.h>
#include <functions/hypergeometric/parabolic_cylinder_d_asymp.h>
#include <functions/integer/integer.h>
#include <utility/util_digit_scale.h>
#include <utility/util_interpolate.h>


namespace ParabolicCylinderD_Series
{
  static INT32 TransitionBandLower(const double& xd);
  static INT32 TransitionBandUpper(const double& xd);

  static const INT32& AsymptoticLimit(void);

  static e_float Reflection                     (const e_float& v, const e_float& x);
  static e_float AtZero                         (const e_float& v, const e_float& x);
  static e_float AtTransition                   (const e_float& v, const e_float& x);
  static e_float AtTransitionRecurUpForPositiveV(const e_float& v, const e_float& x);
  static e_float AtTransitionRecurUpForNegativeV(const e_float& v, const e_float& x);
  static e_float AtInfinity                     (const e_float& v, const e_float& x);
  static e_float AtInfinityRecurUpForPositiveV  (const e_float& v, const e_float& x);
  static e_float RecurUp_V_PositiveInteger      (const INT32 n,    const e_float& x);
}

static INT32 ParabolicCylinderD_Series::TransitionBandLower(const double& xd)
{
  static const std::tr1::array<Util::point<double>, static_cast<std::size_t>(6u)> lo_data =
  {{
    Util::point<double>(static_cast<double>( 4.0), static_cast<double>(-6000.0)),
    Util::point<double>(static_cast<double>( 8.0), static_cast<double>(-4000.0)),
    Util::point<double>(static_cast<double>(12.0), static_cast<double>(-2000.0)),
    Util::point<double>(static_cast<double>(20.0), static_cast<double>(-1000.0)),
    Util::point<double>(static_cast<double>(30.0), static_cast<double>(-2000.0)),
    Util::point<double>(static_cast<double>(40.0), static_cast<double>(-8000.0)),
  }};

  static const std::vector<Util::point<double> > lo(lo_data.begin(), lo_data.end());

  const double y = Util::linear_interpolate<double>::interpolate(xd, lo);

  return static_cast<INT32>(static_cast<INT64>(y * Util::DigitScale()));
}

static INT32 ParabolicCylinderD_Series::TransitionBandUpper(const double& xd)
{
  static const std::tr1::array<Util::point<double>, static_cast<std::size_t>(6u)> lo_data =
  {{
    Util::point<double>(static_cast<double>( 4.0), static_cast<double>(4000.0)),
    Util::point<double>(static_cast<double>( 8.0), static_cast<double>(3000.0)),
    Util::point<double>(static_cast<double>(12.0), static_cast<double>(1000.0)),
    Util::point<double>(static_cast<double>(20.0), static_cast<double>( 500.0)),
    Util::point<double>(static_cast<double>(30.0), static_cast<double>(1000.0)),
    Util::point<double>(static_cast<double>(40.0), static_cast<double>(6000.0)),
  }};

  static const std::vector<Util::point<double> > lo(lo_data.begin(), lo_data.end());

  const double y = Util::linear_interpolate<double>::interpolate(xd, lo);

  return static_cast<INT32>(static_cast<INT64>(y * Util::DigitScale()));
}

static const INT32& ParabolicCylinderD_Series::AsymptoticLimit(void)
{
  static const INT32 asy_lim = static_cast<INT32>(static_cast<INT64>(static_cast<double>(-50000.0) * Util::DigitScale()));
  return asy_lim;
}

static e_float ParabolicCylinderD_Series::Reflection(const e_float& v, const e_float& x)
{
  const e_float minus_x             = -x;
  const e_float x_squared           = minus_x * minus_x;
  const e_float v_half              = v / static_cast<INT32>(2);
  const e_float v_plus_one_over_two = v_half + ef::half();
  const e_float factor_part1        = ((ef::pow(ef::two(), v_plus_one_over_two) * v) * (ef::gamma(v_half) * ef::sin(v_half * ef::pi()))) / ef::sqrt_pi();
  const e_float factor_part2        = minus_x * ef::exp(-x_squared / static_cast<INT32>(4));
  const e_float factor              = factor_part1 * factor_part2;

  return ef::weber_d(v, minus_x) - (factor * ef::hyperg_1f1(ef::half() - v_half, ef::three_half(), x_squared / static_cast<INT32>(2)));
}

static e_float ParabolicCylinderD_Series::AtZero(const e_float& v, const e_float& x)
{
  // Use the Taylor series expansion for |x| small.
  // http://functions.wolfram.com/HypergeometricFunctions/ParabolicCylinderD/06/01/02/

  const e_float v_half             = v / static_cast<INT32>(2);
  const e_float minus_v_half       = -v_half;
  const e_float half_minus_v_half  = ef::half() + minus_v_half;
  const e_float x_squared_over_two = (x * x) / static_cast<INT32>(2);

  const e_float term1 =                       ef::hyperg_1f1(minus_v_half,      ef::half(),       x_squared_over_two)  / ef::gamma(half_minus_v_half);
  const e_float term2 = ((ef::sqrt2() * x) * (ef::hyperg_1f1(half_minus_v_half, ef::three_half(), x_squared_over_two)) / ef::gamma(minus_v_half));

  const e_float factor = (ef::pow(ef::two(), v_half) * ef::sqrt_pi()) * ef::exp(-x_squared_over_two / static_cast<INT32>(2));

  return factor * (term1 - term2);
}

static e_float ParabolicCylinderD_Series::AtTransition(const e_float& v, const e_float& x)
{
  // Use the generalized hypergeometric series which is the difference of two
  // hypergeometric functions. This is the same as the Taylor series expansion
  // near zero.

  return ParabolicCylinderD_Series::AtZero(v, x);
}

static e_float ParabolicCylinderD_Series::AtTransitionRecurUpForPositiveV(const e_float& v, const e_float& x)
{
  // Select an order of parabolic cylinder function which is lower than the desired order
  // and use scaled, upward recursion (originally seeded with 0 and 1) to recur to an
  // order which is high enough to be computed to the desired precision. Then compute
  // the value of the parabolic cylinder function at the high end of the recursion in
  // order to obtain the scale factor.

  // Select an order for recursion which is lower than the order v. Since this subroutine
  // is used for moderate, positive values of v, this order should negative. Therefore
  // the recursion start can be computed assuming that it will be a negative value.
  // More recursions than expected seem to be needed to meet the precision requirements.
  // Therefore empirically calibrated upper and lower bands have been used to extend the
  // recursion start and end points.

  const double xd = ef::to_double(x);

  const INT32   n_v              = static_cast<INT32> (ef::to_int64(v));
  const INT32   n_top            = static_cast<INT32>(ParabolicCylinderD_Series::TransitionBandUpper(xd) + static_cast<INT32>(1));
  const INT32   n_minus          = ParabolicCylinderD_Series::TransitionBandLower(xd);
  const e_float v_frac           = ef::decimal_part(v);
  const e_float one_minus_v_frac = ef::one() - v_frac;

  e_float v_recur_minus_one = -one_minus_v_frac + n_minus;

  e_float Dv_m2 = ef::zero();
  e_float Dv_m1 = ef::one();
  e_float Dv;

  e_float Dvn;

  // Do the upward recursion.

  for(INT32 m = n_minus; m <= n_top; m++)
  {
    Dv = (x * Dv_m1) - (v_recur_minus_one * Dv_m2);

    if(m == n_v)
    {
      Dvn = Dv;
    }

    Dv_m2 = Dv_m1;
    Dv_m1 = Dv;

    ++v_recur_minus_one;
  }

  // Scale the result using the calculable value at the top of the recursion.

  const e_float scale = ef::weber_d(v_frac + n_top, x) / Dv;

  return Dvn * scale;
}

static e_float ParabolicCylinderD_Series::AtTransitionRecurUpForNegativeV(const e_float& v, const e_float& x)
{
  // Use an asymptotic expansion in combination with upward recursion.

  const INT32 n_v = static_cast<INT32>(ef::to_int64(v));

  if(n_v < ParabolicCylinderD_Series::AsymptoticLimit())
  {
    return ParabolicCylinder_Asymp::WeberU_TemmeSiamBook_Eq_12_119_Asymp(-(v + ef::half()), x);
  }

  const INT32 n_v_ofs = static_cast<INT32>(ParabolicCylinderD_Series::AsymptoticLimit() - n_v);
  const INT32 n_bot   = n_v + n_v_ofs;

  const e_float v_frac = ef::decimal_part(v);

  const e_float v_m2 = v_frac + (n_bot - static_cast<INT32>(2));
  const e_float v_m1 = v_frac + (n_bot - static_cast<INT32>(1));

  const e_float a_m2 = -(v_m2 + ef::half());
  const e_float a_m1 = -(v_m1 + ef::half());

  e_float Dv_m2 = ParabolicCylinder_Asymp::WeberU_TemmeSiamBook_Eq_12_119_Asymp(a_m2, x);
  e_float Dv_m1 = ParabolicCylinder_Asymp::WeberU_TemmeSiamBook_Eq_12_119_Asymp(a_m1, x);
  e_float Dv;

  e_float v_recur_minus_one = v_m1;

  for(INT32 m = n_bot; m <= n_v; m++)
  {
    Dv = (x * Dv_m1) - (v_recur_minus_one * Dv_m2);

    Dv_m2 = Dv_m1;
    Dv_m1 = Dv;

    ++v_recur_minus_one;
  }

  return Dv;
}

static e_float ParabolicCylinderD_Series::AtInfinity(const e_float& v, const e_float& x)
{
  // Use the hypergeometric representation for large positive values of x.
  // http://functions.wolfram.com/HypergeometricFunctions/ParabolicCylinderD/06/02/

  const e_float x_squared         = x * x;
  const e_float v_half            = v / static_cast<INT32>(2);
  const e_float minus_v_half      = -v_half;
  const e_float half_minus_v_half = ef::half() + minus_v_half;

  const e_float term                    = ef::hyperg_2f0(half_minus_v_half, minus_v_half, -ef::two() / x_squared);
  const e_float x_pow_v                 = ef::pow(x, v);
  const e_float exp_minus_xsq_over_four = ef::exp(-x_squared / static_cast<INT32>(4));
  const e_float factor                  = x_pow_v * exp_minus_xsq_over_four;

  return term * factor;
}

static e_float ParabolicCylinderD_Series::AtInfinityRecurUpForPositiveV(const e_float& v, const e_float& x)
{
  const INT32 Nm = static_cast<INT32>(ef::to_int64(v));

  const e_float v_frac = ef::decimal_part(v);

  if(Nm == static_cast<INT32>(0))
  {
    return ef::weber_d(v_frac, x);
  }
  else if(Nm == static_cast<INT32>(1))
  {
    return ef::weber_d(v_frac + ef::one(), x);
  }
  else
  {
    e_float v_minus_one = v_frac + ef::one();

    e_float Dv_m2 = ef::weber_d(v_frac, x);
    e_float Dv_m1 = ef::weber_d(v_minus_one, x);
    e_float Dv;

    for(INT32 m = 2; m <= Nm; m++)
    {
      Dv = (x * Dv_m1) - (v_minus_one * Dv_m2);

      Dv_m2 = Dv_m1;
      Dv_m1 = Dv;

      ++v_minus_one;
    }

    return Dv;
  }
}

static e_float ParabolicCylinderD_Series::RecurUp_V_PositiveInteger(const INT32 n, const e_float& x)
{
  const e_float minus_x_squared_over_four = -(x * x) / static_cast<INT32>(4);

  const e_float exp_minus_x_squared_over_four = ef::exp(minus_x_squared_over_four);

  e_float Dn_m2 =     exp_minus_x_squared_over_four;
  e_float Dn_m1 = x * exp_minus_x_squared_over_four;
  e_float Dn;

  if(n == static_cast<INT32>(0))
  {
    return Dn_m2;
  }
  else if(n == static_cast<INT32>(1))
  {
    return Dn_m1;
  }

  for(INT32 m = 2; m <= n; m++)
  {
    const INT32 n_minus_one = static_cast<INT32>(m - static_cast<INT32>(1));

    Dn = (x * Dn_m1) - (n_minus_one * Dn_m2);

    Dn_m2 = Dn_m1;
    Dn_m1 = Dn;
  }

  return Dn;
}

e_float ef::weber_d(const e_float& v, const e_float& x)
{
  const e_float v_half            = v / static_cast<INT32>(2);
  const e_float minus_v_half      = -v_half;
  const e_float half_minus_v_half = ef::half() + minus_v_half;

  if(ef::isint(v) && ef::ispos(v))
  {
    // Handle positive integer order.
    const INT32 n_v = static_cast<INT32>(ef::to_int64(v));

    const e_float Dn = ParabolicCylinderD_Series::RecurUp_V_PositiveInteger(n_v, ef::fabs(x));

    if(ef::isneg(x) && static_cast<INT32>(n_v % static_cast<INT32>(2)) != static_cast<INT32>(0))
    {
      return -Dn;
    }
    else
    {
      return Dn;
    }
  }

  if(ef::iszero(x))
  {
    return (ef::pow(ef::two(), v_half) * ef::sqrt_pi()) / ef::gamma(half_minus_v_half);
  }

  if(ef::isneg(x))
  {
    // Use reflection for negative argument.
    return ParabolicCylinderD_Series::Reflection(v, x);
  }

  static const e_float nine_halves = ef::four() + ef::half();

  const bool bo_converge_small = (x < nine_halves) && (ef::fabs(v) < ef::two_k());

  if(bo_converge_small)
  {
    return ParabolicCylinderD_Series::AtZero(v, x);
  }

  const double xd = ef::to_double(x);

  const std::tr1::array<double, 2u> p_data =
  {{
    ef::to_double(half_minus_v_half),
    ef::to_double(minus_v_half)
  }};

  const bool bo_converge_large = HypergeometricUtil::AsympConverge(std::deque<double>(p_data.begin(), p_data.end()),
                                                                   std::deque<double>(),
                                                                   static_cast<double>(-2.0) / (xd * xd));

  if(bo_converge_large)
  {
    return ParabolicCylinderD_Series::AtInfinity(v, x);
  }

  const INT32 n_v = static_cast<INT32>(ef::to_int64(v));

  if(n_v < ParabolicCylinderD_Series::AsymptoticLimit())
  {
    return ParabolicCylinder_Asymp::WeberU_TemmeSiamBook_Eq_12_119_Asymp(-(v + ef::half()), x);
  }

  if(ispos(v))
  {
    if(   (n_v > -ParabolicCylinderD_Series::AsymptoticLimit())
       && (x   >  ef::sqrt((v + ef::half()) * static_cast<INT32>(4)))
      )
    {
      return ParabolicCylinder_Asymp::WeberU_TemmeSiamBook_Eq_12_119_Asymp(-(v + ef::half()), x);
    }

    if(x > ef::fifty())
    {
      return ParabolicCylinderD_Series::AtInfinityRecurUpForPositiveV(v, x);
    }

    if(n_v < ParabolicCylinderD_Series::TransitionBandUpper(xd))
    {
      return ParabolicCylinderD_Series::AtTransitionRecurUpForPositiveV(v, x);
    }
    else
    {
      return ParabolicCylinderD_Series::AtTransition(v, x);
    }
  }
  else
  {
    return ParabolicCylinderD_Series::AtTransitionRecurUpForNegativeV(v, x);
  }
}

