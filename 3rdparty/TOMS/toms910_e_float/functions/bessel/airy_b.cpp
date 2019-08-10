
#include <functions/bessel/airy.h>
#include <functions/bessel/airy_util.h>
#include <functions/hypergeometric/hypergeometric.h>
#include <functions/bessel/bessel.h>


namespace AiryBi_Series
{
  static e_float AtZero           (const e_float& x);
  static e_float AtTransitionPlus (const e_float& x);
  static e_float AtTransitionMinus(const e_float& xm);
  static e_float AtInfinityPlus   (const e_float& x);
  static e_float AtInfinityMinus  (const e_float& xm);
}

static e_float AiryBi_Series::AtZero(const e_float& x)
{
  // Use the Taylor series expansion for |x| small.
  // http://functions.wolfram.com/Bessel-TypeFunctions/AiryBi/06/01/02/01/

  static const e_float sixth_root_of_three = ef::sqrt(AiryUtil::Three_Pow_1_3and2_3()[0]);

  static const e_float f1 = ef::one() / (sixth_root_of_three * AiryUtil::Gamma_Of_1_3and2_3().back());
  static const e_float f2 = sixth_root_of_three / AiryUtil::Gamma_Of_1_3and2_3().front();

  const e_float x_cubed_over_nine = ef::pown(x, static_cast<INT64>(3)) / ef::nine();

  return   ( f1      * ef::hyperg_0f1(ef::two_third(),  x_cubed_over_nine))
         + ((f2 * x) * ef::hyperg_0f1(ef::four_third(), x_cubed_over_nine));
}

static e_float AiryBi_Series::AtTransitionPlus(const e_float& x)
{
  // Use the representation for x in the positive transition region.

  const e_float sqrt_x = ef::sqrt(x);
  const e_float X      = AiryUtil::make_zeta(sqrt_x);

  return (sqrt_x * (ef::cyl_bessel_i(ef::third(), X) + ef::cyl_bessel_i(-ef::third(), X))) / ef::sqrt3();
}

static e_float AiryBi_Series::AtTransitionMinus(const e_float& xm)
{
  // Use the representation for x in the negative transition region
  // with xm = |x|.

  const e_float sqrt_x = ef::sqrt(xm);

  e_float Jp1_3, Jm1_3, Jp2_3, Jm2_3;

  AiryUtil::BesselJ_Of_1_3and2_3(AiryUtil::make_zeta(sqrt_x), Jp1_3, Jm1_3, Jp2_3, Jm2_3);

  return (sqrt_x / ef::sqrt3()) * (Jm1_3 - Jp1_3);
}

static e_float AiryBi_Series::AtInfinityPlus(const e_float& x)
{
  // Use the hypergeometric representation in exponential form using
  // hyperg_2f0 for large positive values of x.
  // http://functions.wolfram.com/Bessel-TypeFunctions/AiryBi/06/02/01/01/

  static const e_float sixth       = ef::one()  / static_cast<INT32>(6);
  static const e_float five_sixths = ef::five() / static_cast<INT32>(6);
  
  const e_float sqrt_x = ef::sqrt(x);
  const e_float X      = AiryUtil::make_zeta(sqrt_x);
  const e_float sum    = ef::hyperg_2f0(sixth, five_sixths, ef::one() / (X * static_cast<INT32>(2)));

  const e_float factor = ef::exp(X) / (ef::sqrt(sqrt_x) * ef::sqrt_pi());

  return factor * sum;
}

static e_float AiryBi_Series::AtInfinityMinus(const e_float& xm)
{
  // Use the hypergeometric representation in trigonometric form using
  // Hypergeometric4F1 for large negative values of x with xm = |x|.
  // http://functions.wolfram.com/Bessel-TypeFunctions/AiryBi/06/02/01/02/

  static const std::tr1::array<e_float, 4u> p1_data =
  {{
    ef::one()    / static_cast<INT32>(12),
    ef::five()   / static_cast<INT32>(12),
    ef::seven()  / static_cast<INT32>(12),
    e_float(11u) / static_cast<INT32>(12)
  }};

  static const std::deque<e_float> p1(p1_data.begin(), p1_data.end());
  static const std::deque<e_float> q1(static_cast<std::size_t>(1u), ef::half());

  static const std::tr1::array<e_float, 4u> p2_data =
  {{
    p1_data[2u],
    p1_data[3u],
    e_float(13u) / static_cast<INT32>(12),
    e_float(17u) / static_cast<INT32>(12)
  }};

  static const std::deque<e_float> p2(p2_data.begin(), p2_data.end());
  static const std::deque<e_float> q2(static_cast<std::size_t>(1u), ef::three_half());

  const e_float sqrt_x = ef::sqrt(xm);
  const e_float X      = AiryUtil::make_zeta(sqrt_x);
  
  e_float sin_X_plus_pi_quarter;
  e_float cos_X_plus_pi_quarter;
  
  ef::sincos(X + ef::pi_quarter(),
             &sin_X_plus_pi_quarter,
             &cos_X_plus_pi_quarter);

  const e_float minus_one_over_Xsq = ef::one_minus() / (X * X);

  const e_float sum1 = cos_X_plus_pi_quarter * ef::hyperg_pfq(p1, q1, minus_one_over_Xsq);
  const e_float sum2 = sin_X_plus_pi_quarter * ef::hyperg_pfq(p2, q2, minus_one_over_Xsq);

  const e_float sum = sum1 + ((sum2 * static_cast<INT32>(5)) / (X * static_cast<INT32>(72)));
  
  return sum / ef::sqrt(sqrt_x * ef::pi());
}


e_float ef::airy_b(const e_float& x)
{
  if(ef::iszero(x))
  {
    return AiryUtil::Bi_and_BiPrime_OfZero().front();
  }
  else
  {
    static const e_float nine_pow_one_third = ef::cbrt(ef::nine());

    const e_float xx = ef::fabs(x);

    if(xx <= nine_pow_one_third)
    {
      return AiryBi_Series::AtZero(x);
    }
    else if(xx <= ef::hundred())
    {
      return ef::isneg(x) ? AiryBi_Series::AtTransitionMinus(xx)
                          : AiryBi_Series::AtTransitionPlus (x);
    }
    else
    {
      return ef::isneg(x) ? AiryBi_Series::AtInfinityMinus(xx)
                          : AiryBi_Series::AtInfinityPlus (x);
    }
  }
}
