
#include <functions/bessel/airy.h>
#include <functions/bessel/airy_util.h>
#include <functions/bessel/bessel.h>
#include <functions/gamma/gamma.h>
#include <functions/gamma/gamma_util.h>
#include <functions/hypergeometric/hypergeometric.h>


namespace AiryAiPrime_Series
{
  static e_float AtZero           (const e_float& x);
  static e_float AtTransitionPlus (const e_float& x);
  static e_float AtTransitionMinus(const e_float& xm);
  static e_float AtInfinityPlus   (const e_float& x);
  static e_float AtInfinityMinus  (const e_float& xm);
}

static e_float AiryAiPrime_Series::AtZero(const e_float& x)
{
  // Use the Taylor series expansion for |x| small.
  // http://functions.wolfram.com/BesselAiryStruveFunctions/AiryAiPrime/06/01/01/

  static const e_float f1         = AiryUtil::Ai_and_AiPrime_OfZero().front() / static_cast<INT32>(2);
  static const e_float f2         = AiryUtil::Ai_and_AiPrime_OfZero().back();
  static const e_float five_third = ef::five() / static_cast<INT32>(3);

  const e_float x_cubed_over_nine = ef::pown(x, static_cast<INT64>(3)) / ef::nine();

  return   ((f1 * (x * x)) * ef::hyperg_0f1(five_third,  x_cubed_over_nine))
         + ( f2            * ef::hyperg_0f1(ef::third(), x_cubed_over_nine));
}

static e_float AiryAiPrime_Series::AtTransitionPlus(const e_float& x)
{
  // Use the representation for x in the positive transition region.

  static const e_float one_over_sqrt3_pi = ef::one() / (ef::sqrt3() * ef::pi());

  return (-x * ef::cyl_bessel_k(ef::two_third(), AiryUtil::make_zeta(ef::sqrt(x)))) * one_over_sqrt3_pi;
}

static e_float AiryAiPrime_Series::AtTransitionMinus(const e_float& xm)
{
  // Use the representation for x in the negative transition region
  // with xm = |x|.

  e_float Jp1_3, Jm1_3, Jp2_3, Jm2_3;

  AiryUtil::BesselJ_Of_1_3and2_3(AiryUtil::make_zeta(ef::sqrt(xm)), Jp1_3, Jm1_3, Jp2_3, Jm2_3);

  return (ef::third() * xm) * (Jp2_3 - Jm2_3);
}

static e_float AiryAiPrime_Series::AtInfinityPlus(const e_float& x)
{
  // Use the hypergeometric representation in exponential form using
  // hyperg_2f0 for large positive values of x.
  // http://functions.wolfram.com/BesselAiryStruveFunctions/AiryAiPrime/06/02/01/01/

  static const e_float minus_sixth  = ef::one_minus() / static_cast<INT32>(6);
  static const e_float seven_sixths = ef::seven()     / static_cast<INT32>(6);
  
  const e_float sqrt_x = ef::sqrt(x);
  const e_float X      = AiryUtil::make_zeta(sqrt_x);
  const e_float sum    = ef::hyperg_2f0(minus_sixth, seven_sixths, ef::one_minus() / (X * static_cast<INT64>(2)));

  const e_float factor = -(ef::exp(-X) * ef::sqrt(sqrt_x)) / (ef::sqrt_pi() * static_cast<INT32>(2));

  return factor * sum;
}

static e_float AiryAiPrime_Series::AtInfinityMinus(const e_float& xm)
{
  // Use the hypergeometric representation in trigonometric form using
  // Hypergeometric4F1 for large negative values of x with xm = |x|.
  // http://functions.wolfram.com/BesselAiryStruveFunctions/AiryAiPrime/06/02/01/02/

  static const std::tr1::array<e_float, 4u> p1_data =
  {{
    ef::one_minus() / static_cast<INT32>(12),
    ef::five()      / static_cast<INT32>(12),
    ef::seven()     / static_cast<INT32>(12),
    e_float(13u)    / static_cast<INT32>(12)
  }};

  static const std::deque<e_float> p1(p1_data.begin(), p1_data.end());
  static const std::deque<e_float> q1(static_cast<std::size_t>(1u), ef::half());

  static const std::tr1::array<e_float, 4u> p2_data =
  {{
    p1_data[1u],
    e_float(11u) / static_cast<INT32>(12),
    p1_data[3u],
    e_float(19u) / static_cast<INT32>(12)
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

  const e_float sum = sum1 - ((sum2 * static_cast<e_float>(7)) / ((xm * sqrt_x) * static_cast<INT32>(48)));
  
  return -(sum * ef::sqrt(sqrt_x)) / ef::sqrt_pi();
}

e_float ef::airy_a_prime(const e_float& x)
{
  if(ef::iszero(x))
  {
    return AiryUtil::Ai_and_AiPrime_OfZero().back();
  }
  else
  {
    static const e_float nine_pow_one_third = ef::cbrt(ef::nine());

    const e_float xx = ef::fabs(x);

    if(xx <= nine_pow_one_third)
    {
      return AiryAiPrime_Series::AtZero(x);
    }
    else if(xx <= ef::hundred())
    {
      return ef::isneg(x) ? AiryAiPrime_Series::AtTransitionMinus(xx)
                          : AiryAiPrime_Series::AtTransitionPlus (x);
    }
    else
    {
      return ef::isneg(x) ? AiryAiPrime_Series::AtInfinityMinus(xx)
                          : AiryAiPrime_Series::AtInfinityPlus (x);
    }
  }
}
