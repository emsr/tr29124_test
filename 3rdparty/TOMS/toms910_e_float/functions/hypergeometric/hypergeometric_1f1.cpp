
#include <functions/complex/e_float_complex.h>
#include <functions/constants/constants.h>
#include <functions/elementary/elementary.h>
#include <functions/bessel/bessel_recursion.h>
#include <functions/gamma/gamma.h>
#include <functions/hypergeometric/hypergeometric.h>
#include <functions/hypergeometric/hypergeometric_util.h>
#include <functions/hypergeometric/laguerre.h>
#include <functions/integer/integer.h>
#include <utility/util_digit_scale.h>
#include <utility/util_interpolate.h>

namespace Hypergeometric1F1_Series
{
  static e_float AtZero                                              (const e_float& a, const e_float& b, const e_float& x);
  static e_float AtInfinity_Regularized                              (const e_float& a, const e_float& b, const e_float& x);
  static e_float AtRecur_AB_Both_AB_PositiveInteger                  (const INT32 A,    const INT32 B,    const e_float& x);
  static e_float AtRecur_B_down_A_Real_B_NegativeReal                (const e_float& a, const e_float& b, const e_float& x);

  #if defined(NEVER_EXTENDED_H1F1)
  static e_float AtRecur_A_up_Both_AB_PositiveReal                   (const e_float& a, const e_float& b, const e_float& x);
  static e_float AtRecur_B_down_A_NegativeReal_B_PositiveRealExtended(const e_float& a, const e_float& b, const e_float& x);

  static e_float AtAbramowitzStegunEq13_3_6(const e_float& a, const e_float& b, const e_float& x);

  static e_float AtAbramowitzStegunEq13_3_7(const e_float& a,
                                            const e_float& b,
                                                  e_float& an,
                                                  e_float& anm1,
                                                  e_float& anm2,
                                            const INT32 n);
  #endif

  static e_float AtAbramowitzStegunEq13_3_8(const e_float& a,
                                            const e_float& b,
                                            const e_float& h,
                                                  e_float& cn,
                                                  e_float& cnm1,
                                                  e_float& cnm2,
                                            const INT32 n);

  static INT32 AtTransitionOrderAdder(const e_float& x);
}

static e_float Hypergeometric1F1_Series::AtZero(const e_float& a, const e_float& b, const e_float& x)
{
  return ef::hyp1F1(a, b, x);
}

static INT32 Hypergeometric1F1_Series::AtTransitionOrderAdder(const e_float& x)
{
  static const std::tr1::array<Util::point<double>, static_cast<std::size_t>(4u)> scale_x_data =
  {{
    Util::point<double>(static_cast<double>(  10.0), static_cast<double>(21.0)),
    Util::point<double>(static_cast<double>(  20.0), static_cast<double>(18.0)),
    Util::point<double>(static_cast<double>(  50.0), static_cast<double>(15.0)),
    Util::point<double>(static_cast<double>(1000.0), static_cast<double>( 3.0))
  }};

  static const std::vector<Util::point<double> > scale_x(scale_x_data.begin(), scale_x_data.end());

  const double xabs_d = ef::to_double(ef::fabs(x));

  static const double scale    = Util::DigitScale(); ;
  static const double my_scale = ((scale < static_cast<double>(1.0)) ? static_cast<double>(1.0) : scale);

  const double order_adder = (xabs_d * Util::linear_interpolate<double>::interpolate(xabs_d, scale_x)) * my_scale;

  return static_cast<INT32>(static_cast<INT64>(order_adder));
}

#if defined(NEVER_EXTENDED_H1F1)
static e_float Hypergeometric1F1_Series::AtAbramowitzStegunEq13_3_6(const e_float& a, const e_float& b, const e_float& x)
{
  // Compute the series of hyperg_1f1 using a series of Bessel functions In.
  // See Abramowitz and Stegun 13.3.6, page 506.

  const e_float b_minus_a            = b - a;
  const e_float b_minus_a_minus_half = b_minus_a - ef::half();
  const e_float x_over_two           = x / static_cast<INT32>(2);

  std::size_t n_ofs;

  std::deque<e_float> Iv;

  static_cast<void>(BesselRecursion::RecurIv(b_minus_a_minus_half,
                                             x_over_two,
                                             &Iv));

  n_ofs = static_cast<std::size_t>(ef::to_int32(b_minus_a_minus_half));

  // Do the BesselIv series expansion.

  e_float b_minus_a_minus_half_plus_n = b_minus_a_minus_half;

  e_float p1        = (b_minus_a * static_cast<INT32>(2)) - ef::one();
  e_float p2        = b_minus_a - a;
  e_float pb        = b;
  e_float pochham_1 = p1;
  e_float pochham_2 = p2;
  e_float pochham_b = b;
  e_float n_fact    = ef::one();

  e_float sum = Iv[n_ofs] * b_minus_a_minus_half_plus_n;

  ++b_minus_a_minus_half_plus_n;

  bool b_neg_term = true;

  for(INT32 n = static_cast<INT32>(1u); n < static_cast<INT32>(Iv.size()) - n_ofs; n++)
  {
    n_fact *= static_cast<INT32>(n);

    const e_float factor = ((pochham_1 * pochham_2) * b_minus_a_minus_half_plus_n) / (pochham_b * n_fact);
    const e_float term   = Iv[static_cast<std::size_t>(n + n_ofs)] * factor;

    const INT64 order_check = static_cast<INT64>(term.order() - sum.order());

    if((n > static_cast<std::size_t>(20u)) && (order_check < -ef::tol()))
    {
      break;
    }

    !b_neg_term ? sum += term : sum -= term;

    ++b_minus_a_minus_half_plus_n;

    pochham_1 *= ++p1;
    pochham_2 *= ++p2;
    pochham_b *= ++pb;

    b_neg_term = !b_neg_term;
  }

  const e_float scale = (ef::exp(x_over_two) * ef::gamma(b_minus_a_minus_half)) * ef::pow(x / static_cast<INT32>(4), ef::half() - b_minus_a);

  return scale * sum;
}

static e_float Hypergeometric1F1_Series::AtAbramowitzStegunEq13_3_7(const e_float& a,
                                                                    const e_float& b,
                                                                          e_float& an,
                                                                          e_float& anm1,
                                                                          e_float& anm2,
                                                                    const INT32 n)
{
  const e_float term_anm1 = (b + static_cast<INT32>(n - static_cast<INT32>(1))) * anm1;
  const e_float term_anm2 = ((a * static_cast<INT32>(2)) - b) * anm2;
  const e_float anp1      = (term_anm1 + term_anm2) / static_cast<INT32>(n + static_cast<INT32>(1));

  anm2 = anm1;
  anm1 = an;
  an   = anp1;

  return anp1;
}

namespace Hypergeometric1F1_Series
{
  template<typename T> static inline e_float AtTransitionAbramowitzStegunEq13_3_7_Regularized_Template(const e_float& a, const e_float& b, const e_float& x)
  {
    using ef::abs;
    using efz::abs;
    using ef::pow;
    using efz::pow;
    using ef::sqrt;
    using efz::sqrt;
    using ef::real;
    using efz::real;

    const e_float b_minus_one               = b - static_cast<INT32>(1);
    const e_float two_bx_minus_four_ax      = x * (static_cast<INT32>(2) * (b - (a * static_cast<INT32>(2))));
    const T       sqrt_two_bx_minus_four_ax = sqrt(T(two_bx_minus_four_ax));

    std::size_t n_ofs;

    std::deque<T> Jv;

    if(ef::isneg(b_minus_one))
    {
      const e_float v = ef::one() + ef::decimal_part(b_minus_one);

      static_cast<void>(BesselRecursion::RecurJv(v + AtTransitionOrderAdder(abs(x)),
                                                 sqrt_two_bx_minus_four_ax,
                                                 &Jv));

      n_ofs = static_cast<std::size_t>(0u);

      // Extend the list of Jv down into the negative orders using downward recursion.
      T jvp2 = Jv[static_cast<std::size_t>(1u)];
      T jvp1 = Jv[static_cast<std::size_t>(0u)];

      // Compute the number of orders which need to be recurred. Use a rounding
      // correction of 1/2 here.
      const INT32 N = static_cast<INT32>(static_cast<INT32>(1) + static_cast<INT32>(ef::to_int64(ef::half() + ef::integer_part(-b_minus_one))));

      e_float n_plus_one_plus_v = v;

      const T two_over_arg = ef::two() / sqrt_two_bx_minus_four_ax;

      for(INT32 m = static_cast<INT32>(0); m < N; m++)
      {
        Jv.push_front(((jvp1 * two_over_arg) * n_plus_one_plus_v) - jvp2);
        jvp2 = jvp1;
        jvp1 = Jv.front();

        --n_plus_one_plus_v;
      }
    }
    else
    {
      static_cast<void>(BesselRecursion::RecurJv(b_minus_one + AtTransitionOrderAdder(x),
                                                 sqrt_two_bx_minus_four_ax,
                                                 &Jv));

      n_ofs = static_cast<std::size_t>(static_cast<INT32>(ef::to_int64(b_minus_one)));
    }


    e_float x_pow_n                         = x * x;
    T       sqrt_two_bx_minus_four_ax_pow_n = two_bx_minus_four_ax;

    e_float anm2 = ef::one();
    e_float anm1 = ef::zero();
    e_float an   = b / static_cast<INT32>(2);

    T factor = an * (x_pow_n / sqrt_two_bx_minus_four_ax_pow_n);
    T term   = Jv[static_cast<std::size_t>(n_ofs + static_cast<std::size_t>(2u))] * factor;

    T sum = Jv[n_ofs] + term;

    for(std::size_t n = static_cast<std::size_t>(3u); n < static_cast<std::size_t>(Jv.size() - n_ofs); n++)
    {
      an = AtAbramowitzStegunEq13_3_7(a, b, an, anm1, anm2, static_cast<INT32>(n - static_cast<std::size_t>(1u)));

      x_pow_n                         *= x;
      sqrt_two_bx_minus_four_ax_pow_n *= sqrt_two_bx_minus_four_ax;
    
      factor = an * (x_pow_n / sqrt_two_bx_minus_four_ax_pow_n);
      term   = Jv[static_cast<std::size_t>(n + n_ofs)] * factor;

      const INT64 order_check = static_cast<INT64>(abs(term).order() - abs(sum).order());

      if((n > static_cast<std::size_t>(20u)) && (order_check < -ef::tol()))
      {
        break;
      }

      sum += term;
    }

    const T scale = ef::exp(x / static_cast<INT32>(2)) * pow(T(two_bx_minus_four_ax) / static_cast<INT32>(4), T(-b_minus_one) / static_cast<INT32>(2));

    return real(scale * sum);
  }
}
#endif

static e_float Hypergeometric1F1_Series::AtAbramowitzStegunEq13_3_8(const e_float& a,
                                                                    const e_float& b,
                                                                    const e_float& h,
                                                                          e_float& cn,
                                                                          e_float& cnm1,
                                                                          e_float& cnm2,
                                                                    const INT32 n)
{
  const e_float one_minus_two_h = ef::one() - (h * static_cast<std::size_t>(2u));
  const e_float h_minus_one     = h - ef::one();

  const e_float term_cn   = ((one_minus_two_h * n) - (b * h)) * cn;
  const e_float term_cnm1 = ((one_minus_two_h * a) - ((h * h_minus_one) * (b + static_cast<INT32>(n - static_cast<INT32>(1))))) * cnm1;
  const e_float term_cnm2 = ((-h * h_minus_one) * a) * cnm2;

  const e_float cnp1 = ((term_cn + term_cnm1) + term_cnm2) / static_cast<INT32>(n + static_cast<INT32>(1));

  cnm2 = cnm1;
  cnm1 = cn;
  cn   = cnp1;

  return cnp1;
}

namespace Hypergeometric1F1_Series
{
  template<typename T> static inline e_float AtTransitionAbramowitzStegunEq13_3_8_Regularized_Template(const e_float& a, const e_float& b, const e_float& x)
  {
    using ef::abs;
    using efz::abs;
    using ef::pow;
    using efz::pow;
    using ef::sqrt;
    using efz::sqrt;
    using ef::real;
    using efz::real;

    // Compute the series of hyperg_1f1 using a series of Bessel functions In.
    // See Abramowitz and Stegun 13.3.8, page 506.

    const e_float minus_ax      = -a * x;
    const T       sqrt_minus_ax = sqrt(T(minus_ax));
    const e_float b_minus_one   = b - ef::one();

    std::size_t n_ofs;

    std::deque<T> Jv;

    if(ef::isneg(b_minus_one))
    {
      const e_float v = ef::one() + ef::decimal_part(b_minus_one);

      static_cast<void>(BesselRecursion::RecurJv(v + AtTransitionOrderAdder(x),
                                                 sqrt_minus_ax * static_cast<INT32>(2),
                                                 &Jv));

      n_ofs = static_cast<std::size_t>(0u);

      // Extend the list of Jv down into the negative orders using downward recursion.
      T jvp2 = Jv[static_cast<std::size_t>(1u)];
      T jvp1 = Jv[static_cast<std::size_t>(0u)];

      // Compute the number of orders which need to be recurred. Use a rounding
      // correction of 1/2 here.
      const INT32 N = static_cast<INT32>(static_cast<INT32>(1) + static_cast<INT32>(ef::to_int64(ef::half() + ef::integer_part(-b_minus_one))));

      e_float n_plus_one_plus_v = v;

      const T two_over_arg = ef::one() / sqrt_minus_ax;

      for(INT32 m = static_cast<INT32>(0); m < N; m++)
      {
        Jv.push_front(((jvp1 * two_over_arg) * n_plus_one_plus_v) - jvp2);
        jvp2 = jvp1;
        jvp1 = Jv.front();

        --n_plus_one_plus_v;
      }
    }
    else
    {
      static_cast<void>(BesselRecursion::RecurJv(b_minus_one + AtTransitionOrderAdder(x),
                                                 sqrt_minus_ax * static_cast<INT32>(2),
                                                 &Jv));

      n_ofs = static_cast<std::size_t>(static_cast<INT32>(ef::to_int64(b_minus_one)));
    }

    static const e_float h = -ef::pi() / ef::ten();

    // Do the BesselJv series expansion.

    e_float cnm2 = ef::one();
    e_float cnm1 = -b * h;
    e_float cn   = (((ef::one() - (h * static_cast<INT32>(2))) * a) + ((b * (b + ef::one())) * (h * h))) / static_cast<INT32>(2);

    e_float x_pow_n             = x * x;
    T       sqrt_minus_ax_pow_n = minus_ax;

    T sum =   Jv[n_ofs]
            + Jv[static_cast<std::size_t>(n_ofs + static_cast<std::size_t>(1u))] * ((cnm1 * x)       / sqrt_minus_ax)
            + Jv[static_cast<std::size_t>(n_ofs + static_cast<std::size_t>(2u))] * ((cn   * x_pow_n) / sqrt_minus_ax_pow_n);

    for(std::size_t n = static_cast<std::size_t>(3u); n < static_cast<std::size_t>(Jv.size() - n_ofs); n++)
    {
      cn = AtAbramowitzStegunEq13_3_8(a, b, h, cn, cnm1, cnm2, static_cast<INT32>(n - static_cast<std::size_t>(1u)));

      x_pow_n             *= x;
      sqrt_minus_ax_pow_n *= sqrt_minus_ax;

      const T factor = (cn * x_pow_n) / sqrt_minus_ax_pow_n;
      const T term   = Jv[static_cast<std::size_t>(n + n_ofs)] * factor;

      const INT64 order_check = static_cast<INT64>(abs(term).order() - abs(sum).order());

      if((n > static_cast<std::size_t>(20u)) && (order_check < -ef::tol()))
      {
        break;
      }

      sum += term;
    }

    return real((sum * ef::exp(h * x)) / pow(sqrt_minus_ax, T(b_minus_one)));
  }
}

static e_float Hypergeometric1F1_Series::AtInfinity_Regularized(const e_float& a, const e_float& b, const e_float& x)
{
  const e_float b_minus_a = b - a;

  const e_float factor = (ef::pow(x, -b_minus_a) * ef::exp(x)) / ef::gamma(a);
  const e_float term   = ef::hyperg_2f0(b_minus_a, ef::one() - a, ef::one() / x);

  return factor * term;
}

static e_float Hypergeometric1F1_Series::AtRecur_AB_Both_AB_PositiveInteger(const INT32 A, const INT32 B, const e_float& x)
{
  // Recur across in the direction of increasing b and then up in the direction
  // of increasing a.

  const e_float exp_x = ef::exp(x);

  e_float Mb_m2 =  exp_x;                  // Abramowitz and Stegun 13.6.12 (Special Cases), page 509
  e_float Mb_m1 = (exp_x - ef::one()) / x; // Abramowitz and Stegun 13.6.14 (Special Cases), page 509
  e_float Mb;

  if(B == static_cast<INT32>(1))
  {
    Mb = Mb_m2;
  }
  else if(B == static_cast<INT32>(2))
  {
    Mb = Mb_m1;
  }
  else
  {
    for(INT32 b = static_cast<INT32>(3); b <= B; b++)
    {
      // Abramowitz and Stegun Equation 13.4.2.
      const INT32 b_minus_one = static_cast<INT32>(b - static_cast<INT32>(1));
      const INT32 b_minus_two = static_cast<INT32>(b - static_cast<INT32>(2));

      Mb = (((b_minus_one * (b_minus_two + x)) * Mb_m1) - ((b_minus_one * b_minus_two) * Mb_m2)) / (x * b_minus_two);

      Mb_m2 = Mb_m1;
      Mb_m1 = Mb;
    }
  }

  e_float Ma_m2 = ef::one();
  e_float Ma_m1 = Mb;
  e_float Ma;

  if(A == static_cast<INT32>(1))
  {
    Ma = Ma_m1;
  }
  else
  {
    for(INT32 a = static_cast<INT32>(2); a <= A; a++)
    {
      // Abramowitz and Stegun Equation 13.4.1.
      const INT32     a_minus_one = static_cast<INT32>(a - static_cast<INT32>(1));
      const INT32 two_a_minus_one = static_cast<INT32>(a_minus_one * static_cast<INT32>(2));

      Ma = ((((two_a_minus_one - B) + x) * Ma_m1) + ((B - a_minus_one) * Ma_m2)) / a_minus_one;

      Ma_m2 = Ma_m1;
      Ma_m1 = Ma;
    }
  }

  return Ma;
}

static e_float Hypergeometric1F1_Series::AtRecur_B_down_A_Real_B_NegativeReal(const e_float& a, const e_float& b, const e_float& x)
{
  // Recur down in the direction of decreasing b, in other words increasingly negative b.

  const e_float b_decimal = ef::decimal_part(b);

  const INT32 B = ef::to_int32(b);

  e_float Mb_p2 = ef::conf_hyperg(a, b_decimal, x);
  e_float Mb_p1 = ef::conf_hyperg(a, b_decimal - ef::one(), x);
  e_float Mb;

  for(INT32 bn = static_cast<INT32>(-2); bn >= B; bn--)
  {
    // Abramowitz and Stegun Equation 13.4.2.
    const e_float b_plus_zero = b_decimal + bn;
    const e_float b_plus_one  = b_decimal + (bn + static_cast<INT32>(1));

    Mb = (((b_plus_one * (b_plus_zero + x)) * Mb_p1) - ((x * (b_plus_one - a)) * Mb_p2)) / (b_plus_zero * b_plus_one);

    Mb_p2 = Mb_p1;
    Mb_p1 = Mb;
  }

  return Mb;
}

#if defined(NEVER_EXTENDED_H1F1)
static e_float Hypergeometric1F1_Series::AtRecur_A_up_Both_AB_PositiveReal(const e_float& a, const e_float& b, const e_float& x)
{
  // Recur up in the direction of increasing a.

  const e_float a_decimal = ef::decimal_part(a);

  const INT32 A = ef::to_int32(a);

  e_float Ma_m2 = ef::conf_hyperg(a_decimal, b, x);
  e_float Ma_m1 = ef::conf_hyperg(a_decimal + ef::one(), b, x);
  e_float Ma;

  if(A == static_cast<INT32>(1))
  {
    Ma = Ma_m1;
  }
  else
  {
    for(INT32 an = static_cast<INT32>(2); an <= A; an++)
    {
      // Abramowitz and Stegun Equation 13.4.1.
      const e_float     a_minus_one = a_decimal + (an - static_cast<INT32>(1));
      const e_float two_a_minus_one = a_minus_one * static_cast<INT32>(2);

      Ma = ((((two_a_minus_one - b) + x) * Ma_m1) + ((b - a_minus_one) * Ma_m2)) / a_minus_one;

      Ma_m2 = Ma_m1;
      Ma_m1 = Ma;
    }
  }

  return Ma;
}

static e_float Hypergeometric1F1_Series::AtRecur_B_down_A_NegativeReal_B_PositiveRealExtended(const e_float& a, const e_float& b, const e_float& x)
{
  // Extend b to a positive real number which is greater than 2 * (|a| + x).
  const e_float a_abs      = ef::fabs(a);
  const e_float b_extended = ((a_abs + x) * static_cast<INT32>(2)) + static_cast<INT32>(10);

  if(b > b_extended)
  {
    return Hypergeometric1F1_Series::AtTransitionAbramowitzStegunEq13_3_8_Regularized_Template<e_float>(a, b, x) * ef::gamma(b);
  }

  // Recur down from b_extended to the desired value of b.
  const e_float b_decimal  = ef::decimal_part(b);
  const INT32   B          = ef::to_int32(b_extended);
  const e_float b_plus_two = (B + b_decimal) + static_cast<INT32>(2);
        e_float b_plus_one = (B + b_decimal) + static_cast<INT32>(1);

  e_float Mb_p2 = Hypergeometric1F1_Series::AtTransitionAbramowitzStegunEq13_3_8_Regularized_Template<e_float>(a, b_plus_two, x) * ef::gamma(b_plus_two);
  e_float Mb_p1 = Hypergeometric1F1_Series::AtTransitionAbramowitzStegunEq13_3_8_Regularized_Template<e_float>(a, b_plus_one, x) * ef::gamma(b_plus_one);
  e_float Mb;

  const INT32 Bb = ef::to_int32(b);

  for(INT32 bn = B; bn >= Bb; bn--)
  {
    // Abramowitz and Stegun Equation 13.4.2.
    const e_float b_plus_zero = b_decimal + bn;
                  b_plus_one  = b_decimal + (bn + static_cast<INT32>(1));

    Mb = (((b_plus_one * (b_plus_zero + x)) * Mb_p1) - ((x * (b_plus_one - a)) * Mb_p2)) / (b_plus_zero * b_plus_one);

    Mb_p2 = Mb_p1;
    Mb_p1 = Mb;
  }

  return Mb;
}
#endif

e_float ef::hyperg_1f1(const e_float& a, const e_float& b, const e_float& x)
{
  // Check for various special values.

  if(ef::iszero(a) || ef::iszero(x))
  {
    return ef::one();
  }

  if(ef::iszero(b))
  {
    return std::numeric_limits<e_float>::quiet_NaN();
  }

  const e_float b_minus_a = b - a;

  if(ef::iszero(b_minus_a))
  {
    return ef::exp(x);
  }

  if(ef::isone(-a))
  {
    return ef::one() - (x / b);
  }

  if(ef::isone(-b_minus_a))
  {
    return (ef::one() + (x / b)) * ef::exp(x);
  }

  // Check for singularities and complex infinities in case a or b is a negative
  // integer or if both a and b are negative integers.

  const bool a_is_integer = ef::isint(a);
  const bool b_is_integer = ef::isint(b);

  const bool a_is_neg = ef::isneg(a);
  const bool b_is_neg = ef::isneg(b);

  const bool a_is_positive_integer = (a_is_integer && !a_is_neg);
  const bool a_is_negative_integer = (a_is_integer &&  a_is_neg);
  const bool b_is_negative_integer = (b_is_integer &&  b_is_neg);

  if(b_is_negative_integer && !a_is_negative_integer)
  {
    return std::numeric_limits<e_float>::quiet_NaN();
  }

  if(b_is_negative_integer && a_is_negative_integer)
  {
    const INT32 m = ef::to_int32(a);
    const INT32 n = ef::to_int32(b);

    if(m < n)
    {
      return std::numeric_limits<e_float>::quiet_NaN();
    }
    else
    {
      // Abramowitz and Stegun 13.6.9 (Special Cases), page 509.
      const UINT32  N              = static_cast<INT32>(-m);
      const INT32   alpha          = static_cast<INT32>(n - static_cast<INT32>(1));
      const e_float alpha_plus_one = e_float(n);
      const e_float Lna            = ef::laguerre(static_cast<INT32>(N), alpha, x);

      return (ef::factorial(N) * Lna) / ef::pochhammer(alpha_plus_one, N);
    }
  }

  if(a_is_negative_integer)
  {
    const INT32 m = -ef::to_int32(a);

    return ((ef::factorial(static_cast<UINT32>(m)) * ef::gamma(b)) * ef::laguerre(m, b - ef::one(), x)) / ef::gamma(b + m);
  }

  if(ef::isneg(x))
  {
    return ef::exp(x) * ef::conf_hyperg(b_minus_a, b, -x);
  }

  if(a_is_positive_integer && (b_is_integer && !b_is_neg))
  {
    // Both a and b are positive integers. Use forward recursion schemes
    // in a and b to recur up and across to the desired value.
    const INT32 Am = ef::to_int32(a);
    const INT32 Bn = ef::to_int32(b);

    return Hypergeometric1F1_Series::AtRecur_AB_Both_AB_PositiveInteger(Am, Bn, x);
  }

  const double xd = ef::to_double(x);
  const double ad = ef::to_double(a);
  const double bd = ef::to_double(b);

  const bool bo_converge_small = HypergeometricUtil::AsympConverge(std::deque<double>(static_cast<std::size_t>(1u), ad),
                                                                   std::deque<double>(static_cast<std::size_t>(1u), bd),
                                                                   xd);

  if(bo_converge_small)
  {
    return Hypergeometric1F1_Series::AtZero(a, b, x);
  }

  static const double one_d                = static_cast<double>(1.0);
         const double bd_minus_ad          = static_cast<double>(bd - ad);
         const double one_minus_ad         = static_cast<double>(one_d - ad);

  const std::tr1::array<double, 2u> a_data = {{ bd_minus_ad, one_minus_ad }};

  const bool bo_converge_large = HypergeometricUtil::AsympConverge(std::deque<double>(a_data.begin(), a_data.end()),
                                                                   std::deque<double>(),
                                                                   static_cast<double>(one_d / xd));                                                                     

  if(bo_converge_large && !a_is_positive_integer)
  {
    return Hypergeometric1F1_Series::AtInfinity_Regularized(a, b, x) * ef::gamma(b);
  }

  if(b < -ef::ten())
  {
    return Hypergeometric1F1_Series::AtRecur_B_down_A_Real_B_NegativeReal(a, b, x);
  }

  if(ef::isneg(a) != ef::isneg(x))
  {
    // The non-transformed arguments are appropriate for real-valued expansion
    // using Abramowitz and Stegun Equation 13.3.8.
    return Hypergeometric1F1_Series::AtTransitionAbramowitzStegunEq13_3_8_Regularized_Template<e_float>(a, b, x) * ef::gamma(b);
  }
  else
  {
    // The non-transformed arguments are appropriate for complex-valued expansion
    // using Abramowitz and Stegun Equation 13.3.8.
    return Hypergeometric1F1_Series::AtTransitionAbramowitzStegunEq13_3_8_Regularized_Template<ef_complex>(a, b, x) * ef::gamma(b);
  }
}

e_float ef::hyperg_1f1_reg(const e_float& a, const e_float& b, const e_float& x)
{
  const INT32 Bn = ef::to_int32(b);

  if(ef::isint(b) && (Bn < static_cast<INT32>(1)))
  {
    const INT32 n          = static_cast<INT32>(-Bn);
    const INT32 n_plus_one = static_cast<INT32>(n + static_cast<INT32>(1));

    const e_float factor = (ef::pochhammer(a, static_cast<UINT32>(n_plus_one)) * ef::pown(x, static_cast<INT64>(n_plus_one))) / ef::factorial(static_cast<UINT32>(n_plus_one));

    return factor * ef::hyperg_1f1(a + n_plus_one, n + ef::two(), x);
  }
  else
  {
    return ef::hyperg_1f1(a, b, x) / ef::gamma(b);
  }
}
