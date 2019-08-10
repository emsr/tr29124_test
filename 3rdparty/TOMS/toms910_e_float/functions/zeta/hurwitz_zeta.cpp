
#include <e_float/e_float.h>
#include <functions/constants/constants.h>
#include <functions/elementary/elementary.h>
#include <functions/gamma/gamma.h>
#include <functions/integer/integer.h>
#include <functions/zeta/zeta.h>
#include <utility/util_digit_scale.h>
#include <utility/util_interpolate.h>
#include <utility/util_power_j_pow_x.h>

namespace HurwitzZeta_Series
{
  template<typename T> inline static bool AtInfinityConvergenceCheckTemplate(const T& s, const T& a)
  {
    using ef::abs;
    using efz::abs;
    using ef::pow;
    using efz::pow;

    // The series representation which defines the Hurwitz Zeta function can be
    // used if the value of the denominator in the series is relatively large.
    // Select a maximum number of allowed series terms, say 167...500 based on
    // the number of digits.

    static const double max_scaled = static_cast<double>(Util::DigitScale() * static_cast<double>(500.0));
    static const double max_lower  = static_cast<double>(167.0);

    static const double number_of_series_terms_d = std::max(max_scaled, max_lower);
    static const INT32  number_of_series_terms   = static_cast<INT32>(static_cast<INT64>(number_of_series_terms_d));

    // Check if the order of the asymptotic series calculation will be adequate for
    // the precision of the final result.

    const T a_test = ef::one() / pow(a, s);
    const T s_test = ef::one() / pow(a + number_of_series_terms, s);

    const INT64 order_check = static_cast<INT64>(abs(s_test).order() - abs(a_test).order());

    return order_check < -ef::tol();
  }

  static bool AtExponentialFourierSeriesConvergenceCheck(const e_float& s, const e_float& a)
  {
    // This series is valid for large negative s and small to moderate a.
    // Although the actual series is only valid for 0 < a < 1, recursion
    // will be used for larger values of a.

    // Check for large negative s and a < 200. In order to qualify as large s,
    // the value of s should satisfy s < -20 ... -60, based on the digit range
    // of the calculation. 

    static const double  conv_check_s_scaled = static_cast<double>(Util::DigitScale() * static_cast<double>(-60.0));
    static const double  conv_check_s_upper  = static_cast<double>(-20.0);
    static const double  conv_check_s        = std::min(conv_check_s_scaled, conv_check_s_upper);

    static const e_float convergence_check_s(conv_check_s);

    return ((s < convergence_check_s) && (a < ef::two_hundred()));
  }

  template<typename T> inline static T AtInfinityTemplate(const T& s, const T& a)
  {
    // Use the series representation which defines the Hurwitz Zeta function.
    // http://functions.wolfram.com/ZetaFunctionsandPolylogarithms/Zeta2/02/

    using ef::abs;
    using efz::abs;
    using ef::pow;
    using efz::pow;
    using ef::pown;
    using efz::pown;
    using ef::real;
    using efz::real;

    const bool s_is_integer = ef::isint(s);

    const INT64 sn = s_is_integer ? ef::to_int64(real(s)) : static_cast<INT64>(0);

    T k_plus_a = a;
    T sum      = ef::zero();

    for(INT32 k = static_cast<INT32>(0); k < ef::max_iteration(); k++)
    {
      const T k_plus_a_pow_s = !s_is_integer ? pow (k_plus_a, s)
                                             : pown(k_plus_a, sn);

      const T term = ef::one() / k_plus_a_pow_s;

      const INT64 order_check = static_cast<INT64>(abs(term).order() - abs(sum).order());

      if((k > static_cast<INT32>(20)) && (order_check < -ef::tol()))
      {
        break;
      }

      sum += term;

      ++k_plus_a;
    }

    return sum;
  }

  template<typename T> inline static T AtEulerMaclaurinTemplate(const T& s, const T& a)
  {
    using ef::abs;
    using efz::abs;
    using ef::pow;
    using efz::pow;
    using ef::pown;
    using efz::pown;
    using ef::real;
    using efz::real;
    using ef::imag;
    using efz::imag;

    if(abs(imag(s)) > (ef::million() + ef::hundred()))
    {
      // Return NaN if s has a large imaginary part.
      return std::numeric_limits<e_float>::quiet_NaN();
    }

    // Use the Euler-Maclaurin summation formula with an appropriate choice of N.
    // http://functions.wolfram.com/ZetaFunctionsandPolylogarithms/Zeta2/06/05/01/01/

    // Use N = {(digits * 0.5) + 10} + {|imag(s)| * 1.1}
    static const double dfrac = static_cast<double>(0.5) * static_cast<double>(static_cast<INT32>(ef::tol()));
    static const double nd    = static_cast<double>(dfrac + static_cast<double>(10.0));
           const double ni    = static_cast<double>(static_cast<double>(1.10) * ef::to_double(imag(s)));
    static const INT32  N     = static_cast<INT32>(static_cast<INT64>(static_cast<double>(nd + ni)));

    const bool s_is_integer = (ef::isint(s) && (real(s) > ef::int64min()) && (real(s) < ef::int64max()));

    const INT64 sn = s_is_integer ? ef::to_int64(real(s)) : static_cast<INT64>(0);

    T k_plus_a = a;
    T sum0     = ef::zero();

    for(INT32 k = static_cast<INT32>(0); k < N; k++)
    {
      const T k_plus_a_pow_s = (!s_is_integer ? pow (k_plus_a, s)
                                              : pown(k_plus_a, sn));

      sum0 += ef::one() / k_plus_a_pow_s;

      ++k_plus_a;
    }

    const T N_plus_a         = k_plus_a;
    const T N_plus_a_pow_s   = pow(N_plus_a, s);
    const T term_0_1         = ((N_plus_a / (s - ef::one())) + ef::half()) / N_plus_a_pow_s;
    const T N_plus_a_squared = N_plus_a * N_plus_a;

    T N_plus_a_pow_s_plus_two_k_plus_one = N_plus_a_pow_s * N_plus_a;

    T sk         = s;
    T sk_term    = sk;
    T two_k_fact = ef::two();
    T sum1       = ef::zero();

    const e_float abs_sum0_plus_term_0_1 = abs(sum0 + term_0_1);

    for(INT32 k = static_cast<INT32>(1); k < ef::max_iteration(); k++)
    {
      const INT32 two_k = static_cast<INT32>(k * static_cast<INT32>(2));

      const T term =   (sk_term * ef::bernoulli(static_cast<UINT32>(two_k)))
                     / (two_k_fact * N_plus_a_pow_s_plus_two_k_plus_one);

      const INT64 order_check = static_cast<INT64>(abs(term).order() - (abs_sum0_plus_term_0_1 + abs(sum1)).order());

      if((k > static_cast<INT32>(20)) && (order_check < -ef::tol()))
      {
        break;
      }

      sum1 += term;

      two_k_fact *= static_cast<INT32>(static_cast<INT32>(two_k + static_cast<INT32>(1)) * static_cast<INT32>(two_k + static_cast<INT32>(2)));

      sk_term *= ++sk;
      sk_term *= ++sk;

      N_plus_a_pow_s_plus_two_k_plus_one *= N_plus_a_squared;
    }

    return (sum0 + term_0_1) + sum1;
  }

  template<typename T> inline static T AtExponentialFourierSeriesTemplate(const T& s, const T& a)
  {
    using ef::abs;
    using efz::abs;
    using ef::gamma;
    using efz::gamma;
    using ef::pow;
    using efz::pow;
    using ef::real;
    using efz::real;
    using ef::sincos;
    using efz::sincos;

    if(real(a) > ef::one())
    {
      // Use recursion for a slightly larger than 1.
      const e_float a_int_part = ef::integer_part(real(a));

      const T a_frac = a - a_int_part;

      const T result = AtExponentialFourierSeriesTemplate(s, a_frac);

      const INT32 J = static_cast<INT32>(ef::to_int64(a_int_part));

      T sum = ef::zero();

      const T minus_s = -s;

      for(INT32 j = static_cast<INT32>(0); j < J; j++)
      {
        sum += pow(j + a_frac, minus_s);
      }

      return result - sum;
    }

    std::map<UINT32, T> n_pow_x_prime_factor_map;

    const T one_minus_s = ef::one() - s;
    const T two_pi_a    = ef::two_pi() * a;

    T sin_two_pi_ak;
    T cos_two_pi_ak;

    // Calculate the k = 1 summation terms manually, outside of the loop.
    sincos(two_pi_a, &sin_two_pi_ak, &cos_two_pi_ak);

    const T C_x = cos_two_pi_ak;

    // Define the (n - 2) recursion variables for sin / cos recursions.
    T S_m2 = sin_two_pi_ak;
    T C_m2 = cos_two_pi_ak;

    T sum_sin = sin_two_pi_ak;
    T sum_cos = cos_two_pi_ak;

    // Calculate the k = 2 summation terms manually, outside of the loop.
    sincos(two_pi_a * static_cast<INT32>(2), &sin_two_pi_ak, &cos_two_pi_ak);

    // Define the (n - 1) recursion variables for sin / cos recursions.
    T S_m1 = sin_two_pi_ak;
    T C_m1 = cos_two_pi_ak;

    T          k_pow_one_minus_s = Util::j_pow_x(static_cast<UINT32>(2u), one_minus_s, n_pow_x_prime_factor_map);
    T one_over_k_pow_one_minus_s = ef::one() / k_pow_one_minus_s;

    sum_sin += (sin_two_pi_ak * one_over_k_pow_one_minus_s);
    sum_cos += (cos_two_pi_ak * one_over_k_pow_one_minus_s);

    for(INT32 k = static_cast<INT32>(3); k < ef::max_iteration(); k++)
    {
               k_pow_one_minus_s = Util::j_pow_x(static_cast<UINT32>(k), one_minus_s, n_pow_x_prime_factor_map);
      one_over_k_pow_one_minus_s = ef::one() / k_pow_one_minus_s;

      cos_two_pi_ak = ((C_m1 * C_x) * static_cast<INT32>(2)) - C_m2;
      sin_two_pi_ak = ((S_m1 * C_x) * static_cast<INT32>(2)) - S_m2;

      C_m2 = C_m1;
      C_m1 = cos_two_pi_ak;

      S_m2 = S_m1;
      S_m1 = sin_two_pi_ak;

      const T term_sin = sin_two_pi_ak * one_over_k_pow_one_minus_s;
      const T term_cos = cos_two_pi_ak * one_over_k_pow_one_minus_s;

      const INT64 order_check_sin = static_cast<INT64>(abs(term_sin).order() - abs(sum_sin).order());
      const INT64 order_check_cos = static_cast<INT64>(abs(term_cos).order() - abs(sum_cos).order());

      const bool term_sin_is_negligible = order_check_sin < -ef::tol();
      const bool term_cos_is_negligible = order_check_cos < -ef::tol();

      if(   (k > static_cast<INT32>(20))
         && term_sin_is_negligible
         && term_cos_is_negligible
        )
      {
        break;
      }

      sum_sin += term_sin;
      sum_cos += term_cos;
    }

    T sin_pi_half_s;
    T cos_pi_half_s;

    sincos(ef::pi_half() * s, &sin_pi_half_s, &cos_pi_half_s);

    T factor = (ef::two() * pow(T(ef::two_pi()), -one_minus_s)) * gamma(one_minus_s);

    return factor * ((cos_pi_half_s * sum_sin) + (sin_pi_half_s * sum_cos));
  }

  template<typename T> inline static T ZetaTemplate(const T& s, const T& a)
  {
    using ef::iszero;
    using ef::imag;
    using efz::imag;
    using ef::real;
    using efz::real;
    using ef::abs;
    using efz::abs;
    using ef::pow;
    using efz::pow;
    using ef::riemann_zeta;
    using efz::riemann_zeta;

    if(iszero(a))
    {
      return riemann_zeta(s);
    }

    const e_float sr = real(s);
    const e_float ar = real(a);

    if(ef::isneg(ar))
    {
      // Use the reflection of a given in one of the transformations at Wolfram's Function Site.
      // http://functions.wolfram.com/ZetaFunctionsandPolylogarithms/Zeta2/16/01/01/01/
      const INT32 n         = static_cast<INT32>(static_cast<INT64>(1) - ef::to_int64(real(a)));
      const T     aa        = T(ef::one()) + ef::decimal_part(real(a)) + (!iszero(imag(a)) ? (a - real(a)) : ef::zero());
      const T     s_half    = s / static_cast<INT32>(2);
      const T     a_minus_n = aa - n;

      T sum = ef::zero();

      for(INT32 k = static_cast<INT32>(0); k < n; k++)
      {
        const T a_plus_k_minus_n         = a_minus_n + k;
        const T a_plus_k_minus_n_squared = a_plus_k_minus_n * a_plus_k_minus_n;

        sum += ef::one() / pow(a_plus_k_minus_n_squared, s_half);
      }

      return ZetaTemplate(s, aa) + sum;
    }

    if(ef::iszero(imag(a)) && AtExponentialFourierSeriesConvergenceCheck(sr, ar))
    {
      // Support for large negative s.
      return AtExponentialFourierSeriesTemplate(s, a);
    }

    if(AtInfinityConvergenceCheckTemplate(s, a))
    {
      // Support for large positive s.
      return AtInfinityTemplate(s, a);
    }
    else
    {
      static const e_float minus_ten = -ef::ten();

      const bool s_is_moderately_negative = sr < minus_ten;

      const INT32 J = (s_is_moderately_negative && (ar < abs(sr))) ? static_cast<INT32>(ef::to_int64(abs(sr) - ar))
                                                                   : static_cast<INT32>(0);

      if(J == static_cast<INT32>(0))
      {
        return AtEulerMaclaurinTemplate(s, a);
      }
      else
      {
        const T result = AtEulerMaclaurinTemplate(s, a + J);

        const T minus_s = -s;

        T sum = ef::zero();

        for(INT32 j = static_cast<INT32>(0); j < J; j++)
        {
          sum += pow(j + a, minus_s);
        }

        return result + sum;
      }
    }
  }
}

e_float ef::hurwitz_zeta(const e_float& s, const INT32 n)
{
  return (n != static_cast<INT32>(0)) ? HurwitzZeta_Series::ZetaTemplate<e_float>(s, e_float(n))
                                      : ef::riemann_zeta(s);
}

ef_complex efz::hurwitz_zeta(const ef_complex& s, const INT32 n)
{
  return (n != static_cast<INT32>(0)) ? HurwitzZeta_Series::ZetaTemplate<ef_complex>(s, ef_complex(n))
                                      : efz::riemann_zeta(s);
}

e_float ef::hurwitz_zeta(const e_float& s, const e_float& a)
{
  if(ef::isint(a))
  {
    return ef::hurwitz_zeta(s, ef::to_int32(a));
  }
  else
  {
    return HurwitzZeta_Series::ZetaTemplate<e_float>(s, a);
  }
}

ef_complex efz::hurwitz_zeta(const ef_complex& s, const ef_complex& a)
{
  if(ef::isint(a))
  {
    return efz::hurwitz_zeta(s, ef::to_int32(a));
  }
  else
  {
    return HurwitzZeta_Series::ZetaTemplate<ef_complex>(s, a);
  }
}
