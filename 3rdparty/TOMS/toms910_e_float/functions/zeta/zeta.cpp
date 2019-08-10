
#include <algorithm>

#include <e_float/e_float.h>
#include <functions/constants/constants.h>
#include <functions/elementary/elementary.h>
#include <functions/integer/integer.h>
#include <functions/integer/prime.h>
#include <functions/gamma/gamma.h>
#include <functions/tables/tables.h>
#include <functions/zeta/zeta.h>
#include <utility/util_power_j_pow_x.h>

namespace Zeta_Series
{
  static bool has_simple_form_for_zeta_n(const INT32 n)
  {
    // Zeta(n) has a simple form for all negative integer values of n as well
    // as for all positive, even values of n. The simple form is related to
    // a Bernoulli number. Therefore the simple form is only easy to calculate
    // as long as the range of n stays within the range of the tabulated values
    // of the Bernoulli numbers.

    // Check the range of the tabulated values of the Bernoulli numbers.
    static const std::size_t sz_A000367 = static_cast<std::size_t>(Tables::A000367().size());
    static const std::size_t sz_A002445 = static_cast<std::size_t>(Tables::A002445().size());

    static const INT32 b2n_table_max = static_cast<INT32>(2u * std::min(sz_A000367, sz_A002445));

    const bool is_neg  =  n < static_cast<INT32>(0);
    const bool is_even = (n % static_cast<INT32>(2)) == static_cast<INT32>(0);

    return (n < b2n_table_max) && ((!is_neg && is_even) || is_neg);
  }

  template<typename T> inline static T Reflection(const T& s)
  {
    using ef::gamma;
    using efz::gamma;
    using ef::pow;
    using efz::pow;
    using ef::sin;
    using efz::sin;
    using ef::riemann_zeta;
    using efz::riemann_zeta;

    const T two_pi_term = pow(T(ef::two_pi()), s) / ef::pi();
    const T chi         = (two_pi_term * sin(ef::pi_half() * s)) * gamma(ef::one() - s);

    return chi * riemann_zeta(ef::one() - s);
  }

  template<typename T> inline static T ZetaTemplate(const T& s)
  {
    using ef::abs;
    using efz::abs;
    using ef::imag;
    using efz::imag;
    using ef::pow;
    using efz::pow;
    using ef::real;
    using efz::real;

    if(ef::isint(s))
    {
      // Support pure-integer arguments according to certain conditions.
      const INT32 n = ef::to_int32(real(s));

      if(Zeta_Series::has_simple_form_for_zeta_n(n))
      {
        return ef::riemann_zeta(n);
      }
    }

    if(ef::isneg(s))
    {
      // Support negative arguments using the reflection formula.
      // Support arguments with a real part which is less than 1/2.
      return Zeta_Series::Reflection(s);
    }

    // The algorithms for calculating the Riemann zeta function below use calculations
    // of the integer j raised to the power s, or in other words j^s. The calculation of
    // j^s is accelerated using tables of the prime factors of j raised to the power s.
    // The calculation is furthermore accelerated by storing the necessary integer
    // logarithms in a static table.

    // Declare a map of prime factors raised to the power of the argument s.
    std::map<UINT32, T> n_pow_s_prime_factor_map;

    // Generate a list of the first 300 prime numbers.
    static std::deque<UINT32> prime_data;

    if(prime_data.empty())
    {
      ef::prime(static_cast<UINT32>(300u), prime_data);
    }

    static const std::vector<UINT32> primes(prime_data.begin(), prime_data.end());

    // Test the conditions for the expansion of the product of primes. Set up a
    // test for the product of primes. The expansion in the product of primes can
    // be used if the number of prime-power terms remains reasonably small, say
    // less than or equal to 300.
    static const double lg10_max_prime = ::log10(static_cast<double>(primes.back()));
    static const double td             = static_cast<double>(static_cast<INT32>(ef::tol()));
           const double sd             = ef::to_double(real(s));

    if(ef::iszero(imag(s)) && (static_cast<double>(td / sd) < lg10_max_prime))
    {
      // Perform the product of primes.
      T product = ef::one();

      for(std::size_t p = static_cast<std::size_t>(0u); p < primes.size(); p++)
      {
        const T prime_p_pow_s = Util::j_pow_x(primes[p], s, n_pow_s_prime_factor_map);

        product *= prime_p_pow_s / (prime_p_pow_s - ef::one());

        if(abs(prime_p_pow_s).order() > ef::tol())
        {
          break;
        }
      }

      return product;
    }
    else
    {
      if(abs(imag(s)) > (ef::million() + ef::hundred()))
      {
        // Return NaN if s has a large imaginary part.
        return std::numeric_limits<e_float>::quiet_NaN();
      }

      // Use the accelerated alternating converging series for Zeta as shown in:
      // http://numbers.computation.free.fr/Constants/Miscellaneous/zetaevaluations.html
      // taken from P. Borwein, "An Efficient Algorithm for the Riemann Zeta Function",
      // January 1995.

      // Compute the coefficients dk in a loop and calculate the zeta function sum
      // within the same loop on the fly.

      // Set up the factorials and powers for the calculation of the coefficients dk.
      // Note that j = n at this stage in the calculation. Also note that the value of
      // dn is equal to the value of d0 at the end of the loop.

      // Use N = (digits * 1.45) + {|imag(s)| * 1.1}
      static const double nd = static_cast<double>(std::numeric_limits<e_float>::digits10) * static_cast<double>(1.45);
             const double ni = static_cast<double>(static_cast<double>(1.10) * ::fabs(ef::to_double(imag(s))));

      const INT32 N        = static_cast<INT32>(static_cast<INT64>(static_cast<double>(nd + ni)));
            bool  neg_term = (N % static_cast<INT32>(2)) == static_cast<INT32>(0);

      e_float n_plus_j_minus_one_fact = ef::factorial(static_cast<UINT32>((N + N) - 1));
      e_float four_pow_j              = ef::pown(ef::four(), static_cast<INT64>(N));
      e_float n_minus_j_fact          = ef::one();
      e_float two_j_fact              = n_plus_j_minus_one_fact * static_cast<INT32>(static_cast<INT32>(2) * N);

      e_float dn = (n_plus_j_minus_one_fact * four_pow_j) / (n_minus_j_fact * two_j_fact);

      T jps = Util::j_pow_x(N, s, n_pow_s_prime_factor_map);

      T zs = (!neg_term ? dn : -dn) / jps;

      for(INT32 j = N - static_cast<INT32>(1); j >= static_cast<INT32>(0); j--)
      {
        const bool j_is_zero = (j == static_cast<INT32>(0));

        const INT32 two_jp1_two_j =   static_cast<INT32>((static_cast<INT32>(2) * j) + static_cast<INT32>(1))
                                    * static_cast<INT32> (static_cast<INT32>(2) * (!j_is_zero ? j : static_cast<INT32>(1)));

        n_plus_j_minus_one_fact /= static_cast<INT32>(N + j);
        four_pow_j              /= static_cast<INT32>(4);
        n_minus_j_fact          *= static_cast<INT32>(N - j);
        two_j_fact              /= two_jp1_two_j;

        dn += ((n_plus_j_minus_one_fact * four_pow_j) / (n_minus_j_fact * two_j_fact));

        if(!j_is_zero)
        {
          // Increment the zeta function sum.
          jps = Util::j_pow_x(static_cast<UINT32>(j), s, n_pow_s_prime_factor_map);

          neg_term = !neg_term;

          zs += (!neg_term ? dn : -dn) / jps;
        }
      }

      const T two_pow_one_minus_s = pow(T(2), ef::one() - s);

      return zs / (dn * (ef::one() - two_pow_one_minus_s));
    }
  }
}

e_float ef::riemann_zeta(const INT32 n)
{
  // Check if the result of the expansion will significantly differ from one.
  if(n > static_cast<INT32>(1))
  {
    static const double log10_of_2 = ::log10(static_cast<double>(2.0));
    static const double dtol       = static_cast<double>(static_cast<INT32>(ef::tol()));

    const double n_log10_of_2 = static_cast<double>(n) * log10_of_2;

    if(n_log10_of_2 > dtol)
    {
      // The result does not significantly differ from one.
      return ef::one();
    }
  }

  // Check if there is a simple form for Zeta(n).
  if(!Zeta_Series::has_simple_form_for_zeta_n(n))
  {
    // There is no simple form for Zeta(n). Do the e_float calculation.
    return riemann_zeta(e_float(n));
  }
  else
  {
    // There is a simple form for Zeta(n). Use the Bernoulli representation.
    if     (n == static_cast<INT32>(0)) { return -ef::half(); }
    else if(n == static_cast<INT32>(1)) { return std::numeric_limits<e_float>::infinity(); }
    else if(n == static_cast<INT32>(2)) { return  ef::pi_squared() / static_cast<INT32>(6); }
    else if(n == static_cast<INT32>(4)) { return (ef::pi_squared() * ef::pi_squared()) / static_cast<INT32>(90); }
    else
    {
      if(n < static_cast<INT32>(0))
      {
        const UINT32 two_n = static_cast<UINT32>(static_cast<INT32>(-n) + static_cast<INT32>(1));

        const bool is_even = (n % static_cast<INT32>(2)) == static_cast<INT32>(0);

        return is_even ? ef::zero() : -ef::bernoulli(two_n) / static_cast<INT32>(two_n);
      }
      else
      {
        const UINT32 two_n = static_cast<UINT32>(n);

        const e_float two_pi_pow_2n = ef::pown(ef::two_pi(), static_cast<INT64>(two_n));

        return ((two_pi_pow_2n * ef::fabs(ef::bernoulli(two_n))) / ef::factorial(two_n)) / static_cast<INT32>(2);
      }
    }
  }
}

e_float ef::riemann_zeta(const e_float& s)
{
  if(s.has_its_own_riemann_zeta())
  {
    return e_float::my_riemann_zeta(s);
  }
  else
  {
    return Zeta_Series::ZetaTemplate<e_float>(s);
  }
}

ef_complex efz::riemann_zeta(const ef_complex& s)
{
  return Zeta_Series::ZetaTemplate<ef_complex>(s);
}
