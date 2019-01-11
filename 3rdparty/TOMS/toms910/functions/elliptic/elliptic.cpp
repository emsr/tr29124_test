
#include <deque>
#include <numeric>

#include <e_float/e_float.h>
#include <functions/constants/constants.h>
#include <functions/elementary/elementary.h>
#include <functions/elliptic/elliptic.h>

namespace Elliptic_Series
{
  static void AGM(const e_float& phi,
                  const e_float& m,
                        e_float& Fpm,
                        e_float& Km,
                        e_float* const pEm,
                        e_float* const pEpm);
}

static void Elliptic_Series::AGM(const e_float& phi,
                                 const e_float& m,
                                       e_float& Fpm,
                                       e_float& Km,
                                       e_float* const pEm  = static_cast<e_float* const>(0u),
                                       e_float* const pEpm = static_cast<e_float* const>(0u))
{
  // Use the AGM algorithm as described in Computation of Special Functions,
  // Zhang & Jin, 18.3.2, pages 663-665. The implementation is based on the
  // sample code therein. However, the Mathematica argument convention with
  // (k^2 --> m) is used, as described in Stephen Wolfram's Mathematica Book,
  // 4th Ed., Ch. 3.2.11, Page 773.
  
  // Make use of the following properties:
  // F(z + pi*k | m) = F(z | m) + 2k pi K(m)
  // E(z + pi*k | m) = E(z | m) + 2k pi E(m)

  // as well as:
  // F(-z | m) = -F(z | m)
  // E(-z | m) = -E(z | m)

  // The calculations which are needed for EllipticE(...) are only performed if
  // the results from these will actually be used, in other words only if non-zero
  // pointers pEm or pEpm have been supplied to this subroutine.

  // Note that there is special handling for the angular argument phi if this
  // argument is equal to pi/2.

  const bool phi_is_pi_half = (phi == ef::pi_half());

  if(ef::iszero(m))
  {
    Fpm = phi;
    Km  = ef::pi_half();

    if(pEpm != static_cast<e_float* const>(0u)) { *pEpm = phi; }
    if(pEm  != static_cast<e_float* const>(0u)) { *pEm  = ef::pi_half(); }
  }
  else if(ef::isone(m))
  {
    if(pEm != static_cast<e_float* const>(0u)) { *pEm = ef::one(); }

    Km = std::numeric_limits<e_float>::quiet_NaN();

    const e_float sp = ef::sin(phi);

    Fpm = phi_is_pi_half ? std::numeric_limits<e_float>::quiet_NaN()
                         : ef::log((ef::one() + sp) / (ef::one() - sp)) / static_cast<INT32>(2);

    if(pEpm != static_cast<e_float* const>(0u)) { *pEpm = phi_is_pi_half ? ef::one() : sp; }
  }
  else
  {
    e_float a0    = ef::one();
    e_float b0    = ef::sqrt(ef::one() - m); // Mathematica argument convention
    e_float phi_n = phi;
    e_float p2    = ef::one();

    e_float an;

    const bool    m_neg = ef::isneg(m);
    const e_float mk    = ef::sqrt(ef::fabs(m));

    std::deque<e_float> cn          (static_cast<std::size_t>(1u), !m_neg ? mk : -mk);
    std::deque<e_float> two_pow_n_cn(static_cast<std::size_t>(1u), mk);
    std::deque<e_float> sin_phi_n   (static_cast<std::size_t>(1u), ef::sin(phi));
    
    const bool has_e =    (pEm  != static_cast<e_float* const>(0u))
                       || (pEpm != static_cast<e_float* const>(0u));

    for(INT32 n = static_cast<INT32>(1); n < static_cast<INT32>(64); n++)
    {
      an = (a0 + b0) / static_cast<INT32>(2);

      p2 *= static_cast<INT32>(2);

      if(!phi_is_pi_half) { phi_n += ef::atan((b0 / a0) * ef::tan(phi_n)); }

      const e_float cn_term = (a0 - b0) / static_cast<INT32>(2);

      if(has_e)
      {
        cn.push_back(cn_term);
        two_pow_n_cn.push_back(cn.back() * p2);
        sin_phi_n.push_back(!phi_is_pi_half ? ef::sin(phi_n) : ef::zero());
      }

      if(cn_term.order() < -static_cast<INT64>(ef::tol() / static_cast<INT64>(2)))
      {
        break;
      }

      b0 = ef::sqrt(a0 * b0);
      a0 = an;

      if(!phi_is_pi_half)
      {
        phi_n += ef::pi() * static_cast<INT32>(static_cast<INT64>(ef::to_double((phi_n / ef::pi()) + ef::half())));
      }
    }

    const e_float one_over_an = ef::one() / an;

    Fpm = phi_n * one_over_an;

    if(!phi_is_pi_half) { Fpm /= p2; }

    Km = ef::pi_half() * one_over_an;

    if(has_e)
    {
      const e_float one_minus_cn_2ncn_inner_prod_half = ef::one() - (std::inner_product(cn.begin(),
                                                                                        cn.end(),
                                                                                        two_pow_n_cn.begin(),
                                                                                        ef::zero()) / static_cast<INT32>(2));

      if(pEm != static_cast<e_float* const>(0u))
      {
        *pEm = Km * one_minus_cn_2ncn_inner_prod_half;
      }

      if(pEpm != static_cast<e_float* const>(0u))
      {
        *pEpm =   (Fpm * one_minus_cn_2ncn_inner_prod_half)
                +  std::inner_product(sin_phi_n.begin() + static_cast<std::size_t>(1u),
                                      sin_phi_n.end(),
                                      cn.begin() + static_cast<std::size_t>(1u),
                                      ef::zero());
      }
    }
  }
}

e_float ef::ellint_2(const e_float& m, const e_float& phi)
{
  if(ef::fabs(m) > ef::one())
  {
    return std::numeric_limits<e_float>::quiet_NaN();
  }
  else
  {
    if(ef::isneg(phi))
    {
      return -ellint_2(-phi, m);
    }
    else
    {
      e_float k_pi       = ef::integer_part(phi / ef::pi());
      e_float phi_scaled = phi - (k_pi * ef::pi());
      bool    b_neg      = false;

      if(phi_scaled > ef::pi_half())
      {
        ++k_pi;
        phi_scaled = -(phi_scaled - ef::pi());
        b_neg      = true;
      }

      e_float Fpm, Km, Em, Epm;

      Elliptic_Series::AGM(phi_scaled, m, Fpm, Km, &Em, &Epm);

      if(b_neg)
      {
        Epm = -Epm;
      }

      return Epm + ((k_pi * Em) * static_cast<INT32>(2));
    }
  }
}

e_float ef::comp_ellint_2(const e_float& m)
{
  if(ef::fabs(m) > ef::one())
  {
    return std::numeric_limits<e_float>::quiet_NaN();
  }
  else
  {
    e_float Fpm, Km, Em;

    Elliptic_Series::AGM(ef::zero(), m, Fpm, Km, &Em);

    return Em;
  }
}

e_float ef::ellint_1(const e_float& m, const e_float& phi)
{
  if(ef::fabs(m) > ef::one())
  {
    return std::numeric_limits<e_float>::quiet_NaN();
  }
  else
  {
    if(ef::isneg(phi))
    {
      return -ellint_1(-phi, m);
    }
    else
    {
      e_float k_pi       = ef::integer_part(phi / ef::pi());
      e_float phi_scaled = phi - (k_pi * ef::pi());
      bool    b_neg      = false;

      if(phi_scaled > ef::pi_half())
      {
        ++k_pi;
        phi_scaled = -(phi_scaled - ef::pi());
        b_neg      = true;
      }

      e_float Fpm, Km;

      Elliptic_Series::AGM(phi_scaled, m, Fpm, Km);

      if(b_neg)
      {
        Fpm = -Fpm;
      }

      return Fpm + ((k_pi * Km) * static_cast<INT32>(2));
    }
  }
}

e_float ef::comp_ellint_1(const e_float& m)
{
  if(ef::fabs(m) > ef::one())
  {
    return std::numeric_limits<e_float>::quiet_NaN();
  }
  else
  {
    e_float Fpm, Km;

    Elliptic_Series::AGM(ef::zero(), m, Fpm, Km);

    return Km;
  }
}
