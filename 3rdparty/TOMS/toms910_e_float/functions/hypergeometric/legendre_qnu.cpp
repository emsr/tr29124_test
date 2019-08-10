
#include <e_float/e_float.h>
#include <functions/constants/constants.h>
#include <functions/elementary/elementary.h>
#include <functions/gamma/gamma.h>
#include <functions/hypergeometric/hypergeometric.h>
#include <functions/hypergeometric/legendre.h>
#include <functions/hypergeometric/legendre_xnu.h>

namespace Legendre_Series
{
  class LegendreQnu : public LegendreXnu
  {
  public:

    LegendreQnu() { }

    virtual ~LegendreQnu() { }

  private:

    virtual e_float AtReflectNegativeDegree  (const INT32 n, const e_float& u, const e_float& x) const;
    virtual e_float AtReflectNegativeOrder   (const INT32 n, const e_float& u, const e_float& x) const;
    virtual e_float AtReflectNegativeArgument(const INT32 n, const e_float& u, const e_float& x) const;
    virtual e_float AtOnePlus                (const INT32 n, const e_float& u, const e_float& x) const;

    virtual e_float Legendre_nm(const INT32 n, const INT32 m, const e_float& x) const { return ef::legendre_q(n, m, x); } // NOCOVER_LINE
  };
}

e_float Legendre_Series::LegendreQnu::AtReflectNegativeDegree(const INT32 n, const e_float& u, const e_float& x) const
{
  // Handle negative degree with n < 0 using an expression from Wolfram's Function Site
  // http://functions.wolfram.com/HypergeometricFunctions/LegendreQ2General/17/02/02/
  const INT32 nn = static_cast<INT32>(static_cast<INT32>(-n) - static_cast<INT32>(1));

  const bool b_negate_p_term = static_cast<INT32>(nn % static_cast<INT32>(2)) != static_cast<INT32>(0);

  const e_float pi_u = ef::pi() * u;
  const e_float pi_n = ef::pi() * nn;

  const e_float cos_pi_u = ef::cos(pi_u);

  const e_float sin_plus_pi_u_plus_pi_n  = ef::sin(pi_u + pi_n);
  const e_float sin_plus_pi_u_minus_pi_n = ef::sin(pi_u - pi_n);

        e_float p_term = (ef::pi() * cos_pi_u) * ef::legendre_p(nn, u, x);
  const e_float q_term = sin_plus_pi_u_plus_pi_n * MyLegendre(nn, u, x);

  if(b_negate_p_term) { p_term = -p_term; }

  return (p_term - q_term) / sin_plus_pi_u_minus_pi_n;
}

e_float Legendre_Series::LegendreQnu::AtReflectNegativeOrder(const INT32 n, const e_float& u, const e_float& x) const
{
  // Reflect negative order using an expression from Wolfram's Function Site.
  // http://functions.wolfram.com/HypergeometricFunctions/LegendreQ2General/17/02/02/

  const e_float uu = -u;

  const e_float pi_u = ef::pi() * uu;

  e_float sin_pi_u;
  e_float cos_pi_u;
  ef::sincos(pi_u, &sin_pi_u, &cos_pi_u);

  const e_float gamma_plus_v_minus_u_plus_one = ef::gamma((n - uu) + ef::one());
  const e_float gamma_plus_v_plus_u_plus_one  = ef::gamma((n + uu) + ef::one());

  const e_float term_p = ((ef::pi() * sin_pi_u) * ef::legendre_p(n, uu, x)) / static_cast<INT32>(2);
  const e_float term_q = cos_pi_u * MyLegendre(n, uu, x);

  const e_float g_factor = gamma_plus_v_minus_u_plus_one / gamma_plus_v_plus_u_plus_one;

  return g_factor * (term_p + term_q);
}

e_float Legendre_Series::LegendreQnu::AtReflectNegativeArgument(const INT32 n, const e_float& u, const e_float& x) const
{
  const e_float csc_term   = ef::pi() * ef::csc(ef::pi() * u);
  const e_float gamma_term = ef::gamma(-u - n) * ef::gamma((n + static_cast<INT32>(1)) - u);

  return (csc_term / gamma_term) * ef::legendre_q(n, -u, -x);
}

e_float Legendre_Series::LegendreQnu::AtOnePlus(const INT32 n, const e_float& u, const e_float& x) const
{
  // See Computation of Special Functions, Zhang & Jin, Equation 4.6.5, page 105.
  const e_float v_pi_plus_m_pi = (n + u) * ef::pi();
  e_float sin_n_pi_plus_u_pi;
  e_float cos_n_pi_plus_u_pi;

  ef::sincos(v_pi_plus_m_pi, &sin_n_pi_plus_u_pi, &cos_n_pi_plus_u_pi);

  return ef::pi_half() * (((ef::legendre_p(n, u, x) * cos_n_pi_plus_u_pi) - ef::legendre_p(n, u, -x)) / sin_n_pi_plus_u_pi);
}

e_float ef::legendre_q(const INT32 n, const e_float& u, const e_float& x)
{
  static const Legendre_Series::LegendreQnu L_Qnu;
  return L_Qnu.MyLegendre(n, u, x);
}
