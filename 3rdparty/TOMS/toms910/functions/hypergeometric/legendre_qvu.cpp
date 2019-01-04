
#include <e_float/e_float.h>
#include <functions/constants/constants.h>
#include <functions/elementary/elementary.h>
#include <functions/gamma/gamma.h>
#include <functions/gamma/gamma_util.h>
#include <functions/hypergeometric/hypergeometric.h>
#include <functions/hypergeometric/hypergeometric_util.h>
#include <functions/hypergeometric/legendre.h>
#include <functions/hypergeometric/legendre_xvu.h>
#include <functions/integer/integer.h>

namespace Legendre_Series
{
  class LegendreQvu : public LegendreXvu
  {
  public:

    LegendreQvu() { }

    virtual ~LegendreQvu() { }

  private:

    virtual e_float AtReflectNegativeDegree  (const e_float& v, const e_float& u, const e_float& x) const;
    virtual e_float AtReflectNegativeOrder   (const e_float& v, const e_float& u, const e_float& x) const;
    virtual e_float AtReflectNegativeArgument(const e_float& v, const e_float& u, const e_float& x) const;
    virtual e_float AtOnePlus                (const e_float& v, const e_float& u, const e_float& x) const;

    virtual e_float Legendre_vm(const e_float& v, const INT32    m, const e_float& x) const { return ef::legendre_q(v, m, x); } // NOCOVER_LINE
    virtual e_float Legendre_nu(const INT32    n, const e_float& u, const e_float& x) const { return ef::legendre_q(n, u, x); } // NOCOVER_LINE
  };
}

e_float Legendre_Series::LegendreQvu::AtReflectNegativeDegree(const e_float& v, const e_float& u, const e_float& x) const
{
  // Handle negative degree with v < -1 using an expression from Wolfram's Function Site
  // http://functions.wolfram.com/HypergeometricFunctions/LegendreQ2General/17/02/02/
  const e_float vv = -v - ef::one();

  const e_float pi_u = ef::pi() * u;
  const e_float pi_v = ef::pi() * vv;

  const e_float cos_pi_u = ef::cos(pi_u);
  const e_float cos_pi_v = ef::cos(pi_v);

  const e_float sin_plus_pi_u_plus_pi_v  = ef::sin(pi_u + pi_v);
  const e_float sin_plus_pi_u_minus_pi_v = ef::sin(pi_u - pi_v);

  const e_float p_term = (ef::pi() * (cos_pi_u * cos_pi_v)) * ef::legendre_p(vv, u, x);
  const e_float q_term = sin_plus_pi_u_plus_pi_v * MyLegendre(vv, u, x);

  return (p_term - q_term) / sin_plus_pi_u_minus_pi_v;
}

e_float Legendre_Series::LegendreQvu::AtReflectNegativeOrder(const e_float& v, const e_float& u, const e_float& x) const
{
  // Reflect negative order using an expression from Wolfram's Function Site.
  // http://functions.wolfram.com/HypergeometricFunctions/LegendreQ2General/17/02/02/

  const e_float uu = -u;

  const e_float pi_u = ef::pi() * uu;

  e_float sin_pi_u;
  e_float cos_pi_u;
  ef::sincos(pi_u, &sin_pi_u, &cos_pi_u);

  const e_float gamma_plus_v_minus_u_plus_one = ef::gamma((v - uu) + ef::one());
  const e_float gamma_plus_v_plus_u_plus_one  = ef::gamma((v + uu) + ef::one());

  const e_float term_p = ((ef::pi() * sin_pi_u) * ef::legendre_p(v, uu, x)) / static_cast<INT32>(2);
  const e_float term_q = cos_pi_u * MyLegendre(v, uu, x);

  const e_float g_factor = gamma_plus_v_minus_u_plus_one / gamma_plus_v_plus_u_plus_one;

  return g_factor * (term_p + term_q);
}

e_float Legendre_Series::LegendreQvu::AtReflectNegativeArgument(const e_float& v, const e_float& u, const e_float& x) const
{
  const e_float xx = -x;

  e_float sin_pi_v_plus_pi_u;
  e_float cos_pi_v_plus_pi_u;
  ef::sincos(ef::pi() * (v + u), &sin_pi_v_plus_pi_u, &cos_pi_v_plus_pi_u);

  const e_float pvu_term = ef::legendre_p(v, u, xx) * sin_pi_v_plus_pi_u;
  const e_float qvu_term = ef::legendre_q(v, u, xx) * cos_pi_v_plus_pi_u;

  return -(ef::pi_half() * pvu_term) - qvu_term;
}

e_float Legendre_Series::LegendreQvu::AtOnePlus(const e_float& v, const e_float& u, const e_float& x) const
{
  // This series expansion is documented at Wolfram's Function Site.
  // http://functions.wolfram.com/HypergeometricFunctions/LegendreQ2General/06/01/02/
  // The expansion involves the regularized hypergeometric function 2F1.

  e_float sin_u_pi;
  e_float cos_u_pi;
  ef::sincos(u * ef::pi(), &sin_u_pi, &cos_u_pi);

  const e_float one_minus_x          = ef::one() - x;
  const e_float one_plus_x           = ef::one() + x;
  const e_float u_half               = u / static_cast<INT32>(2);
  const e_float one_minus_x_over_two = one_minus_x / static_cast<INT32>(2);

  const e_float one_plus_x_over_one_minus_x_pow_u_half = ef::pow(one_plus_x / one_minus_x, u_half);

  const e_float minus_v    = -v;
  const e_float v_plus_one =  v + ef::one();

  const e_float H2F1_1 = ef::hyperg_2f1_reg(minus_v, v_plus_one, ef::one() - u, one_minus_x_over_two);
  const e_float H2F1_2 = ef::hyperg_2f1_reg(minus_v, v_plus_one, ef::one() + u, one_minus_x_over_two);

  const e_float term1 = (H2F1_1 * one_plus_x_over_one_minus_x_pow_u_half) * cos_u_pi;
  const e_float term2 = (H2F1_2 / one_plus_x_over_one_minus_x_pow_u_half) * ef::pochhammer((v - u) + ef::one(), u * static_cast<INT32>(2));

  return (ef::pi_half() / sin_u_pi) * (term1 - term2);
}

e_float ef::legendre_q(const e_float& v, const e_float& u, const e_float& x)
{
  static const Legendre_Series::LegendreQvu L_Qvu;
  return L_Qvu.MyLegendre(v, u, x);
}
