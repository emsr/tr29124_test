
#include <e_float/e_float.h>
#include <functions/constants/constants.h>
#include <functions/elementary/elementary.h>
#include <functions/gamma/gamma.h>
#include <functions/gamma/gamma_util.h>
#include <functions/hypergeometric/hypergeometric.h>
#include <functions/hypergeometric/legendre.h>
#include <functions/hypergeometric/legendre_xvu.h>
#include <functions/integer/integer.h>

namespace Legendre_Series
{
  class LegendrePvu : public LegendreXvu
  {
  public:

    LegendrePvu() { }

    virtual ~LegendrePvu() { }

  private:

    virtual e_float AtReflectNegativeDegree  (const e_float& v, const e_float& u, const e_float& x) const;
    virtual e_float AtReflectNegativeArgument(const e_float& v, const e_float& u, const e_float& x) const;
    virtual e_float AtOnePlus                (const e_float& v, const e_float& u, const e_float& x) const;

    virtual bool NeedsReflectNegativeOrder(const e_float& u) const { static_cast<void>(u); return false; }

    virtual e_float Legendre_vm(const e_float& v, const INT32    m, const e_float& x) const { return ef::legendre_p(v, m, x); } // NOCOVER_LINE
    virtual e_float Legendre_nu(const INT32    n, const e_float& u, const e_float& x) const { return ef::legendre_p(n, u, x); } // NOCOVER_LINE
  };
}

e_float Legendre_Series::LegendrePvu::AtReflectNegativeDegree(const e_float& v, const e_float& u, const e_float& x) const
{
  // Handle negative degree with v < -1 using
  // Abramowitz & Stegun 8.2.1, P_[-v-1, u](z) = P_[v, u](z).
  const e_float vv = -v - ef::one();

  return MyLegendre(vv, u, x);
}

e_float Legendre_Series::LegendrePvu::AtReflectNegativeArgument(const e_float& v, const e_float& u, const e_float& x) const
{
  const e_float xx = -x;

  e_float sin_pi_v_plus_pi_u;
  e_float cos_pi_v_plus_pi_u;
  ef::sincos(ef::pi() * (v + u), &sin_pi_v_plus_pi_u, &cos_pi_v_plus_pi_u);

  const e_float pvu_term = ef::legendre_p(v, u, xx) * cos_pi_v_plus_pi_u;
  const e_float qvu_term = ef::legendre_q(v, u, xx) * sin_pi_v_plus_pi_u;

  return pvu_term - (qvu_term / ef::pi_half());
}

e_float Legendre_Series::LegendrePvu::AtOnePlus(const e_float& v, const e_float& u, const e_float& x) const
{
  // This series expansion is documented shown at Wolfram's Function Site.
  // http://functions.wolfram.com/HypergeometricFunctions/LegendreP2General/06/01/02/
  // The expansion involves the regularized hypergeometric function 2F1.

  const e_float u_half      = u / static_cast<INT32>(2);
  const e_float one_minus_x = ef::one() - x;

  return   (ef::pow(ef::one() + x, u_half) / ef::pow(one_minus_x, u_half))
         *  ef::hyperg_2f1_reg(-v,
                                ef::one() + v,
                                ef::one() - u,
                                one_minus_x / static_cast<INT32>(2));
}

e_float ef::legendre_p(const e_float& v, const e_float& u, const e_float& x)
{
  static const Legendre_Series::LegendrePvu L_Pvu;

  return L_Pvu.MyLegendre(v, u, x);
}
