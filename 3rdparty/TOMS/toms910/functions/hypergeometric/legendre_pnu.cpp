
#include <e_float/e_float.h>
#include <functions/constants/constants.h>
#include <functions/elementary/elementary.h>
#include <functions/gamma/gamma.h>
#include <functions/hypergeometric/hypergeometric.h>
#include <functions/hypergeometric/legendre.h>
#include <functions/hypergeometric/legendre_xnu.h>

namespace Legendre_Series
{
  class LegendrePnu : public LegendreXnu
  {
  public:

    LegendrePnu() { }

    virtual ~LegendrePnu() { }

  private:

    virtual e_float AtReflectNegativeDegree  (const INT32 n, const e_float& u, const e_float& x) const;
    virtual e_float AtReflectNegativeArgument(const INT32 n, const e_float& u, const e_float& x) const;
    virtual e_float AtOnePlus                (const INT32 n, const e_float& u, const e_float& x) const;

    virtual bool NeedsReflectNegativeOrder(const e_float& u) const { static_cast<void>(u); return false; }

    virtual e_float Legendre_nm(const INT32 n, const INT32 m, const e_float& x) const { return ef::legendre_p(n, m, x); } // NOCOVER_LINE
  };
}

e_float Legendre_Series::LegendrePnu::AtReflectNegativeDegree(const INT32 n, const e_float& u, const e_float& x) const
{
  // Handle negative degree with n < -1 using
  // Abramowitz & Stegun 8.2.1, P_[-n-1, u](z) = P_[n, u](z).
  const INT32 nn = static_cast<INT32>(-n);

  return MyLegendre(static_cast<INT32>(nn - 1), u, x);
}

e_float Legendre_Series::LegendrePnu::AtReflectNegativeArgument(const INT32 n, const e_float& u, const e_float& x) const
{
  const e_float csc_term   = -ef::pi() * ef::csc(ef::pi() * u);
  const e_float gamma_term = ef::gamma(-u - n) * ef::gamma((n + static_cast<INT32>(1)) - u);

  return (csc_term * ef::legendre_p(n, -u, -x)) / gamma_term;
}

e_float Legendre_Series::LegendrePnu::AtOnePlus(const INT32 n, const e_float& u, const e_float& x) const
{
  // Implement the series expansion of Pnu(x) for values of x near +1.

  // This series expansion is documented at Wolfram's Function Site.
  // http://functions.wolfram.com/HypergeometricFunctions/LegendreP2General/06/01/02/
  // The expansion involves the the gamma function.

  const e_float one_minus_x                       = ef::one() - x;
  const e_float one_minus_x_half                  = one_minus_x / static_cast<INT32>(2);
        e_float one_minus_x_half_pow_k_over_kfact = ef::one();
        e_float k_munus_u_plus_one                = ef::one() - u;
        e_float gamma_k_minus_u_plus_one          = ef::gamma(k_munus_u_plus_one);
  const e_float nf                                = e_float(n);

  e_float pochham_minus_n    = -nf;
  e_float pochham_n_plus_one =  nf + ef::one();

  e_float minus_n_p          = pochham_minus_n;
  e_float n_plus_one_p       = pochham_n_plus_one;

  e_float sum  = ef::one() / gamma_k_minus_u_plus_one;

  for(INT32 k = static_cast<INT32>(1); k <= n; k++)
  {
    one_minus_x_half_pow_k_over_kfact *= one_minus_x_half;
    one_minus_x_half_pow_k_over_kfact /= k;

    gamma_k_minus_u_plus_one *= k_munus_u_plus_one++;

    sum += ((pochham_minus_n * pochham_n_plus_one) * one_minus_x_half_pow_k_over_kfact) / gamma_k_minus_u_plus_one;

    pochham_minus_n    *= ++minus_n_p;
    pochham_n_plus_one *= ++n_plus_one_p;
  }

  return ef::pow((ef::one() + x) / one_minus_x, u / static_cast<INT32>(2)) * sum;
}

e_float ef::legendre_p(const INT32 n, const e_float& u, const e_float& x)
{
  static const Legendre_Series::LegendrePnu L_Pnu;

  return L_Pnu.MyLegendre(n, u, x);
}
