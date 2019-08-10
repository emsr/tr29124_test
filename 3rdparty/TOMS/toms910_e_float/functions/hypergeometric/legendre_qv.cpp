
#include <e_float/e_float.h>
#include <functions/constants/constants.h>
#include <functions/elementary/elementary.h>
#include <functions/gamma/gamma.h>
#include <functions/gamma/polygamma.h>
#include <functions/hypergeometric/hypergeometric.h>
#include <functions/hypergeometric/legendre.h>
#include <functions/hypergeometric/legendre_xv.h>
#include <functions/polynomials/polynomials.h>

namespace Legendre_Series
{
  class LegendreQv : public LegendreXv
  {
  public:

    LegendreQv() { }

    virtual ~LegendreQv() { }

  private:

    void TwoSumExpansion(const e_float& v, const e_float& x, e_float& sum1, e_float& sum2) const;

  private:

    virtual e_float AtIdenticallyZero      (const e_float& v) const;
    virtual e_float AtReflectNegativeDegree(const e_float& v, const e_float& x) const;
    virtual e_float AtOnePlus              (const e_float& v, const e_float& x) const;
    virtual e_float AtOneMinus             (const e_float& v, const e_float& x) const;

    virtual e_float Legendre_n(const INT32 n, const e_float& x) const { return ef::legendre_q(n, x); } // NOCOVER_LINE
  };
}

void Legendre_Series::LegendreQv::TwoSumExpansion(const e_float& v, const e_float& x, e_float& sum1, e_float& sum2) const
{
  // Compute the two main parts of the series expansion for Qv.
  // These series expansions are documented at Wolfram's Function Site.
  // http://functions.wolfram.com/HypergeometricFunctions/LegendreQGeneral/06/01/02/
  // http://functions.wolfram.com/HypergeometricFunctions/LegendreQGeneral/06/01/03/
  // The expansions involves pochhammer symbols and the digamma (psi) function.

  e_float v_minus              = -v;
  e_float v_plus1              =  v + ef::one();
  e_float pochham_v_minus      =  v_minus;
  e_float pochham_v_plus1      =  v_plus1;
  e_float x_pow_k_over_k_fact2 =  ef::one();
  INT32   k_plus1              =  static_cast<INT32>(1);
  e_float psi_k_plus1          = -ef::euler_gamma();

  sum1 = ef::one();
  sum2 = psi_k_plus1;

  for(INT32 k = static_cast<INT32>(1); k < ef::max_iteration(); k++)
  {
    // Multiply with [(1 + x) / 2] and divide with the next term of (k!)^2.
    x_pow_k_over_k_fact2 *= x;
    x_pow_k_over_k_fact2 /= static_cast<INT32>(k * k);

    // Do the forward recursion of the psi functions using Psi(z + 1) = Psi(z) + 1/z.
    psi_k_plus1 += (ef::one() / k_plus1++);

    const e_float pvm_pvp = pochham_v_minus * pochham_v_plus1;
    const e_float term1   =  pvm_pvp * x_pow_k_over_k_fact2;
    const e_float term2   = (pvm_pvp * psi_k_plus1) * x_pow_k_over_k_fact2;
    
    if((term1.order() < -ef::tol()) && (term2.order() < -ef::tol()))
    {
      break;
    }

    // Increment and multiply the pochhammer symbols. Note that this is done after the
    // evaluation of the terms. This is because the series summation begins with 1 and
    // the pochhammer symbols have already been initialized to the values at k = 1.
    pochham_v_minus *= ++v_minus;
    pochham_v_plus1 *= ++v_plus1;

    sum1 += term1;
    sum2 += term2;
  }
}

e_float Legendre_Series::LegendreQv::AtIdenticallyZero(const e_float& v) const
{
  const e_float v_half = v / static_cast<INT32>(2);

  const e_float g_factor = ef::gamma(v_half + ef::half()) / ef::gamma(v_half + ef::one());

  return -(((g_factor * ef::sin(ef::pi_half() * v)) * ef::sqrt_pi()) / static_cast<INT32>(2));
}

e_float Legendre_Series::LegendreQv::AtReflectNegativeDegree(const e_float& v, const e_float& x) const
{
  // Handle negative degree with v < -1 using
  // Abramowitz & Stegun 8.2.2, Q_[-v-1](z) = - {pi * cot(v * pi) * P_[v](z)} + Q_[v](z)
  // with u = 0.
  const e_float vv = -v - ef::one();

  return -((ef::pi() * ef::cot(vv * ef::pi())) * ef::legendre_p(vv, x)) + MyLegendre(vv, x);
}

e_float Legendre_Series::LegendreQv::AtOnePlus(const e_float& v, const e_float& x) const
{
  // This series expansion is documented at Wolfram's Function Site.
  // http://functions.wolfram.com/HypergeometricFunctions/LegendreQGeneral/06/01/02/
  // The expansion involves pochhammer symbols and the digamma (psi) function,
  // as well as logarithmic functions. The core of the expansion is carried out with
  // the subroutine CoreExpansion above.

  const e_float one_minus_x = ef::one() - x;

  e_float sum1;
  e_float sum2;

  TwoSumExpansion(v, one_minus_x / static_cast<INT32>(2), sum1, sum2);

  const e_float log_term = ef::small_arg(x) ?  ef::log1p1m2(x)
                                            : (ef::log((ef::one() + x) / (one_minus_x)) / static_cast<INT32>(2));

  const e_float factor1 = log_term - ef::poly_gamma(v + ef::one());

  return (sum1 * factor1) + sum2;
}

e_float Legendre_Series::LegendreQv::AtOneMinus(const e_float& v, const e_float& x) const
{
  // This series expansion is documented shown at Wolfram's Function Site.
  // http://functions.wolfram.com/HypergeometricFunctions/LegendreQGeneral/06/01/03/
  // The expansion involves pochhammer symbols and the digamma (psi) function,
  // as well as trigonometric and logarithmic functions. The core of the expansion
  // is carried out with the subroutine CoreExpansion above.

  const e_float one_plus_x = ef::one() + x;

  e_float sum1;
  e_float sum2;

  TwoSumExpansion(v, one_plus_x / static_cast<INT32>(2), sum1, sum2);

  const e_float pi_v = ef::pi() * v;
        e_float cos_pi_v;
        e_float sin_pi_v;
        
  ef::sincos(pi_v, &sin_pi_v, &cos_pi_v);
  
  const e_float pi_half_csc_pi_v = ef::pi_half() / sin_pi_v;

  const e_float log_term = ef::small_arg(x) ? ef::log1p1m2(x)
                                            : ef::log(one_plus_x / (ef::one() - x)) / static_cast<INT32>(2);

  const e_float factor1 =   (cos_pi_v * ((ef::poly_gamma(-v) - (pi_half_csc_pi_v * cos_pi_v)) +  log_term))
                          - pi_half_csc_pi_v;

  return (sum1 * factor1) - (sum2 * cos_pi_v);
}

e_float ef::legendre_q(const e_float& v, const e_float& x)
{
  static const Legendre_Series::LegendreQv L_Qv;

  return L_Qv.MyLegendre(v, x);
}
