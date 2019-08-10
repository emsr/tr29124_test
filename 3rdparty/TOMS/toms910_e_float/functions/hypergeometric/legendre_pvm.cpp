
#include <e_float/e_float.h>
#include <functions/constants/constants.h>
#include <functions/elementary/elementary.h>
#include <functions/gamma/gamma.h>
#include <functions/gamma/gamma_util.h>
#include <functions/hypergeometric/hypergeometric.h>
#include <functions/hypergeometric/hypergeometric_util.h>
#include <functions/hypergeometric/legendre.h>
#include <functions/hypergeometric/legendre_xvm.h>
#include <functions/integer/integer.h>

namespace Legendre_Series
{
  class LegendrePvm : public LegendreXvm
  {
  public:

    LegendrePvm() { }

    virtual ~LegendrePvm() { }

  private:

    virtual e_float AtIdenticallyZero        (const e_float& v, const INT32 m) const;
    virtual e_float AtReflectNegativeDegree  (const e_float& v, const INT32 m, const e_float& x) const;
    virtual e_float AtReflectNegativeOrder   (const e_float& v, const INT32 m, const e_float& x) const;
    virtual e_float AtReflectNegativeArgument(const e_float& v, const INT32 m, const e_float& x) const;
    virtual e_float AtOnePlus                (const e_float& v, const INT32 m, const e_float& x) const;

    virtual bool NeedsRecurOfM(const e_float& x) const { return ef::isneg(x); }

    virtual e_float Legendre_nm(const INT32 n, const INT32 m, const e_float& x) const { return ef::legendre_p(n, m, x); } // NOCOVER_LINE
    virtual e_float Legendre_v (const e_float& v, const e_float& x) const             { return ef::legendre_p(v, x); };
  };
}

e_float Legendre_Series::LegendrePvm::AtIdenticallyZero(const e_float& v, const INT32 m) const
{
  // See Computation of Special Functions, Zhang & Jin, Equation 4.6.6, page 107.
  const e_float v_plus_m = v + m;

  const e_float v_plus_m_plus_one_over_two  = (v_plus_m + ef::one()) / static_cast<INT32>(2);
  const e_float v_minus_m_over_two_plus_one = (((v - m) / static_cast<INT32>(2)) + ef::one());

  const e_float g_factor = ef::gamma(v_plus_m_plus_one_over_two) / ef::gamma(v_minus_m_over_two_plus_one);

  return ((ef::pow2(static_cast<INT64>(m)) / ef::sqrt_pi()) * ef::cos(ef::pi_half() * v_plus_m)) * g_factor;
}

e_float Legendre_Series::LegendrePvm::AtReflectNegativeDegree(const e_float& v, const INT32 m, const e_float& x) const
{
  // Handle negative degree v < -1 using
  // Abramowitz & Stegun 8.2.1, P_[-v-1, m](z) = P_[v, m](z).
  const e_float vv = -v - ef::one();

  return MyLegendre(vv, m, x);
}

e_float Legendre_Series::LegendrePvm::AtReflectNegativeOrder(const e_float& v, const INT32 m, const e_float& x) const
{
  // Handle negative order m < 0 with Computation of Special Functions, Zhang & Jin,
  // Equation 4.6.10, page 108.

  const e_float v_plus_one = v + ef::one();
  const INT32   mm         = static_cast<INT32>(-m);

  const e_float Lp_vm = (MyLegendre(v, mm, x) * ef::gamma(v_plus_one - mm)) / ef::gamma(v_plus_one + mm);

  const bool b_negate = static_cast<INT32>(mm % static_cast<INT32>(2)) != static_cast<INT32>(0);

  return !b_negate ? Lp_vm : -Lp_vm;
}

e_float Legendre_Series::LegendrePvm::AtReflectNegativeArgument(const e_float& v, const INT32 m, const e_float& x) const
{
  const e_float xx = -x;

  e_float sin_pi_v_plus_pi_m;
  e_float cos_pi_v_plus_pi_m;
  ef::sincos(ef::pi() * (v + m), &sin_pi_v_plus_pi_m, &cos_pi_v_plus_pi_m);

  const e_float pvm_term = ef::legendre_p(v, m, xx) * cos_pi_v_plus_pi_m;
  const e_float qvm_term = ef::legendre_q(v, m, xx) * sin_pi_v_plus_pi_m;

  return pvm_term - (qvm_term / ef::pi_half());
}

e_float Legendre_Series::LegendrePvm::AtOnePlus(const e_float& v, const INT32 m, const e_float& x) const
{
  // Implement the series expansion of Pvm(x) for values of x near +1.

  // This series expansion is documented at Wolfram's Function Site.
  // http://functions.wolfram.com/HypergeometricFunctions/LegendreP2General/06/01/02/
  // The expansion involves the the gamma function.

  const e_float one_minus_x_half                  = (ef::one() - x) / static_cast<INT32>(2);
        e_float one_minus_x_half_pow_k_over_kfact = one_minus_x_half;
        e_float mf(m);
  const e_float one_minus_xsq_pow_mhalf           = ef::pow(ef::one() - x * x, mf / static_cast<INT32>(2));
  const e_float minus_v                           = -v;
  const e_float v_plus_one                        =  v + ef::one();
  
  e_float pochham_m_minus_v         = mf + minus_v;
  e_float pochham_m_plus_v_plus_one = mf + v_plus_one;

  e_float m_minus_v_p               = pochham_m_minus_v;
  e_float m_plus_v_plus_one_p       = pochham_m_plus_v_plus_one;
  
  e_float k_plus_m_fact             = ef::factorial(static_cast<UINT32>(m));

  e_float sum    = ef::one() / k_plus_m_fact;
  k_plus_m_fact *= ++mf;

  sum += ((pochham_m_minus_v * pochham_m_plus_v_plus_one) * one_minus_x_half_pow_k_over_kfact) / k_plus_m_fact;

  for(INT32 k = static_cast<INT32>(2); k < ef::max_iteration(); k++)
  {
    one_minus_x_half_pow_k_over_kfact *= one_minus_x_half;
    one_minus_x_half_pow_k_over_kfact /= k;

    pochham_m_minus_v         *= ++m_minus_v_p;
    pochham_m_plus_v_plus_one *= ++m_plus_v_plus_one_p;

    k_plus_m_fact *= ++mf;

    const e_float term = ((pochham_m_minus_v * pochham_m_plus_v_plus_one) * one_minus_x_half_pow_k_over_kfact) / k_plus_m_fact;

    const INT64 order_check = static_cast<INT64>(term.order() - sum.order());

    if((k > static_cast<INT32>(20)) && (order_check < -ef::tol()))
    {
      break;
    }

    sum += term;
  }

  const e_float factor  =  (ef::pow2(static_cast<INT64>(-m)) * one_minus_xsq_pow_mhalf)
                         * (  ef::pochhammer(minus_v,    static_cast<UINT32>(m))
                            * ef::pochhammer(v_plus_one, static_cast<UINT32>(m)));

  return factor * sum;
}

e_float ef::legendre_p(const e_float& v, const INT32 m, const e_float& x)
{
  static const Legendre_Series::LegendrePvm L_Pvm;

  return L_Pvm.MyLegendre(v, m, x);
}
