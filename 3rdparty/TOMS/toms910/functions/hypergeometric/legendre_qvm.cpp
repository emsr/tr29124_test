
#include <e_float/e_float.h>
#include <functions/constants/constants.h>
#include <functions/elementary/elementary.h>
#include <functions/gamma/gamma.h>
#include <functions/gamma/gamma_util.h>
#include <functions/gamma/polygamma.h>
#include <functions/hypergeometric/hypergeometric.h>
#include <functions/hypergeometric/hypergeometric_util.h>
#include <functions/hypergeometric/legendre.h>
#include <functions/hypergeometric/legendre_xvm.h>
#include <functions/integer/integer.h>

namespace Legendre_Series
{
  class LegendreQvm : public LegendreXvm
  {
  public:

    LegendreQvm() { }

    virtual ~LegendreQvm() { }

  private:

    virtual e_float AtIdenticallyZero        (const e_float& v, const INT32 m) const;
    virtual e_float AtReflectNegativeDegree  (const e_float& v, const INT32 m, const e_float& x) const;
    virtual e_float AtReflectNegativeOrder   (const e_float& v, const INT32 m, const e_float& x) const;
    virtual e_float AtReflectNegativeArgument(const e_float& v, const INT32 m, const e_float& x) const;
    virtual e_float AtOnePlus                (const e_float& v, const INT32 m, const e_float& x) const;

    virtual bool NeedsRecurOfM(const e_float& x) const { static_cast<void>(x); return true; }

    virtual e_float Legendre_nm(const INT32 n, const INT32 m, const e_float& x) const { return ef::legendre_q(n, m, x); } // NOCOVER_LINE
    virtual e_float Legendre_v (const e_float& v, const e_float& x) const             { return ef::legendre_q(v, x); };   // NOCOVER_LINE
  };
}

e_float Legendre_Series::LegendreQvm::AtIdenticallyZero(const e_float& v, const INT32 m) const
{
  // See Computation of Special Functions, Zhang & Jin, Equation 4.6.8, page 107.
  const e_float v_plus_m = v + m;

  const e_float v_plus_m_plus_one_over_two  = (v_plus_m + ef::one()) / static_cast<INT32>(2);
  const e_float v_minus_m_over_two_plus_one = (((v - m) / static_cast<INT32>(2)) + ef::one());

  const e_float g_factor = ef::gamma(v_plus_m_plus_one_over_two) / ef::gamma(v_minus_m_over_two_plus_one);

  const INT64 m_minus_one = static_cast<INT64>(m - static_cast<INT32>(1));

  return -((ef::pow2(m_minus_one) * ef::sqrt_pi()) * ef::sin(ef::pi_half() * v_plus_m)) * g_factor;
}

e_float Legendre_Series::LegendreQvm::AtReflectNegativeDegree(const e_float& v, const INT32 m, const e_float& x) const
{
  // Handle negative degree v < -1 with Computation of Special Functions,
  // Zhang & Jin, Equation 4.6.13, page 108.
  const e_float vv   = -v - ef::one();
  const e_float pi_v = ef::pi() * vv;

  return MyLegendre(vv, m, x) - ((ef::pi() * ef::cot(pi_v)) * ef::legendre_p(vv, m, x));
}

e_float Legendre_Series::LegendreQvm::AtReflectNegativeOrder(const e_float& v, const INT32 m, const e_float& x) const
{
  // Handle negative order m < 0 with Computation of Special Functions, Zhang & Jin,
  // Equation 4.6.10, page 108.

  const e_float v_plus_one = v + ef::one();
  const INT32   mm         = static_cast<INT32>(-m);

  const e_float Qp_vm = (MyLegendre(v, mm, x) * ef::gamma(v_plus_one - mm)) / ef::gamma(v_plus_one + mm);

  const bool b_negate = static_cast<INT32>(mm % static_cast<INT32>(2)) != static_cast<INT32>(0);

  return !b_negate ? Qp_vm : -Qp_vm;
}

e_float Legendre_Series::LegendreQvm::AtReflectNegativeArgument(const e_float& v, const INT32 m, const e_float& x) const
{
  const e_float xx = -x;

  e_float sin_pi_v_plus_pi_m;
  e_float cos_pi_v_plus_pi_m;
  ef::sincos(ef::pi() * (v + m), &sin_pi_v_plus_pi_m, &cos_pi_v_plus_pi_m);

  const e_float pvm_term = ef::legendre_p(v, m, xx) * sin_pi_v_plus_pi_m;
  const e_float qvm_term = ef::legendre_q(v, m, xx) * cos_pi_v_plus_pi_m;

  return -(ef::pi_half() * pvm_term) - qvm_term;
}

e_float Legendre_Series::LegendreQvm::AtOnePlus(const e_float& v, const INT32 m, const e_float& x) const
{
  // This series expansion is documented shown at Wolfram's Function Site.
  // http://functions.wolfram.com/HypergeometricFunctions/LegendreQ2General/06/01/02/

  e_float minus_v_pk     = -v;
  e_float minus_v_pkm    =  minus_v_pk + m;
  e_float v_plus_one_pk  =  v + ef::one();
  e_float v_plus_one_pkm =  v_plus_one_pk + m;
  INT32   m_plus_one_pk  =  static_cast<INT32>(m + static_cast<INT32>(1));

  const e_float v_minus_m_plus_one       = v_plus_one_pk - m;
  const e_float gamma_v_minus_m_plus_one = ef::gamma(v_minus_m_plus_one);
  const e_float v_plus_m_plus_one        = v_plus_one_pk + m;
  const e_float gamma_v_plus_m_plus_one  = ef::gamma(v_plus_m_plus_one);
  const e_float g_factor                 = gamma_v_plus_m_plus_one / gamma_v_minus_m_plus_one;

  const e_float one_minus_x              = ef::one() - x;
  const e_float one_minus_x_over_two     = one_minus_x / static_cast<INT32>(2);

  const e_float m_fact =  ef::factorial(static_cast<UINT32>(m));

  e_float pochham_minus_v_k        =  minus_v_pk;
  e_float pochham_minus_v_km       =  ef::pochhammer(minus_v_pk, m);
  e_float pochham_v_plus_one_k     =  v_plus_one_pk;
  e_float pochham_v_plus_one_km    =  ef::pochhammer(v_plus_one_pk, m);
  e_float pochham_m_plus_one_k     =  e_float(m_plus_one_pk);
  INT32       k_plus_one           =  static_cast<INT32>(1);
  e_float psi_k_plus_one           = -ef::euler_gamma();
  INT32       k_plus_m_plus_one    =  m_plus_one_pk;
  e_float psi_k_plus_m_plus_one    =  ef::poly_gamma(e_float(k_plus_m_plus_one));
  e_float k_fact                   =  ef::one();
  e_float k_plus_m_fact            =  m_fact;
  e_float m_minus_k_minus_one_fact =  m_fact / m;

  e_float one_minus_x_over_two_pow_k = ef::one();

  const bool b_negate_for_m_odd = static_cast<INT32>(m % static_cast<INT32>(2)) != static_cast<INT32>(0);
        bool b_negate_sum1_term = true;

  e_float sum1       =  m_minus_k_minus_one_fact;
  e_float sum2       =  ef::one();
  e_float sum3_part1 = ((pochham_minus_v_km * pochham_v_plus_one_km) * psi_k_plus_one) / m_fact;
  e_float sum3_part2 = psi_k_plus_m_plus_one / m_fact;

  pochham_minus_v_km    *=  minus_v_pkm;
  pochham_v_plus_one_km *=  v_plus_one_pkm;

  for(INT32 k = static_cast<INT32>(1); k < ef::max_iteration(); k++)
  {
    const INT32 k_plus_m  = static_cast<INT32>(k + m);

    k_fact        *= k;
    k_plus_m_fact *= k_plus_m;

    one_minus_x_over_two_pow_k *= one_minus_x_over_two;

    // Do the forward recursion of the psi functions using Psi(z + 1) = Psi(z) + 1/z.
    psi_k_plus_one        += (ef::one() / k_plus_one++);
    psi_k_plus_m_plus_one += (ef::one() / k_plus_m_plus_one++);

    const e_float pvm_pvp = pochham_minus_v_k * pochham_v_plus_one_k;

    // Increment the first sum.
    if(k < m)
    {
      const INT32 m_minus_k = static_cast<INT32>(m - k);

      m_minus_k_minus_one_fact /= m_minus_k;

      const e_float term1 = ((m_minus_k_minus_one_fact * pvm_pvp) * one_minus_x_over_two_pow_k) / k_fact;

      !b_negate_sum1_term ? sum1 += term1 : sum1 -= term1;

      b_negate_sum1_term = !b_negate_sum1_term;
    }

    // Increment the second sum.
    const e_float term2 = (pvm_pvp * one_minus_x_over_two_pow_k) / (pochham_m_plus_one_k * k_fact);

    // Increment the third sum.
    const e_float term3_scale = one_minus_x_over_two_pow_k / (k_fact * k_plus_m_fact);
    const e_float term3_part1 = ((pochham_minus_v_km * pochham_v_plus_one_km) * psi_k_plus_one) * term3_scale;
    const e_float term3_part2 = (psi_k_plus_m_plus_one * pvm_pvp) * term3_scale;

    const INT64 order_check2       = static_cast<INT64>(term2.order()       - sum2.order());
    const INT64 order_check3_part1 = static_cast<INT64>(term3_part1.order() - sum3_part1.order());
    const INT64 order_check3_part2 = static_cast<INT64>(term3_part2.order() - sum3_part2.order());

    const bool term2_is_negligible       = order_check2       < -ef::tol();
    const bool term3_part1_is_negligible = order_check3_part1 < -ef::tol();
    const bool term3_part2_is_negligible = order_check3_part2 < -ef::tol();

    const bool k_can_break_check = (k > static_cast<INT32>(20)) && (k > m);

    if(   k_can_break_check
       && term2_is_negligible
       && term3_part1_is_negligible
       && term3_part2_is_negligible
      )
    {
      break;
    }

    // Increment and multiply the pochhammer symbols. Note that this is done after the
    // evaluation of the terms. This is because the series summation begins with 1 and
    // the pochhammer symbols have already been initialized to the values at k = 1.
    pochham_minus_v_k     *=  ++minus_v_pk;
    pochham_v_plus_one_k  *=  ++v_plus_one_pk;
    pochham_m_plus_one_k  *=  ++m_plus_one_pk;
    pochham_minus_v_km    *=  ++minus_v_pkm;
    pochham_v_plus_one_km *=  ++v_plus_one_pkm;

    sum2       += term2;
    sum3_part1 += term3_part1;
    sum3_part2 += term3_part2;
  }

  static const e_float nine_tenths = ef::nine() / static_cast<INT32>(10);

  if(b_negate_for_m_odd)
  {
    sum3_part2 = -sum3_part2;
    sum2       = -sum2;
  }

  const e_float one_plus_x                = ef::one() + x;
  const e_float one_plus_x_over_two_pow_m = ef::pown(one_plus_x / static_cast<INT32>(2), static_cast<INT64>(m));

  sum3_part1 *= one_plus_x_over_two_pow_m;
  sum3_part2 *= g_factor;

  const e_float one_plus_x_over_one_minus_x = one_plus_x / one_minus_x;

  const e_float log_term_two = ef::small_arg(x) ? (ef::log1p1m2(x) * static_cast<INT32>(2))
                                                :  ef::log(one_plus_x_over_one_minus_x);

  const e_float log_psi_term = log_term_two - (ef::poly_gamma(v_plus_m_plus_one) + ef::poly_gamma(v_minus_m_plus_one));

  const e_float result_part1 = ef::pown(-one_plus_x_over_one_minus_x, static_cast<INT64>(m)) * sum1;
  const e_float result_part2 = ((g_factor * log_psi_term) * sum2) / m_fact;
  const e_float factor       = ef::pow(one_plus_x_over_one_minus_x, e_float(-m) / static_cast<INT32>(2)) / static_cast<INT32>(2);

  const e_float result = (result_part1 + result_part2) + (sum3_part1 + sum3_part2);

  return factor * result;
}

e_float ef::legendre_q(const e_float& v, const INT32 m, const e_float& x)
{
  static const Legendre_Series::LegendreQvm L_Qvm;
  return L_Qvm.MyLegendre(v, m, x);
}
