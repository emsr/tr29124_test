
#include <e_float/e_float.h>
#include <functions/constants/constants.h>
#include <functions/elementary/elementary.h>
#include <functions/gamma/gamma.h>
#include <functions/gamma/gamma_util.h>
#include <functions/hypergeometric/hypergeometric.h>
#include <functions/hypergeometric/legendre.h>
#include <functions/hypergeometric/legendre_xv.h>
#include <functions/polynomials/polynomials.h>
#include <functions/tables/tables.h>

namespace Legendre_Series
{
  class LegendrePv : public LegendreXv
  {
  public:

    LegendrePv() { }

    virtual ~LegendrePv() { }

  private:

    virtual e_float AtIdenticallyZero      (const e_float& v) const;
    virtual e_float AtReflectNegativeDegree(const e_float& v, const e_float& x) const;
    virtual e_float AtOnePlus              (const e_float& v, const e_float& x) const;
    virtual e_float AtOneMinus             (const e_float& v, const e_float& x) const;

    virtual e_float Legendre_n(const INT32 n, const e_float& x) const { return ef::legendre_p(n, x); } // NOCOVER_LINE
  };
}

e_float Legendre_Series::LegendrePv::AtIdenticallyZero(const e_float& v) const
{
  const e_float v_half = v / static_cast<INT32>(2);

  const e_float g_factor = ef::gamma(v_half + ef::half()) / ef::gamma(v_half + ef::one());

  return (g_factor * ef::cos(ef::pi_half() * v)) / ef::sqrt_pi();
}

e_float Legendre_Series::LegendrePv::AtReflectNegativeDegree(const e_float& v, const e_float& x) const
{
  // Handle negative degree with v < -1 using
  // Abramowitz & Stegun 8.2.1, P_[-v-1](z) = P_[v](z).
  const e_float vv = -v;
  return MyLegendre(vv - ef::one(), x);
}

e_float Legendre_Series::LegendrePv::AtOnePlus(const e_float& v, const e_float& x) const
{
  // Implement the series expansion of Pv(x) for values of x near +1.

  // See Computation of Special Functions, Zhang & Jin, 4.6.15, page 108.
  // Also see Abramowitz & Stegun 8.1.2, page 332 for u = 0.

  return ef::hyperg_2f1(-v,
                         ef::one() + v,
                         ef::one(),
                        (ef::one() - x) / static_cast<INT32>(2));
}

e_float Legendre_Series::LegendrePv::AtOneMinus(const e_float& v, const e_float& x) const
{
  // Implement the series expansion of Pv(x) for values of x near -1.

  // See Computation of Special Functions, Zhang & Jin, 4.6.15, page 108.
  // The series representation can also be calculated from Mathematica.
  // For example, the Mathematica Version 4.1 command for obtaining the
  // series up to order 4 is:
  // In[1]:= Series[LegendreP[v, x], {x, -1, 4}]

  // Calculate the initial values of the digamma (psi) functions.
  e_float psi_v;
  e_float psi_minus_v;

  GammaUtil::DiGammaOfPlusXMinusX(v,
                                  psi_v,
                                  psi_minus_v);

  e_float psi_k_minus_v      =  psi_minus_v;             // Psi(-v)
  e_float psi_k_plus_v_plus1 =  psi_v + (ef::one() / v); // Psi(v + 1)
  e_float psi_k_plus1        = -ef::euler_gamma();       // Psi(1)

  // Prepare the initial values of the power terms.
  const e_float one_plus_x_half = (ef::one() + x) / static_cast<INT32>(2);
  const e_float ln_x_plus1_half = (ef::small_arg(x) ? (ef::log1p(x) - ef::ln2())
                                                    :  ef::log(one_plus_x_half));

  e_float one_plus_x_half_pow_k_over_k_fact2 = ef::one();

  // Prepare the initial values of the pochhammer symbols.
  e_float pochham_v_minus = -v;
  e_float pochham_v_plus1 =  v + ef::one();

  e_float v_minus_p       = pochham_v_minus;
  e_float v_plus1_p       = pochham_v_plus1;

  // Prepare the iterations for the digamma (psi) functions.
  e_float k_minus_v      = v_minus_p;
  e_float k_plus_v_plus1 = v_plus1_p;
  e_float k_plus1        = ef::one();

  e_float Pv = ((psi_k_minus_v + psi_k_plus_v_plus1) - (psi_k_plus1 * static_cast<INT32>(2))) + ln_x_plus1_half;

  for(INT32 k = static_cast<INT32>(1); k < ef::max_iteration(); k++)
  {
    // Multiply with [(1 + x) / 2] and divide with the next term of (k!)^2.
    one_plus_x_half_pow_k_over_k_fact2 *= one_plus_x_half;
    one_plus_x_half_pow_k_over_k_fact2 /= k;
    one_plus_x_half_pow_k_over_k_fact2 /= k;

    const e_float term_v = ((pochham_v_minus * pochham_v_plus1) * one_plus_x_half_pow_k_over_k_fact2);

    // Increment and multiply the pochhammer symbols. Note that this is done after the
    // evaluation of term_nu. This is because the series summation begins with 1 and
    // the pochhammer symbols have already been initialized to the values at k = 1.
    pochham_v_minus *= ++v_minus_p;
    pochham_v_plus1 *= ++v_plus1_p;

    // Do the forward recursion of the psi functions using Psi(z + 1) = Psi(z) + 1/z.
    psi_k_minus_v      += (ef::one() / k_minus_v++);
    psi_k_plus_v_plus1 += (ef::one() / k_plus_v_plus1++);
    psi_k_plus1        += (ef::one() / k_plus1++);

    const e_float term_psi = ((psi_k_minus_v + psi_k_plus_v_plus1) - (psi_k_plus1 * static_cast<INT32>(2))) + ln_x_plus1_half;

    const e_float term = term_v * term_psi;

    if(term_v.order() < -ef::tol())
    {
      break;
    }

    Pv += term;
  }

  return (ef::sin(ef::pi() * v) / ef::pi()) * Pv;
}

e_float ef::legendre_p(const e_float& v, const e_float& x)
{
  static const Legendre_Series::LegendrePv L_Pv;

  return L_Pv.MyLegendre(v, x);
}
