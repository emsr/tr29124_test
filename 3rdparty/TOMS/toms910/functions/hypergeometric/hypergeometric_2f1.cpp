
#include <functions/constants/constants.h>
#include <functions/elementary/elementary.h>
#include <functions/gamma/gamma.h>
#include <functions/gamma/gamma_util.h>
#include <functions/hypergeometric/hypergeometric.h>

namespace Hypergeometric2F1_Series
{
  static e_float AtNegativeIntegerA(const INT32 an, const e_float& b, const e_float& c, const e_float& x);
}

static e_float Hypergeometric2F1_Series::AtNegativeIntegerA(const INT32 an, const e_float& b, const e_float& c, const e_float& x)
{
  // Implement Abramowitz & Stegun 15.4.1 and 15.4.2.
  e_float x_pow_n_div_n_fact(x);
  e_float pochham_m         (an);
  e_float pochham_b         (b);
  e_float pochham_c         (c);
  e_float mp                (an);
  e_float bp                (b);
  e_float cp                (c);

  e_float H2F1 = ef::one() + (((pochham_m * pochham_b) / pochham_c) * x_pow_n_div_n_fact);

  // Series expansion of hyperg_1f1(a, b ; x).
  for(INT32 n = static_cast<INT32>(2); n <= static_cast<INT32>(-an); n++)
  {
    x_pow_n_div_n_fact *= x;
    x_pow_n_div_n_fact /= n;
    
    pochham_m *= ++mp;
    pochham_b *= ++bp;
    pochham_c *= ++cp;

    const e_float term = ((pochham_m * pochham_b) / pochham_c) * x_pow_n_div_n_fact;

    const INT64 order_check = static_cast<INT64>(term.order() - H2F1.order());

    if((n > static_cast<INT32>(20)) && (order_check < -ef::tol()))
    {
      break;
    }

    H2F1 += term;
  }
  
  return H2F1;
}

e_float ef::hyperg_2f1(const e_float& a, const e_float& b, const e_float& c, const e_float& x)
{
  // See the description near Abramowitz and Stegun 15.1.1, page 556.

  if(ef::iszero(c))
  {
    return std::numeric_limits<e_float>::quiet_NaN();
  }
  
  if(ef::isone(x))
  {
    return ef::zero();
  }

  if(ef::isone(-x))
  {
    return ef::zero();
  }

  // Check for integer parameters and take the necessary actions for such parameters.

  const bool a_is_negative_integer = (ef::isint(a) && ef::isneg(a));
  const bool b_is_negative_integer = (ef::isint(b) && ef::isneg(b));
  const bool c_is_negative_integer = (ef::isint(c) && ef::isneg(c));

  const INT32 an = ef::to_int32(a);
  const INT32 bn = ef::to_int32(b);

  if(c_is_negative_integer)
  {
    const INT32 cn = static_cast<INT32>(ef::to_int64(c));

    if(!a_is_negative_integer && !b_is_negative_integer)
    {
      return std::numeric_limits<e_float>::quiet_NaN();
    }

    if(a_is_negative_integer && b_is_negative_integer)
    {
      if(cn > std::max(an, bn))
      {
        return std::numeric_limits<e_float>::quiet_NaN();
      }
    }

    if(a_is_negative_integer)
    {
      return cn > an ? std::numeric_limits<e_float>::quiet_NaN()
                     : Hypergeometric2F1_Series::AtNegativeIntegerA(an, b, c, x);
    }

    if(b_is_negative_integer)
    {
      return cn > bn ? std::numeric_limits<e_float>::quiet_NaN()
                     : Hypergeometric2F1_Series::AtNegativeIntegerA(bn, a, c, x);
    }
  }

  if(a_is_negative_integer && b_is_negative_integer)
  {
    return an > bn ? Hypergeometric2F1_Series::AtNegativeIntegerA(an, b, c, x)
                   : Hypergeometric2F1_Series::AtNegativeIntegerA(bn, a, c, x);
  }

  if(a_is_negative_integer)
  {
    return Hypergeometric2F1_Series::AtNegativeIntegerA(an, b, c, x);
  }

  if(b_is_negative_integer)
  {
    return Hypergeometric2F1_Series::AtNegativeIntegerA(bn, a, c, x);
  }

  return ef::hyp2F1(a, b, c, x);
}

e_float ef::hyperg_2f1_reg(const e_float& a, const e_float& b, const e_float& c, const e_float& x)
{
  // As described in The Mathematica Book, Fourth Edition, page 1158.
  // The series converges for all x, |x| < 1, real and complex.
  // Also include the limit c --> -m.

  if(ef::isint(c))
  {
    const INT32 m = ef::to_int32(c);

    if(m <= static_cast<INT32>(0))
    {
      // Implement Abramowitz & Stegun 15.1.2, page 556.
      // See also:
      // http://functions.wolfram.com/HypergeometricFunctions/Hypergeometric2F1Regularized/17/02/04/
      const INT32   m_plus1       = static_cast<INT32>(-m) + static_cast<INT32>(1u);
      const e_float factor        =   (ef::pochhammer(a, m_plus1) * ef::pochhammer(b, m_plus1))
                                    *  ef::pown(x, static_cast<INT64>(m_plus1));
      const e_float ef_m_plus1    = e_float(m_plus1);

      return factor * ef::hyperg_2f1_reg(a + ef_m_plus1,
                                         b + ef_m_plus1,
                                         e_float(static_cast<INT32>(-m) + static_cast<UINT32>(2u)),
                                         x);
    }
  }

  return ef::hyperg_2f1(a, b, c, x) / ef::gamma(c);
}
