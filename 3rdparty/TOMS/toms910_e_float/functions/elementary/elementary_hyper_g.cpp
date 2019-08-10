
#include <vector>
#include <algorithm>
#include <numeric>
#include <functional>

#include <e_float/e_float.h>
#include <functions/constants/constants.h>
#include <functions/elementary/elementary.h>

e_float ef::hyp0F0(const e_float& x)
{
  // Compute the series representation of Hypergeometric0F0 taken from
  // http://functions.wolfram.com/HypergeometricFunctions/Hypergeometric0F0/06/01/
  // There are no checks on input range or parameter boundaries.

  e_float x_pow_n_div_n_fact(x);
  
  e_float H0F0 = ef::one() + x_pow_n_div_n_fact;

  INT32 n;

  // Series expansion of hyperg_0f0(; ; x).
  for(n = static_cast<INT32>(2); n < ef::max_iteration(); n++)
  {
    x_pow_n_div_n_fact *= x;
    x_pow_n_div_n_fact /= n;

    const INT64 order_check = static_cast<INT64>(x_pow_n_div_n_fact.order() - H0F0.order());

    if((n > static_cast<INT32>(20)) && (order_check < -ef::tol()))
    {
      break;
    }

    H0F0 += x_pow_n_div_n_fact;
  }

  return n < ef::max_iteration() ? H0F0 : std::numeric_limits<e_float>::quiet_NaN();
}

e_float ef::hyp0F1(const e_float& b, const e_float& x)
{
  // Compute the series representation of Hypergeometric0F1 taken from
  // http://functions.wolfram.com/HypergeometricFunctions/Hypergeometric0F1/06/01/01/
  // There are no checks on input range or parameter boundaries.

  e_float x_pow_n_div_n_fact(x);
  e_float pochham_b         (b);
  e_float bp                (b);
  
  e_float H0F1 = ef::one() + (x_pow_n_div_n_fact / pochham_b);

  INT32 n;

  // Series expansion of hyperg_0f1(; b; x).
  for(n = static_cast<INT32>(2); n < ef::max_iteration(); n++)
  {
    x_pow_n_div_n_fact *= x;
    x_pow_n_div_n_fact /= n;

    pochham_b *= ++bp;

    const e_float term = x_pow_n_div_n_fact / pochham_b;

    const INT64 order_check = static_cast<INT64>(term.order() - H0F1.order());

    if((n > static_cast<INT32>(20)) && (order_check < -ef::tol()))
    {
      break;
    }

    H0F1 += term;
  }

  return n < ef::max_iteration() ? H0F1 : std::numeric_limits<e_float>::quiet_NaN();
}

e_float ef::hyp1F0(const e_float& a, const e_float& x)
{
  // Compute the series representation of Hypergeometric1F0 taken from
  // http://functions.wolfram.com/HypergeometricFunctions/Hypergeometric1F0/06/01/01/
  // and also see the corresponding section for the power function (i.e. x^a).
  // There are no checks on input range or parameter boundaries.

  e_float x_pow_n_div_n_fact(x);
  e_float pochham_a         (a);
  e_float ap                (a);

  e_float H1F0 = ef::one() + (pochham_a * x_pow_n_div_n_fact);

  INT32 n;

  // Series expansion of hyperg_1f0(a; ; x).
  for(n = static_cast<INT32>(2); n < ef::max_iteration(); n++)
  {
    x_pow_n_div_n_fact *= x;
    x_pow_n_div_n_fact /= n;

    pochham_a *= ++ap;

    const e_float term = pochham_a * x_pow_n_div_n_fact;

    const INT64 order_check = static_cast<INT64>(term.order() - H1F0.order());

    if((n > static_cast<INT32>(20)) && (order_check < -ef::tol()))
    {
      break;
    }

    H1F0 += term;
  }

  return n < ef::max_iteration() ? H1F0 : std::numeric_limits<e_float>::quiet_NaN();
}

e_float ef::hyp1F1(const e_float& a, const e_float& b, const e_float& x)
{
  // Compute the series representation of hyperg_1f1 taken from
  // Abramowitz and Stegun 13.1.2, page 504.
  // There are no checks on input range or parameter boundaries.

  e_float x_pow_n_div_n_fact(x);
  e_float pochham_a         (a);
  e_float pochham_b         (b);
  e_float ap                (a);
  e_float bp                (b);
  
  e_float H1F1 = ef::one() + ((pochham_a / pochham_b) * x_pow_n_div_n_fact);

  INT32 n;

  // Series expansion of hyperg_1f1(a, b ; x).
  for(n = static_cast<INT32>(2); n < ef::max_iteration(); n++)
  {
    x_pow_n_div_n_fact *= x;
    x_pow_n_div_n_fact /= n;

    pochham_a *= ++ap;
    pochham_b *= ++bp;

    const e_float term = (pochham_a / pochham_b) * x_pow_n_div_n_fact;

    const INT64 order_check = static_cast<INT64>(term.order() - H1F1.order());

    if((n > static_cast<INT32>(20)) && (order_check < -ef::tol()))
    {
      break;
    }

    H1F1 += term;
  }

  return n < ef::max_iteration() ? H1F1 : std::numeric_limits<e_float>::quiet_NaN();
}

e_float ef::hyp2F0(const e_float& a, const e_float& b, const e_float& x)
{
  // Compute the series representation of hyperg_2f0.
  // There are no checks on input range or parameter boundaries.

  e_float x_pow_n_div_n_fact(x);
  e_float pochham_a         (a);
  e_float pochham_b         (b);
  e_float ap                (a);
  e_float bp                (b);

  e_float H2F0 = ef::one() + ((pochham_a * pochham_b) * x_pow_n_div_n_fact);

  INT32 n;

  // Series expansion of hyperg_2f0(a, b; ; x).
  for(n = static_cast<INT32>(2); n < ef::max_iteration(); n++)
  {
    x_pow_n_div_n_fact *= x;
    x_pow_n_div_n_fact /= n;

    pochham_a *= ++ap;
    pochham_b *= ++bp;

    const e_float term = (pochham_a * pochham_b) * x_pow_n_div_n_fact;

    const INT64 order_check = static_cast<INT64>(term.order() - H2F0.order());

    if((n > static_cast<INT32>(20)) && (order_check < -ef::tol()))
    {
      break;
    }

    H2F0 += term;
  }

  return n < ef::max_iteration() ? H2F0 : std::numeric_limits<e_float>::quiet_NaN();
}

e_float ef::hyp2F1(const e_float& a, const e_float& b, const e_float& c, const e_float& x)
{
  // Compute the series representation of hyperg_2f1 taken from
  // Abramowitz and Stegun 15.1.1.
  // There are no checks on input range or parameter boundaries.

  e_float x_pow_n_div_n_fact(x);
  e_float pochham_a         (a);
  e_float pochham_b         (b);
  e_float pochham_c         (c);
  e_float ap                (a);
  e_float bp                (b);
  e_float cp                (c);

  e_float H2F1 = ef::one() + (((pochham_a * pochham_b) / pochham_c) * x_pow_n_div_n_fact);

  INT32 n;

  // Series expansion of hyperg_2f1(a, b; c; x).
  for(n = static_cast<INT32>(2); n < ef::max_iteration(); n++)
  {
    x_pow_n_div_n_fact *= x;
    x_pow_n_div_n_fact /= n;

    pochham_a *= ++ap;
    pochham_b *= ++bp;
    pochham_c *= ++cp;

    const e_float term = ((pochham_a * pochham_b) / pochham_c) * x_pow_n_div_n_fact;

    const INT64 order_check = static_cast<INT64>(term.order() - H2F1.order());

    if((n > static_cast<INT32>(20)) && (order_check < -ef::tol()))
    {
      break;
    }

    H2F1 += term;
  }

  return n < ef::max_iteration() ? H2F1 : std::numeric_limits<e_float>::quiet_NaN();
}

e_float ef::hypPFQ(const std::deque<e_float>& a, const  std::deque<e_float>& b, const e_float& x)
{
  // Compute the series representation of hyperg_pfq.
  // There are no checks on input range or parameter boundaries.

  e_float x_pow_n_div_n_fact(x);

  // The pochhammer symbols for the multiplications in the series expansion
  // will be stored in STL-containers.
  std::vector<e_float> ap(a.begin(), a.end());
  std::vector<e_float> bp(b.begin(), b.end());

  // Initialize the pochhammer product terms with the products of the form:
  // [(a0)_1 * (a1)_1 * (a2)_1 * ...], or [(b0)_1 * (b1)_1 * (b2)_1 * ...].
  e_float pochham_a = std::accumulate(ap.begin(), ap.end(), ef::one(), std::multiplies<e_float>());
  e_float pochham_b = std::accumulate(bp.begin(), bp.end(), ef::one(), std::multiplies<e_float>());

  e_float HPFQ = ef::one() + ((pochham_a / pochham_b) * x_pow_n_div_n_fact);

  INT32 n;

  // Series expansion of hyperg_pfq[{a0, a1, a2, ...}; {b0, b1, b2, ...}; x].
  for(n = static_cast<INT32>(2); n < ef::max_iteration(); n++)
  {
    x_pow_n_div_n_fact *= x;
    x_pow_n_div_n_fact /= n;

    // Increment each of the pochhammer elements in a and b.
    std::transform(ap.begin(), ap.end(), ap.begin(), std::bind1st(std::plus<e_float>(), ef::one()));
    std::transform(bp.begin(), bp.end(), bp.begin(), std::bind1st(std::plus<e_float>(), ef::one()));

    // Multiply the pochhammer product terms with the products of the incremented
    // pochhammer elements. These are products of the form:
    // [(a0)_k * (a1)_k * (a2)_k * ...], or [(b0)_k * (b1)_k * (b2)_k * ...].
    pochham_a *= std::accumulate(ap.begin(), ap.end(), ef::one(), std::multiplies<e_float>());
    pochham_b *= std::accumulate(bp.begin(), bp.end(), ef::one(), std::multiplies<e_float>());

    const e_float term = (pochham_a / pochham_b) * x_pow_n_div_n_fact;

    const INT64 order_check = static_cast<INT64>(term.order() - HPFQ.order());

    if((n > static_cast<INT32>(20)) && (order_check < -ef::tol()))
    {
      break;
    }

    HPFQ += term;
  }

  return n < ef::max_iteration() ? HPFQ : std::numeric_limits<e_float>::quiet_NaN();
}
