
#include <e_float/e_float.h>
#include <functions/constants/constants.h>
#include <functions/elementary/elementary.h>
#include <functions/integer/integer.h>
#include <functions/gamma/gamma.h>
#include <functions/tables/tables.h>

e_float ef::euler(const UINT32 n)
{
  if(static_cast<UINT32>(n % static_cast<UINT32>(2u)) != static_cast<UINT32>(0u))
  {
    return ef::zero();
  }

  const std::size_t n_half = static_cast<std::size_t>(n / static_cast<UINT32>(2u));

  if(n_half < Tables::A000364().size())
  {
    return Tables::A000364()[n_half]();
  }

  // Do a loop calculation for higher numbered Euler numbers.
  // See Abramowitz & Stegun 23.1.18, page 805, for x = 0.
  // See Computation of Special Functions, Zhang & Jin, 1.2.13, page 9.

  e_float sum = ef::one();

  bool b_neg_term = true;

  UINT32 k;
  for(k = static_cast<UINT32>(3u); k < static_cast<UINT32>(ef::max_iteration()); k += static_cast<UINT32>(2u))
  {
    const e_float term = ef::pown(ef::one() / static_cast<INT32>(k), static_cast<INT64>(n + static_cast<UINT32>(1u)));
    
    const INT64 order_check = static_cast<INT64>(term.order() - sum.order());

    if((k > static_cast<UINT32>(20)) && (order_check < -ef::tol()))
    {
      break;
    }

    !b_neg_term ? sum += term : sum -= term;

    b_neg_term = !b_neg_term;
  }
  
  if(k >= static_cast<UINT32>(ef::max_iteration()))
  {
    return std::numeric_limits<e_float>::quiet_NaN();
  }

  static const e_float two_over_pi = ef::two() / ef::pi();

  const e_float en = sum * ((ef::pown(two_over_pi, static_cast<INT64>(n + static_cast<UINT32>(1u))) * ef::factorial(n)) * static_cast<UINT32>(2u));

  const bool b_neg = static_cast<UINT32>(static_cast<UINT32>(n % static_cast<UINT32>(4u)) / static_cast<UINT32>(2u)) != static_cast<UINT32>(0u);

  return !b_neg ? en : -en;
}
