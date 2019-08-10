
#include <e_float/e_float.h>
#include <functions/constants/constants.h>
#include <functions/elementary/elementary.h>
#include <functions/integer/integer.h>
#include <functions/gamma/gamma.h>
#include <functions/tables/tables.h>

e_float ef::bernoulli(const UINT32 n)
{
  if(static_cast<UINT32>(n % static_cast<UINT32>(2u)) != static_cast<UINT32>(0u))
  {
    return (n == static_cast<UINT32>(1u) ? -ef::half() : ef::zero());
  }
  else
  {
    static const std::size_t sz_A000367 = Tables::A000367().size();
    static const std::size_t sz_A002445 = Tables::A002445().size();
    static const std::size_t sz_max     = std::min(sz_A000367, sz_A002445);

    const std::size_t n_half = static_cast<std::size_t>(n / static_cast<UINT32>(2u));

    if(n_half < sz_max)
    {
      return Tables::A000367()[n_half]() / Tables::A002445()[n_half]();
    }
    else
    {
      // Do a loop calculation for higher numbered Bernoulli numbers.
      // See Abramowitz & Stegun 23.1.18, page 805, for x = 0.
      // See Computation of Special Functions, Zhang & Jin, 1.1.16, page 5.

      e_float sum = ef::one();

      // TBD: Check the power of two using logs and floating point math to see
      //      if the loop is even necessary.

      for(INT32 k = static_cast<INT32>(2); k < ef::max_iteration(); k++)
      {
        const e_float one_over_k = ef::one() / k;
        const e_float term       = ef::pown(one_over_k, static_cast<INT64>(n));
        
        if(term.order() < -ef::tol())
        {
          break;
        }

        sum += term;
      }

      const bool b_neg = static_cast<UINT32>(static_cast<UINT32>(n / static_cast<UINT32>(2u)) & static_cast<UINT32>(1u)) == static_cast<UINT32>(0u);

      const e_float factor = ((ef::factorial(n) / ef::pown(ef::two_pi(), static_cast<INT64>(n))) * static_cast<INT32>(2));
      const e_float bn     = sum * factor;
      
      return !b_neg ? bn : -bn;
    }
  }
}
