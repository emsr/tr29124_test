
#include <functions/integer/integer.h>
#include <functions/tables/tables.h>

e_float ef::stirling2(const UINT32 n, const UINT32 k)
{
  if(n == static_cast<UINT32>(0u))
  {
    return k == static_cast<UINT32>(0u) ? ef::one() : ef::zero();
  }

  if(k == static_cast<UINT32>(0u) || (k > n))
  {
    return ef::zero();
  }

  if((k == n) || (k == static_cast<UINT32>(1u)))
  {
    return ef::one();
  }

  const std::size_t n_minus_one = static_cast<std::size_t>(n - static_cast<UINT32>(1u));
  const std::size_t k_minus_one = static_cast<std::size_t>(k - static_cast<UINT32>(1u));

  if(n_minus_one < static_cast<UINT32>(Tables::A008277().size()))
  {
    return Tables::A008277()[n_minus_one]()[k_minus_one];
  }
  else
  {
    // Use recursion to calculate higher Stirling numbers of the second kind.
    // In order to limit the code complexity, the complete triangle is computed even
    // if this is not needed for the actual column in question. The recursion results
    // are not stored and subsequent calculations must carry out the full recursion.

    static const std::size_t table_size_plus_one = static_cast<std::size_t>(Tables::A008277().size() + static_cast<std::size_t>(1u));

    std::vector<e_float> s2_row_previous = (Tables::A008277().back())();
    std::vector<e_float> s2_row;

    for(std::size_t nn = table_size_plus_one; nn <= static_cast<std::size_t>(n); nn++)
    {
      s2_row.resize(nn);

      s2_row[static_cast<std::size_t>(0u)] = ef::one();

      for(std::size_t kk = static_cast<std::size_t>(1u); kk < static_cast<std::size_t>(nn - static_cast<std::size_t>(1u)); kk++)
      {
        s2_row[kk] = s2_row_previous[kk - 1u] + (s2_row_previous[kk] * static_cast<INT32>(kk + static_cast<std::size_t>(1u)));
      }

      s2_row.back() = ef::one();

      s2_row_previous = s2_row;
    }

    return s2_row[k_minus_one];
  }
}
