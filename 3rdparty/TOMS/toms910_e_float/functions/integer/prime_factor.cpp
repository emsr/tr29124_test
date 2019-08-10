
#include <vector>
#include <algorithm>

#include <functions/integer/prime.h>

namespace Primes
{
  std::deque<UINT32>& Data(void);

  static bool IsPrimeFactor(UINT32& np, const UINT32 p)
  {
    const UINT32 q = static_cast<UINT32>(np / p);
    const UINT32 r = static_cast<UINT32>(np - static_cast<UINT32>(q * p));

    const bool is_prime_factor = (r == static_cast<UINT32>(0u));

    if(is_prime_factor)
    {
      np = q;
    }
    
    return is_prime_factor;
  }

  static void Factors(const UINT32 n, std::deque<Util::point<UINT32> >& pf)
  {
    // Compute the prime factors of the unsigned integer n. Use the divide algorithm of
    // "The Art of Computer Programming Volume 2 Semi-numerical Algorithms Third Edition",
    // Donald Knuth (Algorithm A, Chapter 4.5.4, page 380 and pages 378-417).
    static const std::size_t sz = Data().size();

    pf.clear();

    const UINT32 sqrt_n = static_cast<UINT32>(static_cast<UINT64>(::sqrt(static_cast<double>(n)) + static_cast<double>(0.5)));

    UINT32 np = n;

    for(std::size_t i = static_cast<std::size_t>(0u); i < sz; i++)
    {
      const UINT32 p = Data()[i];

      if(IsPrimeFactor(np, p))
      {
        Util::point<UINT32> ip(p, static_cast<UINT32>(1u));

        while(IsPrimeFactor(np, p))
        {
          ++ip.y;
        }

        pf.push_back(ip);
      }

      if(static_cast<UINT32>(np / p) <= p)
      {
        pf.push_back(Util::point<UINT32>(np, static_cast<UINT32>(1u)));

        break;
      }

      if((np == static_cast<UINT32>(1u)) || (p >= sqrt_n))
      {
        break;
      }
    }
  }
}

void ef::prime_factors(const UINT32 n, std::deque<Util::point<UINT32> >& pf)
{
  // Factor the input integer into a list of primes. For small inputs less than 10,000
  // use the tabulated prime factors list. Calculate the prime factors for larger inputs
  // above 10,000.
  static std::vector<std::deque<Util::point<UINT32> > > prime_factors_list;

  if(prime_factors_list.empty())
  {
    // Generate a table of the sets of the first 10,000 integer prime factorizations.
    prime_factors_list.resize(static_cast<std::size_t>(10000u));

    prime_factors_list[static_cast<std::size_t>(0u)] = std::deque<Util::point<UINT32> >(static_cast<std::size_t>(1u), Util::point<UINT32>(static_cast<UINT32>(0u), static_cast<UINT32>(1u)));
    prime_factors_list[static_cast<std::size_t>(1u)] = std::deque<Util::point<UINT32> >(static_cast<std::size_t>(1u), Util::point<UINT32>(static_cast<UINT32>(1u), static_cast<UINT32>(1u)));
    prime_factors_list[static_cast<std::size_t>(2u)] = std::deque<Util::point<UINT32> >(static_cast<std::size_t>(1u), Util::point<UINT32>(static_cast<UINT32>(2u), static_cast<UINT32>(1u)));
    prime_factors_list[static_cast<std::size_t>(3u)] = std::deque<Util::point<UINT32> >(static_cast<std::size_t>(1u), Util::point<UINT32>(static_cast<UINT32>(3u), static_cast<UINT32>(1u)));

    static const UINT32 n_five = static_cast<UINT32>(5u);

    std::deque<UINT32>::const_iterator it_next_prime = std::find(Primes::Data().begin(), Primes::Data().end(), n_five);

    for(std::size_t i = static_cast<std::size_t>(4u); i < prime_factors_list.size(); i++)
    {
      if((it_next_prime != Primes::Data().end()) && (static_cast<UINT32>(i) == *it_next_prime))
      {
        ++it_next_prime;

        prime_factors_list[i] = std::deque<Util::point<UINT32> >(static_cast<std::size_t>(1u),
                                                                 Util::point<UINT32>(static_cast<UINT32>(i),
                                                                 static_cast<UINT32>(1u)));
      }
      else
      {
        Primes::Factors(static_cast<UINT32>(i), prime_factors_list[i]);
      }
    }
  }

  if(static_cast<std::size_t>(n) < prime_factors_list.size())
  {
    pf = prime_factors_list[static_cast<std::size_t>(n)];
  }
  else
  {
    Primes::Factors(n, pf);
  }
}
