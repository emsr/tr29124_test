
#include <vector>
#include <algorithm>

#include <functions/integer/prime.h>

namespace Primes
{
  struct Inserter
  {
  private:

    mutable UINT32 count;
    mutable std::back_insert_iterator<std::deque<UINT32> > it;

  public:

    static const std::size_t start_index = static_cast<std::size_t>(2u);

    explicit Inserter(std::deque<UINT32>& sequence) : count(static_cast<UINT32>(start_index)),
                                                      it   (std::back_inserter(sequence)) { }

    void operator()(const bool& bo_is_not_prime) const
    {
      const bool bo_is_prime = !bo_is_not_prime;

      if(bo_is_prime)
      {
        *it = count;
      }

      ++count;
    }
  };

  static void Generator(const UINT32 n, std::deque<UINT32>& primes_data)
  {
    // Establish the range of the prime number calculation. Use an approximation
    // related to the prime number theorem to obtain the value of the maximum prime
    // number or a minimum of at least 100. Also be sure to limit this range to
    // within the upper limit of UINT32.

    static const UINT32 min_hundred = static_cast<UINT32>(100u);
    static const double xmax        = static_cast<double>(std::numeric_limits<UINT32>::max());

    const UINT32 N       = std::max(min_hundred, n);
    const double xn      = static_cast<double>(N);
    const double logn    = ::log(xn);
    const double loglogn = ::log(logn);
    const double top     = xn * (((logn + loglogn) - static_cast<double>(1.0)) + ((static_cast<double>(1.8) * loglogn) / logn));
    const double xlim    = std::min(top, xmax);
    const UINT32 nlim    = static_cast<UINT32>(static_cast<UINT64>(xlim));
    const UINT32 limit   = std::max(n, nlim);

    // Use a sieve algorithm to generate a boolean table representation of the primes.

    std::vector<bool> sieve(static_cast<std::size_t>(limit), false);

    UINT32 i = static_cast<UINT32>(Primes::Inserter::start_index);
    UINT32 i2;

    while((i2 = static_cast<UINT32>(i * i)) < limit)
    {
      if(!sieve[i])
      {
        for(UINT32 j = i2; j < limit; j = static_cast<UINT32>(j + i))
        {
          sieve[j] = true;
        }
      }

      ++i;
    }

    // Extract the prime numbers into the data table by inserting them from the sieve.
    primes_data.clear();

    std::for_each(sieve.begin() + Primes::Inserter::start_index,
                  sieve.end(),
                  Primes::Inserter(primes_data));

    primes_data.resize(static_cast<std::size_t>(n), static_cast<UINT32>(0u));
  }

  std::deque<UINT32>& Data(void)
  {
    // Create a static data table of primes and return a reference to it.
    static std::deque<UINT32> primes;

    if(primes.empty())
    {
      // Select a maximum count of prime numbers to be stored in the data table.
      // This number is selected such that the value of the highest prime will slightly
      // exceed 0x10000 (decimal 65,536). This number is significant because it is
      // the maximum value which needs to be tested while computing the prime factors
      // of unsigned 32-bit integers, as done in the subroutine Factors(...).
      Primes::Generator(static_cast<UINT32>(6550u), primes);
    }

    return primes;
  }
}

void ef::prime(const UINT32 n, std::deque<UINT32>& primes)
{
  // For small values of n less than the size of the prime data table, the primes
  // can be copied from the data table. For large values of n, the primes must be
  // generated.
  if(n < static_cast<UINT32>(Primes::Data().size()))
  {
    primes.assign(Primes::Data().begin(), Primes::Data().begin() + static_cast<std::size_t>(n));
  }
  else
  {
    Primes::Generator(n, primes);
  }
}
