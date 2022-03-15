#ifndef SF_PRIME_H
#define SF_PRIME_H 1

#include <emsr/detail/sf_prime.tcc>

namespace emsr
{

  constexpr uint32_t
  num_primes()
  { return emsr::detail::s_num_primes; }

  constexpr uint32_t
  prime(uint16_t i)
  {
    if (i < emsr::detail::s_num_primes)
      return static_cast<uint32_t>(emsr::detail::s_prime[i])
	   + (i < emsr::detail::s_shift_index ? 0 : emsr::detail::s_shift_value);
    else
      return 0;
  }

} // namespace emsr

#endif // SF_PRIME_H
