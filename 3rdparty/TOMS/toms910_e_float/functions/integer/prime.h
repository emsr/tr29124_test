
#ifndef _PRIME_2008_09_11_H_
  #define _PRIME_2008_09_11_H_

  #include <deque>

  #include <e_float/e_float.h>
  #include <utility/util_point.h>

  namespace ef
  {
    void prime        (const UINT32 n, std::deque<UINT32>& primes);
    void prime_factors(const UINT32 n, std::deque<Util::point<UINT32> >& pf);
  }

#endif // _PRIME_2008_09_11_H_
