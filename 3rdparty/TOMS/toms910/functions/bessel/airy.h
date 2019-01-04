
#ifndef _AIRY_2008_02_16_H_
  #define _AIRY_2008_02_16_H_

  #include <deque>

  #include <e_float/e_float.h>

  namespace ef
  {
    e_float airy_a      (const e_float& x);
    e_float airy_a_prime(const e_float& x);
    e_float airy_b      (const e_float& x);
    e_float airy_b_prime(const e_float& x);

    std::deque<e_float> airy_a_zero(const UINT32 k);
    std::deque<e_float> airy_b_zero(const UINT32 k);
  }
  
  namespace AiryZero
  {
    double ai_estimate_sth_zero(const UINT32 s);
  }

#endif // _AIRY_2008_02_16_H_
