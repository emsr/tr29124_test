
#ifndef _BESSEL_RECURSION_2008_02_21_H_
  #define _BESSEL_RECURSION_2008_02_21_H_

  #include <vector>

  #include <functions/complex/e_float_complex.h>

  namespace BesselRecursion
  {
    e_float RecurJn(const INT32 n,
                    const e_float& x,
                    std::deque<e_float>* const pJn = static_cast<std::deque<e_float>* const>(0u));
    e_float RecurJv(const e_float& v,
                    const e_float& x,
                    std::deque<e_float>* const pJv = static_cast<std::deque<e_float>* const>(0u));
    e_float RecurIn(const INT32 n,
                    const e_float& x,
                    std::deque<e_float>* const pIn = static_cast<std::deque<e_float>* const>(0u));
    e_float RecurIv(const e_float& v,
                    const e_float& x,
                    std::deque<e_float>* const pIv = static_cast<std::deque<e_float>* const>(0u));

    ef_complex RecurJv(const e_float& v,
                       const ef_complex& z,
                       std::deque<ef_complex>* const pJv = static_cast<std::deque<ef_complex>* const>(0u));
  }
  
#endif // _BESSEL_RECURSION_2008_02_21_H_
