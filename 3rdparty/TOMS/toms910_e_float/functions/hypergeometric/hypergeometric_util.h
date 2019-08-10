
#ifndef _HYPERGEOMETRIC_UTIL_2008_03_01_H_
  #define _HYPERGEOMETRIC_UTIL_2008_03_01_H_

  #include <deque>

  #include <e_float/e_float.h>

  namespace HypergeometricUtil
  {
    const INT32& AsympConvergeMaxOrder(void);

    bool AsympConverge(const std::deque<double>& a,
                       const std::deque<double>& b,
                       const double x,
                       const INT32  prec = static_cast<INT32>(ef::tol()),
                       const INT32  maxo = AsympConvergeMaxOrder());
  }

#endif // _HYPERGEOMETRIC_UTIL_2008_03_01_H_
