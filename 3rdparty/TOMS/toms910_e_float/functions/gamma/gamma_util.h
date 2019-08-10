
#ifndef _GAMMA_UTIL_2008_01_10_H_
  #define _GAMMA_UTIL_2008_01_10_H_
  
  #include <e_float/e_float.h>

  namespace GammaUtil
  {
    void   GammaOfPlusXMinusX(const e_float& x, e_float& gamma_plus_x, e_float& gamma_minus_x);
    void DiGammaOfPlusXMinusX(const e_float& x, e_float& psi_plus_x,   e_float& psi_minus_x);
  }
  
#endif // _GAMMA_UTIL_2008_01_10_H_
