
#ifndef _ELLIPTIC_2008_03_27_H_
  #define _ELLIPTIC_2008_03_27_H_
  
  #include <e_float/e_float.h>

  namespace ef
  {
    e_float comp_ellint_1(const e_float& m);                     // K(m)      = EllipticK[m]
    e_float      ellint_1(const e_float& m, const e_float& phi); // F(phi, m) = EllipticF[phi, m]
    e_float comp_ellint_2(const e_float& m);                     // E(m)      = EllipticE[m]
    e_float      ellint_2(const e_float& m, const e_float& phi); // E(phi, m) = EllipticE[phi, m]
  }

#endif // _ELLIPTIC_2008_03_27_H_
