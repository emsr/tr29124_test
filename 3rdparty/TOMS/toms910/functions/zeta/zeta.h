
#ifndef _ZETA_2008_09_12_H_
  #define _ZETA_2008_09_12_H_

  #include <functions/complex/e_float_complex.h>

  namespace ef
  {
    e_float riemann_zeta(const INT32 n);
    e_float riemann_zeta(const e_float& s);
    e_float hurwitz_zeta(const e_float& s, const INT32 n);
    e_float hurwitz_zeta(const e_float& s, const e_float& a);
  }

  namespace efz
  {
    ef_complex riemann_zeta(const ef_complex& s);
    ef_complex hurwitz_zeta(const ef_complex& s, const INT32 n);
    ef_complex hurwitz_zeta(const ef_complex& s, const ef_complex& a);
  }

#endif // _ZETA_2008_09_12_H_
