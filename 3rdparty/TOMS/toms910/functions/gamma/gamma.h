
#ifndef _GAMMA_2008_01_01_H_
  #define _GAMMA_2008_01_01_H_

  #include <functions/complex/e_float_complex.h>

  namespace ef
  {
    e_float gamma               (const e_float& x);
    e_float gamma_near_n        (const INT32 n, const e_float& x);
    e_float incomplete_gamma    (const e_float& a, const e_float& x);
    e_float gen_incomplete_gamma(const e_float& a, const e_float& x0, const e_float& x1);
    e_float beta                (const e_float& a, const e_float& b);
    e_float incomplete_beta     (const e_float& x, const e_float& a, const e_float& b);
    e_float factorial           (const UINT32 n);
    e_float factorial2          (const  INT32 n);
    e_float binomial            (const UINT32 n, const UINT32 k);
    e_float binomial            (const UINT32 n, const e_float& y);
    e_float binomial            (const e_float& x, const UINT32 k);
    e_float binomial            (const e_float& x, const e_float& y);
    e_float pochhammer          (const e_float& x, const UINT32 n);
    e_float pochhammer          (const e_float& x, const e_float& a);
  }

  namespace efz
  {
    ef_complex gamma     (const ef_complex& z);
    ef_complex beta      (const ef_complex& a, const ef_complex& b);
    ef_complex pochhammer(const ef_complex& z, const UINT32 n);
    ef_complex pochhammer(const ef_complex& z, const ef_complex& a);
  }

#endif // _GAMMA_2008_01_01_H_
