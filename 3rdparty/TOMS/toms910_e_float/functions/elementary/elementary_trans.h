
#ifndef _ELEMENTARY_TRANS_2009_12_09_H_
  #define _ELEMENTARY_TRANS_2009_12_09_H_

  namespace ef
  {
    e_float pow2    (const INT64 p);
    e_float pown    (const e_float& x, const INT64 p);
    e_float inv     (const e_float& x);
    e_float sqrt    (const e_float& x);
    e_float cbrt    (const e_float& x);
    e_float rootn   (const e_float& x, const INT32 p);
    e_float exp     (const e_float& x);
    e_float log     (const e_float& x);
    e_float log10   (const e_float& x);
    e_float loga    (const e_float& a, const e_float& x);
    e_float log1p   (const e_float& x);
    e_float log1p1m2(const e_float& x);
    e_float pow     (const e_float& x, const e_float& a);
    void    sinhcosh(const e_float& x, e_float* const p_sin, e_float* const p_cos);
    e_float sinh    (const e_float& x);
    e_float cosh    (const e_float& x);
    e_float tanh    (const e_float& x);
    e_float asinh   (const e_float& x);
    e_float acosh   (const e_float& x);
    e_float atanh   (const e_float& x);
  }

#endif // _ELEMENTARY_TRANS_2009_12_09_H_
