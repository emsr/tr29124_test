
#ifndef _LAGUERRE_2009_02_04_H_
  #define _LAGUERRE_2009_02_04_H_

  namespace ef
  {
    e_float laguerre(const e_float& v, const e_float& x);
    e_float laguerre(const e_float& v, const e_float& L, const e_float& x);
    e_float laguerre(const INT32 n,    const e_float& L, const e_float& x);
    e_float laguerre(const INT32 n,    const INT32 m,    const e_float& x);
  }

#endif // _LAGUERRE_2009_02_04_H_
