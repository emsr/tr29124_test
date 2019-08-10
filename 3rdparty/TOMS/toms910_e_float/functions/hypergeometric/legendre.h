
#ifndef _LEGENDRE_2007_12_28_H_
  #define _LEGENDRE_2007_12_28_H_

  namespace ef
  {
    e_float legendre_p(const e_float& v, const e_float& x);
    e_float legendre_p(const e_float& v, const e_float& u, const e_float& x);
    e_float legendre_p(const e_float& v, const INT32    m, const e_float& x);
    e_float legendre_p(const INT32    n, const e_float& u, const e_float& x);
    e_float legendre_p(const INT32    n, const INT32    m, const e_float& x);

    e_float legendre_q(const e_float& v, const e_float& x);
    e_float legendre_q(const e_float& v, const e_float& u, const e_float& x);
    e_float legendre_q(const e_float& v, const INT32    m, const e_float& x);
    e_float legendre_q(const INT32    n, const e_float& u, const e_float& x);
    e_float legendre_q(const INT32    n, const INT32    m, const e_float& x);
  }

#endif // _LEGENDRE_2007_12_28_H_
