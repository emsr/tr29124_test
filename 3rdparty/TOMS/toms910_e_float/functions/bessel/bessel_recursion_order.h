
#ifndef _BESSEL_RECURSION_ORDER_2009_11_01_H_
  #define _BESSEL_RECURSION_ORDER_2009_11_01_H_

  #include <e_float/e_float.h>

  namespace BesselRecursionOrder
  {
    INT32 RecursionStartOrderJ0            (const double x);
    INT32 RecursionStartOrderJn            (const INT32 n,  const double x);
    INT32 RecursionStartOrderI0            (const double x);
    INT32 RecursionStartOrderIn            (const INT32 n,  const double x);
    INT32 RecursionStartOrderKv            (const double v, const double x);
    INT32 RecursionStartOrderKv_v_near_half(const double v, const double x);
  }

#endif // _BESSEL_RECURSION_ORDER_2009_11_01_H_
