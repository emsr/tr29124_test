
#ifndef _BESSEL_2008_02_23_H_
  #define _BESSEL_2008_02_23_H_

  #include <deque>

  #include <e_float/e_float.h>
  #include <functions/functions.h>

  namespace ef
  {
    e_float cyl_bessel_i(const e_float& v, const e_float& x);
    e_float cyl_bessel_i(const INT32 n,    const e_float& x);
    e_float cyl_bessel_j(const e_float& v, const e_float& x);
    e_float cyl_bessel_j(const INT32 n,    const e_float& x);
    e_float cyl_bessel_k(const e_float& v, const e_float& x);
    e_float cyl_bessel_k(const INT32 n,    const e_float& x);
    e_float cyl_bessel_y(const e_float& v, const e_float& x);
    e_float cyl_bessel_y(const INT32 n,    const e_float& x);

    inline e_float cyl_bessel_i_prime(const e_float& v, const e_float& x) { return ((v * cyl_bessel_i(v, x)) / x) - cyl_bessel_i(v + ef::one(), x); }
    inline e_float cyl_bessel_i_prime(const INT32 n,    const e_float& x) { return ((n * cyl_bessel_i(n, x)) / x) - cyl_bessel_i(n + static_cast<INT32>(1), x); }
    inline e_float cyl_bessel_j_prime(const e_float& v, const e_float& x) { return ((v * cyl_bessel_j(v, x)) / x) - cyl_bessel_j(v + ef::one(), x); }
    inline e_float cyl_bessel_j_prime(const INT32 n,    const e_float& x) { return ((n * cyl_bessel_j(n, x)) / x) - cyl_bessel_j(n + static_cast<INT32>(1), x); }
    inline e_float cyl_bessel_k_prime(const e_float& v, const e_float& x) { return ((v * cyl_bessel_k(v, x)) / x) - cyl_bessel_k(v + ef::one(), x); }
    inline e_float cyl_bessel_k_prime(const INT32 n,    const e_float& x) { return ((n * cyl_bessel_k(n, x)) / x) - cyl_bessel_k(n + static_cast<INT32>(1), x); }
    inline e_float cyl_bessel_y_prime(const e_float& v, const e_float& x) { return ((v * cyl_bessel_y(v, x)) / x) - cyl_bessel_y(v + ef::one(), x); }
    inline e_float cyl_bessel_y_prime(const INT32 n,    const e_float& x) { return ((n * cyl_bessel_y(n, x)) / x) - cyl_bessel_y(n + static_cast<INT32>(1), x); }

    std::deque<e_float> cyl_bessel_j_zero(const e_float& v, const UINT32 k);
    std::deque<e_float> cyl_bessel_j_zero(const INT32 n,    const UINT32 k);
  }

#endif // _BESSEL_2008_02_23_H_
