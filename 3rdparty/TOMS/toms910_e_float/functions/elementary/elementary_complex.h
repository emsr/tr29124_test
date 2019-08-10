
#ifndef _ELEMENTARY_COMPLEX_2009_12_09_H_
  #define _ELEMENTARY_COMPLEX_2009_12_09_H_

  #include <complex>

  #include <functions/complex/e_float_complex.h>

  namespace efz
  {
    std::complex<double> to_double(const ef_complex& z);

    inline e_float norm(const ef_complex& z) { return z.norm(); }
           e_float abs (const ef_complex& z);
           e_float arg (const ef_complex& z);
    inline e_float real(const ef_complex& z) { return z.real(); }
    inline e_float imag(const ef_complex& z) { return z.imag(); }

    inline ef_complex conj(const ef_complex& z) { return ef_complex(z.real(), -z.imag()); }
    inline ef_complex iz  (const ef_complex& z) { const e_float tmp(z.real()); return ef_complex(-z.imag(), tmp); }

    ef_complex polar   (const e_float& mod, const e_float& arg);
    ef_complex sin     (const ef_complex& z);
    ef_complex cos     (const ef_complex& z);
    ef_complex tan     (const ef_complex& z);
    void       sincos  (const ef_complex& z, ef_complex* const p_sin, ef_complex* const p_cos);
    ef_complex csc     (const ef_complex& z);
    ef_complex sec     (const ef_complex& z);
    ef_complex cot     (const ef_complex& z);
    ef_complex asin    (const ef_complex& z);
    ef_complex acos    (const ef_complex& z);
    ef_complex atan    (const ef_complex& z);
    ef_complex inv     (const ef_complex& z);
    ef_complex sqrt    (const ef_complex& z);
    ef_complex exp     (const ef_complex& z);
    ef_complex log     (const ef_complex& z);
    ef_complex log10   (const ef_complex& z);
    ef_complex loga    (const ef_complex& a, const ef_complex& z);
    ef_complex pown    (const ef_complex& z, const INT64 p);
    ef_complex pow     (const ef_complex& z, const ef_complex& a);
    ef_complex rootn   (const ef_complex& z, const INT32 p);
    ef_complex sinh    (const ef_complex& z);
    ef_complex cosh    (const ef_complex& z);
    ef_complex tanh    (const ef_complex& z);
    void       sinhcosh(const ef_complex& z, ef_complex* const p_sinh, ef_complex* const p_cosh);
    ef_complex asinh   (const ef_complex& z);
    ef_complex acosh   (const ef_complex& z);
    ef_complex atanh   (const ef_complex& z);
  }

  // Overwrite global template ef_complex operator functions, some of which are overloaded.
  // These operator functions are:
  //   operator!=
  //   operator==
  //   operator*
  //   operator+
  //   operator-
  //   operator/
  //   operator<<
  //   TBD: operator>>
  inline bool operator==(const ef_complex& u, const ef_complex& v) { return (u.real() == v.real()) && (u.imag() == v.imag()); }
  inline bool operator!=(const ef_complex& u, const ef_complex& v) { return (u.real() != v.real()) || (u.imag() != v.imag()); }

  bool operator==(const ef_complex& u, const e_float& v);
  bool operator!=(const ef_complex& u, const e_float& v);

  bool operator==(const e_float& u, const ef_complex& v);
  bool operator!=(const e_float& u, const ef_complex& v);

  bool operator==(const ef_complex& u, const INT32 n);
  bool operator!=(const ef_complex& u, const INT32 n);

  bool operator==(const INT32 n, const ef_complex& u);
  bool operator!=(const INT32 n, const ef_complex& u);

  inline ef_complex operator+ (const ef_complex& u, const ef_complex& v) { return ef_complex(u) += v; }
  inline ef_complex operator- (const ef_complex& u, const ef_complex& v) { return ef_complex(u) -= v; }
  inline ef_complex operator* (const ef_complex& u, const ef_complex& v) { return ef_complex(u) *= v; }
  inline ef_complex operator/ (const ef_complex& u, const ef_complex& v) { return ef_complex(u) /= v; }

  inline ef_complex operator+ (const ef_complex& u, const e_float& v) { return ef_complex(u) += v; }
  inline ef_complex operator- (const ef_complex& u, const e_float& v) { return ef_complex(u) -= v; }
  inline ef_complex operator* (const ef_complex& u, const e_float& v) { return ef_complex(u) *= v; }
  inline ef_complex operator/ (const ef_complex& u, const e_float& v) { return ef_complex(u) /= v; }

  inline ef_complex operator+ (const INT32 n, const ef_complex& v) { return ef_complex(e_float(n)) += v; }
  inline ef_complex operator- (const INT32 n, const ef_complex& v) { return ef_complex(e_float(n)) -= v; }
  inline ef_complex operator* (const INT32 n, const ef_complex& v) { return ef_complex(v)          *= n; }
  inline ef_complex operator/ (const INT32 n, const ef_complex& v) { return ef_complex(e_float(n)) /= v; }

  inline ef_complex operator+ (const ef_complex& z, const INT32 n) { return ef_complex(z) += n; }
  inline ef_complex operator- (const ef_complex& z, const INT32 n) { return ef_complex(z) -= n; }
  inline ef_complex operator* (const ef_complex& z, const INT32 n) { return ef_complex(z) *= n; }
  inline ef_complex operator/ (const ef_complex& z, const INT32 n) { return ef_complex(z) /= n; }

#endif // _ELEMENTARY_COMPLEX_2009_12_09_H_
