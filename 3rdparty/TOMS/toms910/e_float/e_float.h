
#ifndef _E_FLOAT_2004_02_28_H_
  #define _E_FLOAT_2004_02_28_H_

  #include <limits>

  #if defined(E_FLOAT_TYPE_EFX)
    #include <e_float/efx/e_float_efx.h>
    using efx::e_float;
  #elif defined(E_FLOAT_TYPE_GMP)
    #include <e_float/gmp/e_float_gmp.h>
    using gmp::e_float;
  #elif defined(E_FLOAT_TYPE_MPFR)
    #include <e_float/mpfr/e_float_mpfr.h>
    using mpfr::e_float;
  #elif defined(E_FLOAT_TYPE_F90)
    #include <e_float/f90/e_float_f90.h>
    using f90::e_float;
  #else
    #error e_float type undefined!
  #endif

  // Global operators post-increment and post-decrement
  inline e_float operator++(e_float& u, int) { const e_float v(u); ++u; return v; }
  inline e_float operator--(e_float& u, int) { const e_float v(u); --u; return v; }

  // Global unary operators of e_float reference.
  inline       e_float  operator-(const e_float& u) { return e_float(u).negate(); }
  inline       e_float& operator+(      e_float& u) { return u; }
  inline const e_float& operator+(const e_float& u) { return u; }

  // Global add/sub/mul/div of const e_float reference with const e_float reference
  inline e_float operator+(const e_float& u, const e_float& v) { return e_float(u) += v; }
  inline e_float operator-(const e_float& u, const e_float& v) { return e_float(u) -= v; }
  inline e_float operator*(const e_float& u, const e_float& v) { return e_float(u) *= v; }
  inline e_float operator/(const e_float& u, const e_float& v) { return e_float(u) /= v; }

  // Specialization for global add/sub/mul/div of const e_float reference with INT32
  inline e_float operator+(const e_float& u, const INT32 n) { return (e_float(u) += e_float(n)); }
  inline e_float operator-(const e_float& u, const INT32 n) { return (e_float(u) -= e_float(n)); }
  inline e_float operator*(const e_float& u, const INT32 n) { return  e_float(u).mul_by_int(n); }
  inline e_float operator/(const e_float& u, const INT32 n) { return  e_float(u).div_by_int(n); }

  inline e_float operator+(const INT32 n, const e_float& u) { return (e_float(n) += u); }
  inline e_float operator-(const INT32 n, const e_float& u) { return (e_float(n) -= u); }
  inline e_float operator*(const INT32 n, const e_float& u) { return (e_float(n) *= u); }
  inline e_float operator/(const INT32 n, const e_float& u) { return (e_float(n) /= u); }

  // Specializations of global self-add/sub/mul-div of e_float reference with INT32
  inline e_float& operator+=(e_float& u, const INT32 n) { return (u += e_float(n)); }
  inline e_float& operator-=(e_float& u, const INT32 n) { return (u -= e_float(n)); }
  inline e_float& operator*=(e_float& u, const INT32 n) { return u.mul_by_int(n); }
  inline e_float& operator/=(e_float& u, const INT32 n) { return u.div_by_int(n); }

  // Global comparison operators of const e_float reference with const e_float reference
  inline bool operator< (const e_float& u, const e_float& v) { return (u.cmp(v) <  static_cast<INT32>(0)); }
  inline bool operator<=(const e_float& u, const e_float& v) { return (u.cmp(v) <= static_cast<INT32>(0)); }
  inline bool operator==(const e_float& u, const e_float& v) { return (u.cmp(v) == static_cast<INT32>(0)); }
  inline bool operator!=(const e_float& u, const e_float& v) { return (u.cmp(v) != static_cast<INT32>(0)); }
  inline bool operator>=(const e_float& u, const e_float& v) { return (u.cmp(v) >= static_cast<INT32>(0)); }
  inline bool operator> (const e_float& u, const e_float& v) { return (u.cmp(v) >  static_cast<INT32>(0)); }

  // Specializations of global comparison of e_float reference with INT32
  inline bool operator< (const e_float& u, const INT32 n) { return u <  e_float(n); }
  inline bool operator<=(const e_float& u, const INT32 n) { return u <= e_float(n); }
  inline bool operator==(const e_float& u, const INT32 n) { return u == e_float(n); }
  inline bool operator!=(const e_float& u, const INT32 n) { return u != e_float(n); }
  inline bool operator>=(const e_float& u, const INT32 n) { return u >= e_float(n); }
  inline bool operator> (const e_float& u, const INT32 n) { return u >  e_float(n); }

  inline bool operator< (const INT32 n, const e_float& u) { return e_float(n) <  u; }
  inline bool operator<=(const INT32 n, const e_float& u) { return e_float(n) <= u; }
  inline bool operator==(const INT32 n, const e_float& u) { return e_float(n) == u; }
  inline bool operator!=(const INT32 n, const e_float& u) { return e_float(n) != u; }
  inline bool operator>=(const INT32 n, const e_float& u) { return e_float(n) >= u; }
  inline bool operator> (const INT32 n, const e_float& u) { return e_float(n) >  u; }

  namespace ef
  {
    const e_float& zero     (void);
    const e_float& one      (void);
    const e_float& half     (void);
    const e_float& value_min(void);
    const e_float& value_max(void);
    const e_float& value_eps(void);
    const e_float& value_inf(void);
    const e_float& value_nan(void);

    inline INT64 tol(void) { return static_cast<INT64>(e_float::ef_digits10_tol); }

    inline e_float fabs(const e_float& x) { return (x.isneg() ? e_float(x).negate() : x); }
  }

  // Specialization of std::numeric_limits<e_float>.
  namespace std
  {
    template <> class numeric_limits<e_float>
    {
    public: // Implement the "usual" public members for floating point types.
      static const bool  is_specialized    = true;
      static const bool  is_signed         = true;
      static const bool  is_integer        = false;
      static const bool  is_exact          = false;
      static const bool  is_bounded        = true;
      static const bool  is_modulo         = false;
      static const bool  is_iec559         = false;
      static const int   digits10          = e_float::ef_digits10;
      static const int   digits            = ((digits10 * 301) + 500) / 1000;
      static const INT64 max_exponent      = e_float::ef_max_exp;
      static const INT64 min_exponent      = e_float::ef_min_exp;
      static const INT64 max_exponent10    = e_float::ef_max_exp10;
      static const INT64 min_exponent10    = e_float::ef_min_exp10;
      static const int   radix             = 10;
      static const int   round_style       = std::round_to_nearest;
      static const bool  has_infinity      = true;
      static const bool  has_quiet_NaN     = true;
      static const bool  has_signaling_NaN = false;
      static const int   has_denorm        = std::denorm_absent;
      static const bool  has_denorm_loss   = false;
      static const bool  traps             = false;
      static const bool  tinyness_before   = false;

      static const e_float& (min)      (void) throw() { return ef::value_min(); }
      static const e_float& (max)      (void) throw() { return ef::value_max(); }
      static const e_float& epsilon    (void) throw() { return ef::value_eps(); }
      static const e_float& round_error(void) throw() { return ef::half(); }
      static const e_float& infinity   (void) throw() { return ef::value_inf(); }
      static const e_float& quiet_NaN  (void) throw() { return ef::value_nan(); }
    };
  }

#endif // _E_FLOAT_2004_02_28_H_
