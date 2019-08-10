
#ifndef _ELEMENTARY_MATH_2009_12_09_H_
  #define _ELEMENTARY_MATH_2009_12_09_H_

  #include <limits>

  class ef_complex;

  namespace ef
  {
    inline INT32 max_iteration(void) { return static_cast<INT32>(10000); }

    e_float floor(const e_float& x);
    e_float ceil (const e_float& x);
    INT32   sgn  (const e_float& x);

           bool isnan(const double x);
    inline bool isnan(const e_float& x)        { return x.isnan(); }
    inline bool isnan(const ef_complex& z)     { return z.isnan(); }

           bool isfinite(const double x);
    inline bool isfinite(const e_float& x)     { return x.isfinite(); }
    inline bool isfinite(const ef_complex& z)  { return z.isfinite(); }

           bool isinf(const double x);
    inline bool isinf(const e_float& x)        { return x.isinf(); }
    inline bool isinf(const ef_complex& z)     { return z.isinf(); }

    inline bool isneg(const double x)          { return x < static_cast<double>(0.0); }
    inline bool isneg(const e_float& x)        { return x.isneg(); }
    inline bool isneg(const ef_complex& z)     { return z.isneg(); }

    inline e_float abs (const e_float& x)      { return ef::fabs(x); }
    inline e_float real(const e_float& x)      { return x; }
    inline e_float imag(const e_float& x)      { static_cast<void>(x); return ef::zero(); }

    inline bool ispos(const double x)          { return !isneg(x); }
    inline bool ispos(const e_float& x)        { return !x.isneg(); }
    inline bool ispos(const ef_complex& z)     { return !z.isneg(); }

           bool isint(const double x);
    inline bool isint(const e_float& x)        { return x.isint(); }
    inline bool isint(const ef_complex& z)     { return z.isint(); }

    inline bool isone(const double x)          { return ::fabs(static_cast<double>(1.0) - x) < std::numeric_limits<double>::min() * static_cast<double>(2.0); }
    inline bool isone(const e_float& x)        { return x.isone(); }
    inline bool isone(const ef_complex& z)     { return z.isone(); }

    inline bool iszero(const double x)         { return ::fabs(x) < std::numeric_limits<double>::min() * static_cast<double>(2.0); }
    inline bool iszero(const e_float& x)       { return x.iszero(); }
    inline bool iszero(const ef_complex& z)    { return z.iszero(); }

    e_float integer_part(const e_float& x);
    e_float decimal_part(const e_float& x);

    void to_parts(const e_float& x, double& mantissa, INT64& exponent);

    inline double to_double(const double& x)   { return x; }
           double to_double(const e_float& x);
           double to_double(const ef_complex& z);

    INT64 to_int64(const double x);
    INT64 to_int64(const e_float& x);
    INT64 to_int64(const ef_complex& z);

    INT32 to_int32(const double x);
    INT32 to_int32(const e_float& x);
    INT32 to_int32(const ef_complex& z);

    bool small_arg(const double x);
    bool small_arg(const e_float& x);
    bool small_arg(const ef_complex& z);

    bool large_arg(const double x);
    bool large_arg(const e_float& x);
    bool large_arg(const ef_complex& z);

    bool near_one(const double x);
    bool near_one(const e_float& x);
    bool near_one(const ef_complex& z);

    bool near_int(const double x);
    bool near_int(const e_float& x);
    bool near_int(const ef_complex& z);
  }

#endif // _ELEMENTARY_MATH_2009_12_09_H_
