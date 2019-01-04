
#ifndef _E_FLOAT_BASE_2004_02_28_H_
  #define _E_FLOAT_BASE_2004_02_28_H_

  #if defined(__GNUC__)
    #include <tr1/array>
  #else
    #include <array>
  #endif

  #include <iostream>
  #include <string>

  #include <e_float/e_float_types.h>

  #define E_FLOAT_DIGITS10 100

  #if defined(E_FLOAT_TYPE_EFX)
    namespace efx { class e_float; }
    using efx::e_float;
  #elif defined(E_FLOAT_TYPE_GMP)
    namespace gmp { class e_float; }
    using gmp::e_float;
  #elif defined(E_FLOAT_TYPE_MPFR)
    namespace mpfr { class e_float; }
    using mpfr::e_float;
  #elif defined(E_FLOAT_TYPE_F90)
    namespace f90 { class e_float; }
    using f90::e_float;
    #undef E_FLOAT_DIGITS10
    #define E_FLOAT_DIGITS10 30
  #else
    #error e_float type is undefined! Define the e_float type!
  #endif

  class e_float_base
  {
  public:

    // The value of ef_digits10_setting is the desired number of decimal digits
    // of precision in the e_float implementation. It is limited to the range of
    // 30...300 decimal digits. The digit tolerance is invariant and it is set
    // to be 15 percent larger than ef_digits10. All of these quantities are
    // invariant and they are set at compile time.
    static const INT32 ef_digits10_setting = E_FLOAT_DIGITS10;
    static const INT32 ef_digits10         = ((ef_digits10_setting < 30) ? 30 : ((ef_digits10_setting > 300) ? 300 : ef_digits10_setting));
    static const INT32 ef_digits10_extra   = static_cast<INT32>(((ef_digits10 * 15) + 50) / 100);
    static const INT32 ef_digits10_tol     = static_cast<INT32>(ef_digits10 + ((ef_digits10_extra < 15) ? 15 : ef_digits10_extra));

    static const std::string::size_type& width_of_exponent_field(void);

    // This predominantly abstract base class has no constructors.

    // Virtual destructor.
    virtual ~e_float_base() { }

    // Specific special values.
    virtual const e_float_base& my_value_nan(void) const = 0;
    virtual const e_float_base& my_value_inf(void) const = 0;
    virtual const e_float_base& my_value_max(void) const = 0;
    virtual const e_float_base& my_value_min(void) const = 0;

    virtual INT32 cmp(const e_float& v) const = 0;

    virtual void precision(const INT32 prec_digits) = 0;

    // Basic operations.
    virtual e_float_base& operator= (const e_float& v) = 0;
    virtual e_float_base& operator+=(const e_float& v) = 0;
    virtual e_float_base& operator-=(const e_float& v) = 0;
    virtual e_float_base& operator*=(const e_float& v) = 0;
    virtual e_float_base& operator/=(const e_float& v) = 0;
    virtual e_float_base& mul_by_int(const INT32 n) = 0;
    virtual e_float_base& div_by_int(const INT32 n) = 0;

    virtual e_float_base& calculate_inv (void) = 0;
    virtual e_float_base& calculate_sqrt(void) = 0;

    // Comparison functions
    virtual bool isnan   (void) const = 0;
    virtual bool isinf   (void) const = 0;
    virtual bool isfinite(void) const = 0;

    virtual bool iszero (void) const = 0;
    virtual bool isone  (void) const = 0;
    virtual bool isint  (void) const = 0;
    virtual bool isneg  (void) const = 0;
            bool ispos  (void) const { return !isneg(); }

    virtual e_float_base& negate(void) = 0;

    // Operators pre-increment and pre-decrement
    virtual e_float_base& operator++(void) = 0;
    virtual e_float_base& operator--(void) = 0;

    // Argument range and check functions
    virtual INT64 order(void) const = 0;

    // Conversion routines
    virtual void    extract_parts       (double& mantissa, INT64& exponent) const = 0;
    virtual double  extract_double      (void) const = 0;
    virtual INT64   extract_int64       (void) const = 0;
    virtual e_float extract_integer_part(void) const = 0;
    virtual e_float extract_decimal_part(void) const = 0;

    // Formated Output routine.
    virtual void wr_string(std::string& str, std::ostream& os) const = 0;
    virtual bool rd_string(const char* const s) = 0;

    // Specific higher functions which might be present in the MP implementation.
    virtual bool has_its_own_cbrt         (void) const { return false; }
    virtual bool has_its_own_rootn        (void) const { return false; }
    virtual bool has_its_own_exp          (void) const { return false; }
    virtual bool has_its_own_log          (void) const { return false; }
    virtual bool has_its_own_sin          (void) const { return false; }
    virtual bool has_its_own_cos          (void) const { return false; }
    virtual bool has_its_own_tan          (void) const { return false; }
    virtual bool has_its_own_asin         (void) const { return false; }
    virtual bool has_its_own_acos         (void) const { return false; }
    virtual bool has_its_own_atan         (void) const { return false; }
    virtual bool has_its_own_sinh         (void) const { return false; }
    virtual bool has_its_own_cosh         (void) const { return false; }
    virtual bool has_its_own_tanh         (void) const { return false; }
    virtual bool has_its_own_asinh        (void) const { return false; }
    virtual bool has_its_own_acosh        (void) const { return false; }
    virtual bool has_its_own_atanh        (void) const { return false; }
    virtual bool has_its_own_gamma        (void) const { return false; }
    virtual bool has_its_own_riemann_zeta (void) const { return false; }
    virtual bool has_its_own_cyl_bessel_jn(void) const { return false; }
    virtual bool has_its_own_cyl_bessel_yn(void) const { return false; }

    static e_float my_cbrt         (const e_float& x);
    static e_float my_rootn        (const e_float& x, const UINT32 p);
    static e_float my_exp          (const e_float& x);
    static e_float my_log          (const e_float& x);
    static e_float my_sin          (const e_float& x);
    static e_float my_cos          (const e_float& x);
    static e_float my_tan          (const e_float& x);
    static e_float my_asin         (const e_float& x);
    static e_float my_acos         (const e_float& x);
    static e_float my_atan         (const e_float& x);
    static e_float my_sinh         (const e_float& x);
    static e_float my_cosh         (const e_float& x);
    static e_float my_tanh         (const e_float& x);
    static e_float my_asinh        (const e_float& x);
    static e_float my_acosh        (const e_float& x);
    static e_float my_atanh        (const e_float& x);
    static e_float my_gamma        (const e_float& x);
    static e_float my_riemann_zeta (const e_float& x);
    static e_float my_cyl_bessel_jn(const INT32 n, const e_float& x);
    static e_float my_cyl_bessel_yn(const INT32 n, const e_float& x);

    // Utility get-functions to be used for Python exports.
    e_float     get_pos(void) const;
    e_float     get_neg(void) const;
    e_float     get_abs(void) const;
    INT64       get_int(void) const;
    double      get_flt(void) const;
    std::string get_str(void) const;
  };

  std::ostream& operator<<(std::ostream& os, const e_float_base& f);
  std::istream& operator>>(std::istream& is, e_float_base& f);

#endif // _E_FLOAT_BASE_2004_02_28_H_
