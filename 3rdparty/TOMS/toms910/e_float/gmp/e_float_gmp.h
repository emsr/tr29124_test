
#ifndef _GMP_E_FLOAT_2008_05_01_H_
  #define _GMP_E_FLOAT_2008_05_01_H_

  #if defined(__INTEL_COMPILER)
    #pragma warning (disable:193)
    #pragma warning (disable:981)
  #endif

  #if defined(_MSC_VER)
    #pragma warning (disable:4127)
  #endif

  #include <cmath>
  #include <string>

  #if defined(__GNUC__)
    #include <tr1/array>
  #else
    #include <array>
  #endif

  #include <e_float/e_float_base.h>
  #include <e_float/gmp/e_float_gmp_types.h>

  namespace gmp
  {
    class e_float : public ::e_float_base
    {
    public:

      static const INT64 ef_max_exp   = static_cast<INT64>(+2147483640LL);
      static const INT64 ef_min_exp   = static_cast<INT64>(-2147483640LL);
      static const INT64 ef_max_exp10 = static_cast<INT64>(+646456990LL);
      static const INT64 ef_min_exp10 = static_cast<INT64>(-646456990LL);

    private:

      static const INT32 ef_digits2 = static_cast<INT32>(((ef_digits10_tol * 3322) + 500) / 1000);

      typedef enum enum_fpclass
      {
        ef_finite,
        ef_inf,
        ef_NaN
      }
      t_fpclass;

      t_fpclass fpclass;
      INT32     prec_elem;
      ::mpf_t   rop;

    public:
    
      // Constructors
      e_float(void);

      explicit e_float(const  INT32 n);
      explicit e_float(const  INT64 n);
      explicit e_float(const UINT32 u);
      explicit e_float(const UINT64 u);
      explicit e_float(const double d);
      explicit e_float(const long double ld);
      explicit e_float(const char* const s);
      explicit e_float(const std::string& str);

      e_float(const e_float& f);

      // Constructor from mantissa and exponent.
      e_float(const double mantissa, const INT64 exponent);

      virtual ~e_float();

    private:
    
      explicit e_float(const ::mpf_t& op);

    private:

      static const double& d_log2  (void);
      static const INT64&  max_exp2(void);
      static const INT64&  min_exp2(void);

      static void init(void);

      void from_uint64(const UINT64 u);
      void from_uint32(const UINT32 u);

      INT32 cmp_data(const ::mpf_t& v) const;

    public:

      virtual INT32 cmp(const e_float& v) const;

      virtual const e_float& my_value_nan(void) const;
      virtual const e_float& my_value_inf(void) const;
      virtual const e_float& my_value_max(void) const;
      virtual const e_float& my_value_min(void) const;

      virtual void precision(const INT32 prec_digits);

      // Basic operations.
      virtual e_float& operator= (const e_float& v);
      virtual e_float& operator+=(const e_float& v);
      virtual e_float& operator-=(const e_float& v);
      virtual e_float& operator*=(const e_float& v);
      virtual e_float& operator/=(const e_float& v);
      virtual e_float& mul_by_int(const INT32 n);
      virtual e_float& div_by_int(const INT32 n);

      virtual e_float& calculate_inv (void);
      virtual e_float& calculate_sqrt(void);

      // Comparison functions
      virtual bool isnan   (void) const { return fpclass == ef_NaN; }
      virtual bool isinf   (void) const { return fpclass == ef_inf; }
      virtual bool isfinite(void) const { return fpclass == ef_finite; }

      virtual bool iszero (void) const;
      virtual bool isone  (void) const;
      virtual bool isint  (void) const;
      virtual bool isneg  (void) const;

      virtual e_float& negate(void);

      // Operators pre-increment and pre-decrement
      virtual e_float& operator++(void);
      virtual e_float& operator--(void);

      // Conversion routines
      virtual void    extract_parts       (double& mantissa, INT64& exponent) const;
      virtual double  extract_double      (void) const;
      virtual INT64   extract_int64       (void) const;
      virtual e_float extract_integer_part(void) const;
      virtual e_float extract_decimal_part(void) const;

      // Argument range and check functions
      virtual INT64 order(void) const;

    private:

      virtual void wr_string(std::string& str, std::ostream& os) const;
      virtual bool rd_string(const char* const s);
    };
  }

#endif // _GMP_E_FLOAT_2008_05_01_H_
