
#ifndef _E_FLOAT_F90_2010_01_13_H_
  #define _E_FLOAT_F90_2010_01_13_H_

  #include <cmath>
  #include <string>

  #if defined(__GNUC__)
    #include <tr1/array>
  #else
    #include <array>
  #endif

  #include <e_float/e_float_base.h>
  #include <e_float/f90/e_float_f90_types.h>

  namespace f90
  {
    class e_float : public ::e_float_base
    {
    public:

      static const INT64 ef_max_exp   = static_cast<INT64>(+16383LL);
      static const INT64 ef_min_exp   = static_cast<INT64>(-16382LL);
      static const INT64 ef_max_exp10 = static_cast<INT64>(+4931LL);
      static const INT64 ef_min_exp10 = static_cast<INT64>(-4931LL);

    private:

      static const INT32 ef_digits2 = static_cast<INT32>(((ef_digits10_tol * 3322) + 500) / 1000);

    private:

      ::ef_f90_t data;

    public:

      e_float();
      explicit e_float(const INT32 n);
      explicit e_float(const INT64 n);
      explicit e_float(const UINT32 u);
      explicit e_float(const UINT64 u);
      explicit e_float(const double d);
      explicit e_float(const long double ld);
      explicit e_float(const char* s);
      explicit e_float(const std::string& str);

      e_float(const e_float& f);

      e_float(const double mantissa, const INT64 exponent);

      virtual ~e_float() { }

    private:

      static const double& d_log2(void);
      static bool char_is_nonzero_predicate(const char& c) { return (c != static_cast<char>('0')); }

    public:

      virtual INT32 cmp(const e_float& v) const;

      virtual const e_float& my_value_nan(void) const;
      virtual const e_float& my_value_inf(void) const;
      virtual const e_float& my_value_max(void) const;
      virtual const e_float& my_value_min(void) const;

      virtual void precision(const INT32 prec_digits) { static_cast<void>(prec_digits); }

      virtual e_float& operator= (const e_float& v) { data = v.data; return *this; }
      virtual e_float& operator+=(const e_float& v);
      virtual e_float& operator-=(const e_float& v);
      virtual e_float& operator*=(const e_float& v);
      virtual e_float& operator/=(const e_float& v);

      virtual e_float& mul_by_int(const INT32 n);
      virtual e_float& div_by_int(const INT32 n);

      virtual e_float& calculate_inv (void);
      virtual e_float& calculate_sqrt(void);

      virtual bool isnan   (void) const;
      virtual bool isinf   (void) const;
      virtual bool isfinite(void) const;

      virtual bool iszero  (void) const;
      virtual bool isone   (void) const;
      virtual bool isint   (void) const;
      virtual bool isneg   (void) const;

      virtual e_float& negate(void);

      virtual e_float& operator++(void);
      virtual e_float& operator--(void);

      virtual void    extract_parts       (double& mantissa, INT64& exponent) const;
      virtual double  extract_double      (void) const;
      virtual INT64   extract_int64       (void) const;
      virtual e_float extract_integer_part(void) const;
      virtual e_float extract_decimal_part(void) const;

      virtual INT64 order(void) const;

      static e_float my_cbrt        (const e_float& x);
      static e_float my_rootn       (const e_float& x, const UINT32 p);
      static e_float my_exp         (const e_float& x);
      static e_float my_log         (const e_float& x);
      static e_float my_sin         (const e_float& x);
      static e_float my_cos         (const e_float& x);
      static e_float my_tan         (const e_float& x);
      static e_float my_asin        (const e_float& x);
      static e_float my_acos        (const e_float& x);
      static e_float my_atan        (const e_float& x);
      static e_float my_sinh        (const e_float& x);
      static e_float my_cosh        (const e_float& x);
      static e_float my_tanh        (const e_float& x);
      static e_float my_asinh       (const e_float& x);
      static e_float my_acosh       (const e_float& x);
      static e_float my_atanh       (const e_float& x);

    private:

      virtual void wr_string(std::string& str, std::ostream& os) const;
      virtual bool rd_string(const char* const s);

    public:

      virtual bool has_its_own_cbrt (void) const { return true; }
      virtual bool has_its_own_rootn(void) const { return true; }
      virtual bool has_its_own_exp  (void) const { return true; }
      virtual bool has_its_own_log  (void) const { return true; }
      virtual bool has_its_own_sin  (void) const { return true; }
      virtual bool has_its_own_cos  (void) const { return true; }
      virtual bool has_its_own_tan  (void) const { return true; }
      virtual bool has_its_own_asin (void) const { return true; }
      virtual bool has_its_own_acos (void) const { return true; }
      virtual bool has_its_own_atan (void) const { return true; }
      virtual bool has_its_own_sinh (void) const { return true; }
      virtual bool has_its_own_cosh (void) const { return true; }
      virtual bool has_its_own_tanh (void) const { return true; }
      virtual bool has_its_own_asinh(void) const { return true; }
      virtual bool has_its_own_acosh(void) const { return true; }
      virtual bool has_its_own_atanh(void) const { return true; }
    };
  }

#endif // _E_FLOAT_F90_2010_01_13_H_
