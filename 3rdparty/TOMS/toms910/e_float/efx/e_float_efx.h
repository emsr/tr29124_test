#ifndef _E_FLOAT_EFX_2004_02_28_H_
  #define _E_FLOAT_EFX_2004_02_28_H_

  #if defined(__INTEL_COMPILER)
    #pragma warning (disable:981)
  #endif

  #include <cmath>
  #include <string>

  #include <e_float/e_float_base.h>
  #include <e_float/efx/e_float_efx_array.h>

  namespace efx
  {
    class e_float : public ::e_float_base
    {
    public:

      static const INT64 ef_max_exp   = static_cast<INT64>(+9223372036854775795LL);
      static const INT64 ef_min_exp   = static_cast<INT64>(-9223372036854775795LL);
      static const INT64 ef_max_exp10 = static_cast<INT64>(+3063937869882635616LL); // Approx. [ef_max_exp / log10(2)], also an even multiple of 8
      static const INT64 ef_min_exp10 = static_cast<INT64>(-3063937869882635616LL);

    private:

      static const INT32 ef_elem_digits10     = static_cast<INT32>(8);
      static const INT32 ef_digits10_num_base = static_cast<INT32>((ef_digits10_tol / ef_elem_digits10) + (((ef_digits10_tol % ef_elem_digits10) != 0) ? 1 : 0));
      static const INT32 ef_elem_number       = static_cast<INT32>(ef_digits10_num_base + 2);

    private:

      typedef enum enum_fpclass
      {
        ef_finite,
        ef_inf,
        ef_NaN
      }
      t_fpclass;

      static const INT32 ef_elem_mask = static_cast<INT32>(100000000);

      typedef efx::array<UINT32, static_cast<std::size_t>(ef_elem_number)> array_type;

      array_type data;
      INT64      exp;
      bool       neg;
      t_fpclass  fpclass;
      INT32      prec_elem;

    public:

      // Constructors
      e_float();

      explicit e_float(const  INT32 n) : neg(n < static_cast<INT32>(0)) { from_uint32(UINT32(n < static_cast<INT32>(0) ? -n : n)); }
      explicit e_float(const  INT64 n) : neg(n < static_cast<INT64>(0)) { from_uint64(UINT64(n < static_cast<INT64>(0) ? -n : n)); }
      explicit e_float(const UINT32 u) : neg(false) { from_uint32(u); }
      explicit e_float(const UINT64 u) : neg(false) { from_uint64(u); }

      explicit e_float(const double d);
      explicit e_float(const long double ld);

      explicit e_float(const char* const s);
      explicit e_float(const std::string& str);

      e_float(const e_float& f) : data     (f.data),
                                  exp      (f.exp),
                                  neg      (f.neg),
                                  fpclass  (f.fpclass),
                                  prec_elem(f.prec_elem) { }

      // Constructor from mantissa and exponent.
      e_float(const double mantissa, const INT64 exponent);

      virtual ~e_float() { }

    private:

      static bool data_elem_is_nonzero_predicate(const UINT32& d) { return d != static_cast<UINT32>(0u); }
      static bool char_is_nonzero_predicate     (const char& c)   { return (c != static_cast<char>('0')); }

      void from_uint64(const UINT64 u);
      void from_uint32(const UINT32 u);

      INT32 cmp_data(const array_type& vd) const;

      static void mul_loop_uv(const UINT32* const u, const UINT32* const v, UINT32* const w, const INT32 p);
      static UINT32 mul_loop_n(UINT32* const u, UINT32 n, const INT32 p);
      static UINT32 div_loop_n(UINT32* const u, UINT32 n, const INT32 p);

    public:

      virtual INT32 cmp(const e_float& v) const;

      virtual const e_float& my_value_nan(void) const;
      virtual const e_float& my_value_inf(void) const;
      virtual const e_float& my_value_max(void) const;
      virtual const e_float& my_value_min(void) const;

      virtual void precision(const INT32 prec_digits)
      {
        const INT32 elems = static_cast<INT32>(  static_cast<INT32>((prec_digits + (ef_elem_digits10 / 2)) / ef_elem_digits10)
                                               + static_cast<INT32>((prec_digits % ef_elem_digits10) != 0 ? 1 : 0));

        prec_elem = elems > ef_elem_number ? ef_elem_number : elems;
      }

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
      virtual bool isneg  (void) const { return neg; }

      virtual e_float& negate(void) { if(!iszero()) { neg = !neg; } return *this; }

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
      virtual INT64 order(void) const { return (iszero() ? static_cast<INT64>(0)
                                                         : static_cast<INT64>(exp + static_cast<INT64>(::log10(static_cast<double>(data[0])) + 0.5))); }

    private:

      virtual void wr_string(std::string& str, std::ostream& os) const;
      virtual bool rd_string(const char* const s);
    };
  }

#endif // _E_FLOAT_EFX_2004_02_28_H_
