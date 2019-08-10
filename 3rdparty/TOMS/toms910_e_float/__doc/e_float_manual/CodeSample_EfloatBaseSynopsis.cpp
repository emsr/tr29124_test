// Synopsis of the e_float_base class.
class e_float_base
{
public:  // Public interface.

  // Digit settings
  static const INT32 ef_digits10_setting;
  static const INT32 ef_digits10;
  static const INT32 ef_digits10_extra;
  static const INT32 ef_digits10_tol;

  static const std::string::size_type& width_of_exponent_field(void);

  // This predominantly abstract base class has no constructors.

  // Virtual destructor.
  virtual ~e_float_base() { }

  // Specific special values.
  virtual const e_float_base& my_value_nan(void) const = 0;
  virtual const e_float_base& my_value_inf(void) const = 0;
  virtual const e_float_base& my_value_max(void) const = 0;
  virtual const e_float_base& my_value_min(void) const = 0;

  virtual INT32 cmp(const e_float&) const = 0;

  virtual void precision(const INT32) = 0;

  // Basic operations.
  virtual e_float_base& operator= (const e_float&) = 0;
  virtual e_float_base& operator+=(const e_float&) = 0;
  virtual e_float_base& operator-=(const e_float&) = 0;
  virtual e_float_base& operator*=(const e_float&) = 0;
  virtual e_float_base& operator/=(const e_float&) = 0;
  virtual e_float_base& mul_by_int(const INT32) = 0;
  virtual e_float_base& div_by_int(const INT32) = 0;

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
  virtual void    extract_parts       (double&, INT64&) const = 0;
  virtual double  extract_double      (void) const = 0;
  virtual INT64   extract_int64       (void) const = 0;
  virtual e_float extract_integer_part(void) const = 0;
  virtual e_float extract_decimal_part(void) const = 0;

  // Formated input and output routines.
  virtual void wr_string(std::string&, std::ostream&) const = 0;
  virtual bool rd_string(const char* const) = 0;

    // The implementation of the "has_its_own"-mechanism.
  virtual bool has_its_own_cbrt         (void) const;
  virtual bool has_its_own_rootn        (void) const;
  virtual bool has_its_own_exp          (void) const;
  virtual bool has_its_own_log          (void) const;
  virtual bool has_its_own_sin          (void) const;
  virtual bool has_its_own_cos          (void) const;
  virtual bool has_its_own_tan          (void) const;
  virtual bool has_its_own_asin         (void) const;
  virtual bool has_its_own_acos         (void) const;
  virtual bool has_its_own_atan         (void) const;
  virtual bool has_its_own_sinh         (void) const;
  virtual bool has_its_own_cosh         (void) const;
  virtual bool has_its_own_tanh         (void) const;
  virtual bool has_its_own_asinh        (void) const;
  virtual bool has_its_own_acosh        (void) const;
  virtual bool has_its_own_atanh        (void) const;
  virtual bool has_its_own_gamma        (void) const;
  virtual bool has_its_own_riemann_zeta (void) const;
  virtual bool has_its_own_cyl_bessel_jn(void) const;
  virtual bool has_its_own_cyl_bessel_yn(void) const;

  static e_float my_cbrt         (const e_float&);
  static e_float my_rootn        (const e_float&, const UINT32);
  static e_float my_exp          (const e_float&);
  static e_float my_log          (const e_float&);
  static e_float my_sin          (const e_float&);
  static e_float my_cos          (const e_float&);
  static e_float my_tan          (const e_float&);
  static e_float my_asin         (const e_float&);
  static e_float my_acos         (const e_float&);
  static e_float my_atan         (const e_float&);
  static e_float my_sinh         (const e_float&);
  static e_float my_cosh         (const e_float&);
  static e_float my_tanh         (const e_float&);
  static e_float my_asinh        (const e_float&);
  static e_float my_acosh        (const e_float&);
  static e_float my_atanh        (const e_float&);
  static e_float my_gamma        (const e_float&);
  static e_float my_riemann_zeta (const e_float&);
  static e_float my_cyl_bessel_jn(const INT32, const e_float&);
  static e_float my_cyl_bessel_yn(const INT32, const e_float&);

  // Utility get-functions to be used for Python exports.
  e_float     get_pos(void) const;
  e_float     get_neg(void) const;
  e_float     get_abs(void) const;
  INT64       get_int(void) const;
  double      get_flt(void) const;
  std::string get_str(void) const;
};

// Global operators with iostream objects.
std::ostream& operator<<(std::ostream& os, const e_float_base& f);
std::istream& operator>>(std::istream& is, e_float_base& f);
