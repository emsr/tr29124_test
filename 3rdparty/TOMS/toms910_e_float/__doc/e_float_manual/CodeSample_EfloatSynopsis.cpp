// Synopsis of the e_float class.
namespace my_ef_space
{
  class e_float
  {
  private: // Private implementation details.
  public:  // Public interface.

    // Digit characteristics.
    static const INT64 ef_max_exp;
    static const INT64 ef_min_exp;
    static const INT64 ef_max_exp10;
    static const INT64 ef_min_exp10;

    // Constructors and destructor.
    e_float(void);
    explicit e_float(const INT32 n);
    explicit e_float(const INT64 n);
    explicit e_float(const UINT32 u);
    explicit e_float(const UINT64 u);
    explicit e_float(const double d);
    explicit e_float(const char* const s);
    explicit e_float(const std::string& str);
    e_float(const e_float& f);
    e_float(const double mantissa, const INT64 exponent);
    ~e_float();

    // Special values and precision.
    virtual const e_float& my_value_nan(void) const;
    virtual const e_float& my_value_inf(void) const;
    virtual const e_float& my_value_max(void) const;
    virtual const e_float& my_value_min(void) const;

    virtual void precision(const INT32 prec_digits);

    // Basic mathematical operations.
    virtual e_float& operator= (const e_float& v);
    virtual e_float& operator+=(const e_float& v);
    virtual e_float& operator-=(const e_float& v);
    virtual e_float& operator*=(const e_float& v);
    virtual e_float& operator/=(const e_float& v);

    virtual e_float& mul_by_int(const INT32 n);
    virtual e_float& div_by_int(const INT32 n);

    virtual e_float& calculate_inv (void);
    virtual e_float& calculate_sqrt(void);

    // Comparison functions.
    virtual INT32 cmp(const e_float& v) const;

    virtual bool isnan   (void) const;
    virtual bool isinf   (void) const;
    virtual bool isfinite(void) const;

    virtual bool iszero  (void) const;
    virtual bool isone   (void) const;
    virtual bool isint   (void) const;
    virtual bool isneg   (void) const;

    // Operators negate, post increment/decrement.
    virtual e_float& negate(void);
    virtual e_float& operator++(void);
    virtual e_float& operator--(void);

    // Extraction of the integer and decimal parts.
    virtual void    extract_parts       (double& mantissa, INT64& exponent) const;
    virtual double  extract_double      (void) const;
    virtual INT64   extract_int64       (void) const;
    virtual e_float extract_integer_part(void) const;
    virtual e_float extract_decimal_part(void) const;

    // The base-10 order of the e_float.
    virtual INT64 order(void) const;

    // The implementation of the "has_its_own"-mechanism.
    virtual bool has_its_own_cbrt (void) const;
    virtual bool has_its_own_rootn(void) const;
    virtual bool has_its_own_exp  (void) const;
    virtual bool has_its_own_log  (void) const;
    virtual bool has_its_own_sin  (void) const;
    virtual bool has_its_own_cos  (void) const;
    virtual bool has_its_own_tan  (void) const;
    virtual bool has_its_own_asin (void) const;
    virtual bool has_its_own_acos (void) const;
    virtual bool has_its_own_atan (void) const;
    virtual bool has_its_own_sinh (void) const;
    virtual bool has_its_own_cosh (void) const;
    virtual bool has_its_own_tanh (void) const;
    virtual bool has_its_own_asinh(void) const;
    virtual bool has_its_own_acosh(void) const;
    virtual bool has_its_own_atanh(void) const;

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
  };
}
