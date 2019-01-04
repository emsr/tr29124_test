// Synopsis of the ef_complex public interface.
  class ef_complex
  {
  private:
  public:

    explicit ef_complex(const  INT32 n);
    explicit ef_complex(const  INT64 n);
    explicit ef_complex(const UINT32 u);
    explicit ef_complex(const UINT64 u);
    explicit ef_complex(const double d);

    ef_complex(const e_float& re = ef::zero(), const e_float& im = ef::zero());

    ef_complex(const ef_complex& z);

    e_float real(void) const;
    e_float imag(void) const;

    static e_float real(const ef_complex& z);
    static e_float imag(const ef_complex& z);

    e_float norm(void) const;

    ef_complex& operator= (const ef_complex& v);
    ef_complex& operator= (const e_float& v);
    ef_complex& operator+=(const ef_complex& v);
    ef_complex& operator-=(const ef_complex& v);
    ef_complex& operator*=(const ef_complex& v);
    ef_complex& operator/=(const ef_complex& v);

    ef_complex& operator+=(const e_float& v);
    ef_complex& operator-=(const e_float& v);
    ef_complex& operator*=(const e_float& v);
    ef_complex& operator/=(const e_float& v);

    // Operators pre-increment and post-increment
    const ef_complex& operator++(void);
          ef_complex  operator++(int);

    // Operators pre-decrement and post-decrement
    const ef_complex& operator--(void);
          ef_complex  operator--(int);

    // Unary operators.
    ef_complex  operator-(void) const;
    ef_complex& operator+(void);

    // Operators with integer.
    ef_complex& operator+=(const INT32 n);
    ef_complex& operator-=(const INT32 n);
    ef_complex& operator*=(const INT32 n);
    ef_complex& operator/=(const INT32 n);

    bool isnan   (void) const;
    bool isinf   (void) const;
    bool isfinite(void) const;
    bool isneg   (void) const;
    bool ispos   (void) const;
    bool isint   (void) const;
    bool isone   (void) const;
    bool iszero  (void) const;

    // Utility get-functions to be used for Python exports.
    ef_complex  get_pos(void) const;
    ef_complex  get_neg(void) const;
    e_float     get_abs(void) const;
    INT64       get_int(void) const;
    double      get_flt(void) const;
    std::string get_str(void) const;
  };

  inline std::ostream& operator<<(std::ostream& os, const ef_complex& z) { return os << '(' << z.real() << ',' << z.imag() << ')'; }
