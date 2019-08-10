// Specialization of std::numeric_limits<e_float>.
namespace std
{
  template <> class numeric_limits<e_float>
  {
  public: // Implement some "usual" public members for floating point types.
    static const bool  is_specialized;
    static const bool  is_signed;
    static const bool  is_integer;
    static const bool  is_exact;
    static const bool  is_bounded;
    static const bool  is_modulo;
    static const bool  is_iec559;
    static const int   digits10;
    static const int   digits;
    static const INT64 max_exponent;
    static const INT64 min_exponent;
    static const INT64 max_exponent10
    static const INT64 min_exponent10;
    static const int   radix;
    static const int   round_style;
    static const bool  has_infinity;
    static const bool  has_quiet_NaN;
    static const bool  has_signaling_NaN;
    static const int   has_denorm;
    static const bool  has_denorm_loss;
    static const bool  traps;
    static const bool  tinyness_before;

    static const e_float& (min)      (void) throw();
    static const e_float& (max)      (void) throw();
    static const e_float& epsilon    (void) throw();
    static const e_float& round_error(void) throw();
    static const e_float& infinity   (void) throw();
    static const e_float& quiet_NaN  (void) throw();
  };
}
