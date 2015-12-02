
#include <limits>

namespace std
{

  /**
   *  @brief Part of std::numeric_limits.
   */

  template<typename _Tp>
    constexpr bool
    is_specialized(_Tp) noexcept
    { return numeric_limits<_Tp>::is_specialized; }

  template<typename _Tp>
    constexpr _Tp
    min(_Tp) noexcept
    { return numeric_limits<_Tp>::min(); }

  template<typename _Tp>
    constexpr _Tp
    max(_Tp) noexcept
    { return numeric_limits<_Tp>::max(); }

  template<typename _Tp>
    static constexpr _Tp
    lowest(_Tp) noexcept
    { return numeric_limits<_Tp>::lowest(); }

  template<typename _Tp>
    constexpr int
    digits(_Tp) noexcept
    { return numeric_limits<_Tp>::digits; }

  template<typename _Tp>
    constexpr int
    digits10(_Tp) noexcept
    { return numeric_limits<_Tp>::digits10; }

  template<typename _Tp>
    constexpr int
    max_digits10(_Tp) noexcept
    { return numeric_limits<_Tp>::max_digits10; }

  template<typename _Tp>
    constexpr bool
    xxx_is_signed(_Tp) noexcept
    { return numeric_limits<_Tp>::is_signed; }

  template<typename _Tp>
    constexpr bool
    is_integer(_Tp) noexcept
    { return numeric_limits<_Tp>::is_integer; }

  template<typename _Tp>
    constexpr bool
    is_exact(_Tp) noexcept
    { return numeric_limits<_Tp>::is_exact; }

  template<typename _Tp>
    constexpr int
    radix(_Tp) noexcept
    { return numeric_limits<_Tp>::radix; }

  template<typename _Tp>
    constexpr _Tp
    epsilon(_Tp) noexcept
    { return numeric_limits<_Tp>::epsilon(); }

  template<typename _Tp>
    constexpr _Tp
    round_error(_Tp) noexcept
    { return numeric_limits<_Tp>::round_error(); }

  template<typename _Tp>
    constexpr int
    min_exponent(_Tp) noexcept
    { return numeric_limits<_Tp>::min_exponent; }

  template<typename _Tp>
    constexpr int
    min_exponent10(_Tp) noexcept
    { return numeric_limits<_Tp>::min_exponent10; }

  template<typename _Tp>
    constexpr int
    max_exponent(_Tp) noexcept
    { return numeric_limits<_Tp>::max_exponent; }

  template<typename _Tp>
    constexpr int
    max_exponent10(_Tp) noexcept
    { return numeric_limits<_Tp>::max_exponent10; }

  template<typename _Tp>
    constexpr bool
    has_infinity(_Tp) noexcept
    { return numeric_limits<_Tp>::has_infinity; }

  template<typename _Tp>
    constexpr bool
    has_quiet_NaN(_Tp) noexcept
    { return numeric_limits<_Tp>::has_quiet_NaN; }

  template<typename _Tp>
    constexpr bool
    has_signaling_NaN(_Tp) noexcept
    { return numeric_limits<_Tp>::has_signaling_NaN; }

  template<typename _Tp>
    constexpr float_denorm_style
    has_denorm(_Tp) noexcept
    { return numeric_limits<_Tp>::has_denorm; }

  template<typename _Tp>
    constexpr bool
    has_denorm_loss(_Tp) noexcept
    { return numeric_limits<_Tp>::has_denorm_loss; }

  template<typename _Tp>
    constexpr _Tp
    infinity(_Tp) noexcept
    { return numeric_limits<_Tp>::infinity(); }

  template<typename _Tp>
    constexpr _Tp
    quiet_NaN(_Tp) noexcept
    { return numeric_limits<_Tp>::quiet_NaN(); }

  template<typename _Tp>
    constexpr _Tp
    signaling_NaN(_Tp) noexcept
    { return numeric_limits<_Tp>::signaling_NaN(); }

  template<typename _Tp>
    constexpr _Tp
    denorm_min(_Tp) noexcept
    { return numeric_limits<_Tp>::denorm_min(); }

  template<typename _Tp>
    constexpr bool
    is_iec559(_Tp) noexcept
    { return numeric_limits<_Tp>::is_iec559; }

  template<typename _Tp>
    constexpr bool
    is_bounded(_Tp) noexcept
    { return numeric_limits<_Tp>::is_bounded; }

  template<typename _Tp>
    constexpr bool
    is_modulo(_Tp) noexcept
    { return numeric_limits<_Tp>::is_modulo; }

  template<typename _Tp>
    constexpr bool
    traps(_Tp) noexcept
    { return numeric_limits<_Tp>::traps; }

  template<typename _Tp>
    constexpr bool
    tinyness_before(_Tp) noexcept
    { return numeric_limits<_Tp>::tinyness_before; }

  template<typename _Tp>
    constexpr float_round_style
    round_style(_Tp) noexcept
    { return numeric_limits<_Tp>::round_style; }

}
