
#ifndef NUMERIC_LIMITS_MPREAL_H
#define NUMERIC_LIMITS_MPREAL_H 1

#include <limits>
#include <mpreal.h>
#include <math_mpreal.h>
#include <emsr/numeric_limits.h>

namespace emsr
{

  /**
   *  @brief Part of std::numeric_limits.
   *  The idea is that types with, say, non-constexpr or even dynamic epsilon()
   *  can participate in this.
   *  I think variable templates could be specialized with non-constexpr types
   *  but I need something to work in C++11 and variable templates won't allow
   *  extraction of variable max from a mp number.
   */

  // Constexpr function template versions of std::numeric_limits.

  template<>
    inline bool
    is_specialized<mpfr::mpreal>(mpfr::mpreal) noexcept
    { return std::numeric_limits<mpfr::mpreal>::is_specialized; }

  template<>
    inline mpfr::mpreal
    lim_min<mpfr::mpreal>(mpfr::mpreal proto) noexcept
    { return mpfr::minval(proto.getPrecision()); }

  template<>
    inline mpfr::mpreal
    lim_max<mpfr::mpreal>(mpfr::mpreal proto) noexcept
    { return mpfr::maxval(proto.getPrecision()); }

  template<>
    inline mpfr::mpreal
    lowest<mpfr::mpreal>(mpfr::mpreal proto) noexcept
    { return -mpfr::maxval(proto.getPrecision()); }

  template<>
    inline int
    digits<mpfr::mpreal>(mpfr::mpreal proto) noexcept
    { return proto.getPrecision(); }

  template<>
    inline int
    digits10<mpfr::mpreal>(mpfr::mpreal proto) noexcept
    { return mpfr::bits2digits(proto.getPrecision()); }

  template<>
    inline int
    max_digits10<mpfr::mpreal>(mpfr::mpreal proto) noexcept
    { return digits10(proto); }

  template<>
    inline bool
    is_signed<mpfr::mpreal>(mpfr::mpreal) noexcept
    { return std::numeric_limits<mpfr::mpreal>::is_signed; }

  template<>
    inline bool
    is_integer<mpfr::mpreal>(mpfr::mpreal) noexcept
    { return std::numeric_limits<mpfr::mpreal>::is_integer; }

  template<>
    inline bool
    is_exact<mpfr::mpreal>(mpfr::mpreal) noexcept
    { return std::numeric_limits<mpfr::mpreal>::is_exact; }

  template<>
    inline int
    radix<mpfr::mpreal>(mpfr::mpreal) noexcept
    { return std::numeric_limits<mpfr::mpreal>::radix; }

  template<>
    inline mpfr::mpreal
    epsilon<mpfr::mpreal>(mpfr::mpreal proto) noexcept
    { return  mpfr::machine_epsilon(proto.getPrecision()); }

  template<>
    inline mpfr::mpreal
    round_error<mpfr::mpreal>(mpfr::mpreal proto) noexcept
    { return std::numeric_limits<mpfr::mpreal>::round_error(proto.getPrecision()); }

  template<>
    inline int
    min_exponent<mpfr::mpreal>(mpfr::mpreal) noexcept
    { return std::numeric_limits<mpfr::mpreal>::min_exponent; }

  template<>
    inline int
    min_exponent10<mpfr::mpreal>(mpfr::mpreal) noexcept
    { return std::numeric_limits<mpfr::mpreal>::min_exponent10; }

  template<>
    inline int
    max_exponent<mpfr::mpreal>(mpfr::mpreal) noexcept
    { return std::numeric_limits<mpfr::mpreal>::max_exponent; }

  template<>
    inline int
    max_exponent10<mpfr::mpreal>(mpfr::mpreal) noexcept
    { return std::numeric_limits<mpfr::mpreal>::max_exponent10; }

  template<>
    inline bool
    has_infinity<mpfr::mpreal>(mpfr::mpreal) noexcept
    { return std::numeric_limits<mpfr::mpreal>::has_infinity; }

  template<>
    inline bool
    has_quiet_NaN<mpfr::mpreal>(mpfr::mpreal) noexcept
    { return std::numeric_limits<mpfr::mpreal>::has_quiet_NaN; }

  template<>
    inline bool
    has_signaling_NaN<mpfr::mpreal>(mpfr::mpreal) noexcept
    { return std::numeric_limits<mpfr::mpreal>::has_signaling_NaN; }

  template<>
    inline std::float_denorm_style
    has_denorm<mpfr::mpreal>(mpfr::mpreal) noexcept
    { return std::numeric_limits<mpfr::mpreal>::has_denorm; }

  // No std::numeric_limits<mpfr::mpreal>::has_denorm_loss.
  template<>
    inline bool
    has_denorm_loss<mpfr::mpreal>(mpfr::mpreal) noexcept
    { return std::numeric_limits<mpfr::mpreal>::has_denorm; }

  template<>
    inline mpfr::mpreal
    infinity<mpfr::mpreal>(mpfr::mpreal) noexcept
    { return std::numeric_limits<mpfr::mpreal>::infinity(); }

  template<>
    inline mpfr::mpreal
    quiet_NaN<mpfr::mpreal>(mpfr::mpreal) noexcept
    { return std::numeric_limits<mpfr::mpreal>::quiet_NaN(); }

  template<>
    inline mpfr::mpreal
    signaling_NaN<mpfr::mpreal>(mpfr::mpreal) noexcept
    { return std::numeric_limits<mpfr::mpreal>::signaling_NaN(); }

  template<>
    inline mpfr::mpreal
    denorm_min<mpfr::mpreal>(mpfr::mpreal) noexcept
    { return std::numeric_limits<mpfr::mpreal>::denorm_min(); }

  template<>
    inline bool
    is_iec559<mpfr::mpreal>(mpfr::mpreal) noexcept
    { return std::numeric_limits<mpfr::mpreal>::is_iec559; }

  template<>
    inline bool
    is_bounded<mpfr::mpreal>(mpfr::mpreal) noexcept
    { return std::numeric_limits<mpfr::mpreal>::is_bounded; }

  template<>
    inline bool
    is_modulo<mpfr::mpreal>(mpfr::mpreal) noexcept
    { return std::numeric_limits<mpfr::mpreal>::is_modulo; }

  template<>
    inline bool
    traps<mpfr::mpreal>(mpfr::mpreal) noexcept
    { return std::numeric_limits<mpfr::mpreal>::traps; }

  template<>
    inline bool
    tinyness_before<mpfr::mpreal>(mpfr::mpreal) noexcept
    { return std::numeric_limits<mpfr::mpreal>::tinyness_before; }

  template<>
    inline std::float_round_style
    round_style<mpfr::mpreal>(mpfr::mpreal) noexcept
    { return std::numeric_limits<mpfr::mpreal>::round_style(); }

  // Extra bits to help with numerics...
  // These depend math functions which aren't constexpr for mpfr::mpreal.
  // These are specializations of the functions in bits/numeric_limits.h

  template<>
    inline mpfr::mpreal
    max_integer(mpfr::mpreal proto) noexcept
    { return std::ldexp(1, digits(proto)); }

  template<>
    inline mpfr::mpreal
    sqrt_max<mpfr::mpreal>(mpfr::mpreal proto) noexcept
    { return std::sqrt(lim_max(proto)); }

  template<>
    inline mpfr::mpreal
    cbrt_max<mpfr::mpreal>(mpfr::mpreal proto) noexcept
    { return std::cbrt(lim_max(proto)); }

  template<>
    inline mpfr::mpreal
    root_max(mpfr::mpreal root) noexcept
    { return std::pow(lim_max(root), 1 / root); }

  template<>
    inline mpfr::mpreal
    log_max<mpfr::mpreal>(mpfr::mpreal proto) noexcept
    { return std::log(lim_max(proto)); }

  template<>
    inline mpfr::mpreal
    log10_max<mpfr::mpreal>(mpfr::mpreal proto) noexcept
    { return std::log10(lim_max(proto)); }

  template<>
    inline mpfr::mpreal
    sqrt_min<mpfr::mpreal>(mpfr::mpreal proto) noexcept
    { return std::sqrt(lim_min(proto)); }

  template<>
    inline mpfr::mpreal
    cbrt_min<mpfr::mpreal>(mpfr::mpreal proto) noexcept
    { return std::cbrt(lim_min(proto)); }

  template<>
    inline mpfr::mpreal
    root_min(mpfr::mpreal root) noexcept
    { return std::pow(lim_min(root), 1 / root); }

  template<>
    inline mpfr::mpreal
    log_min<mpfr::mpreal>(mpfr::mpreal proto) noexcept
    { return std::log(lim_min(proto)); }

  template<>
    inline mpfr::mpreal
    log10_min<mpfr::mpreal>(mpfr::mpreal proto) noexcept
    { return std::log10(lim_min(proto)); }

  template<>
    inline mpfr::mpreal
    sqrt_eps<mpfr::mpreal>(mpfr::mpreal proto) noexcept
    { return std::sqrt(epsilon(proto)); }

  template<>
    inline mpfr::mpreal
    cbrt_eps<mpfr::mpreal>(mpfr::mpreal proto) noexcept
    { return std::cbrt(epsilon(proto)); }

  template<>
    inline mpfr::mpreal
    root_eps(mpfr::mpreal root) noexcept
    { return std::pow(epsilon(root), 1 / root); }

  template<>
    inline mpfr::mpreal
    log_eps<mpfr::mpreal>(mpfr::mpreal proto) noexcept
    { return std::log(epsilon(proto)); }

  template<>
    inline mpfr::mpreal
    log10_eps<mpfr::mpreal>(mpfr::mpreal proto) noexcept
    { return std::log10(epsilon(proto)); }

} // namespace emsr

#endif // NUMERIC_LIMITS_MPREAL_H
