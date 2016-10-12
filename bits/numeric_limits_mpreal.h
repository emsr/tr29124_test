
#ifndef _GLIBCXX_BITS_NUMERIC_LIMITS_MPREAL_H
#define _GLIBCXX_BITS_NUMERIC_LIMITS_MPREAL_H 1

#include <limits>
#include <mpreal.h>
#include <bits/math_mpreal.h>

namespace __gnu_cxx _GLIBCXX_VISIBILITY(default)
{
_GLIBCXX_BEGIN_NAMESPACE_VERSION

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
    bool
    __is_specialized<mpfr::mpreal>(mpfr::mpreal) _GLIBCXX_USE_NOEXCEPT
    { return std::numeric_limits<mpfr::mpreal>::is_specialized; }

  template<>
    mpfr::mpreal
    __min<mpfr::mpreal>(mpfr::mpreal __proto) _GLIBCXX_USE_NOEXCEPT
    { return mpfr::minval(__proto.getPrecision()); }

  template<>
    mpfr::mpreal
    __max<mpfr::mpreal>(mpfr::mpreal __proto) _GLIBCXX_USE_NOEXCEPT
    { return mpfr::maxval(__proto.getPrecision()); }

  template<>
    mpfr::mpreal
    __lowest<mpfr::mpreal>(mpfr::mpreal __proto) _GLIBCXX_USE_NOEXCEPT
    { return -mpfr::maxval(__proto.getPrecision()); }

  template<>
    int
    __digits<mpfr::mpreal>(mpfr::mpreal __proto) _GLIBCXX_USE_NOEXCEPT
    { return __proto.getPrecision(); }

  template<>
    int
    __digits10<mpfr::mpreal>(mpfr::mpreal __proto) _GLIBCXX_USE_NOEXCEPT
    { return mpfr::bits2digits(__proto.getPrecision()); }

  template<>
    int
    __max_digits10<mpfr::mpreal>(mpfr::mpreal __proto) _GLIBCXX_USE_NOEXCEPT
    { return __digits10(__proto); }

  template<>
    bool
    __is_signed<mpfr::mpreal>(mpfr::mpreal) _GLIBCXX_USE_NOEXCEPT
    { return std::numeric_limits<mpfr::mpreal>::is_signed; }

  template<>
    bool
    __is_integer<mpfr::mpreal>(mpfr::mpreal) _GLIBCXX_USE_NOEXCEPT
    { return std::numeric_limits<mpfr::mpreal>::is_integer; }

  template<>
    bool
    __is_exact<mpfr::mpreal>(mpfr::mpreal) _GLIBCXX_USE_NOEXCEPT
    { return std::numeric_limits<mpfr::mpreal>::is_exact; }

  template<>
    int
    __radix<mpfr::mpreal>(mpfr::mpreal) _GLIBCXX_USE_NOEXCEPT
    { return std::numeric_limits<mpfr::mpreal>::radix; }

  template<>
    mpfr::mpreal
    __epsilon<mpfr::mpreal>(mpfr::mpreal __proto) _GLIBCXX_USE_NOEXCEPT
    { return  mpfr::machine_epsilon(__proto.getPrecision()); }

  template<>
    mpfr::mpreal
    __round_error<mpfr::mpreal>(mpfr::mpreal __proto) _GLIBCXX_USE_NOEXCEPT
    { return std::numeric_limits<mpfr::mpreal>::round_error(__proto.getPrecision()); }

  template<>
    int
    __min_exponent<mpfr::mpreal>(mpfr::mpreal) _GLIBCXX_USE_NOEXCEPT
    { return std::numeric_limits<mpfr::mpreal>::min_exponent; }

  template<>
    int
    __min_exponent10<mpfr::mpreal>(mpfr::mpreal) _GLIBCXX_USE_NOEXCEPT
    { return std::numeric_limits<mpfr::mpreal>::min_exponent10; }

  template<>
    int
    __max_exponent<mpfr::mpreal>(mpfr::mpreal) _GLIBCXX_USE_NOEXCEPT
    { return std::numeric_limits<mpfr::mpreal>::max_exponent; }

  template<>
    int
    __max_exponent10<mpfr::mpreal>(mpfr::mpreal) _GLIBCXX_USE_NOEXCEPT
    { return std::numeric_limits<mpfr::mpreal>::max_exponent10; }

  template<>
    bool
    __has_infinity<mpfr::mpreal>(mpfr::mpreal) _GLIBCXX_USE_NOEXCEPT
    { return std::numeric_limits<mpfr::mpreal>::has_infinity; }

  template<>
    bool
    __has_quiet_NaN<mpfr::mpreal>(mpfr::mpreal) _GLIBCXX_USE_NOEXCEPT
    { return std::numeric_limits<mpfr::mpreal>::has_quiet_NaN; }

  template<>
    bool
    __has_signaling_NaN<mpfr::mpreal>(mpfr::mpreal) _GLIBCXX_USE_NOEXCEPT
    { return std::numeric_limits<mpfr::mpreal>::has_signaling_NaN; }

  template<>
    std::float_denorm_style
    __has_denorm<mpfr::mpreal>(mpfr::mpreal) _GLIBCXX_USE_NOEXCEPT
    { return std::numeric_limits<mpfr::mpreal>::has_denorm; }

  // No std::numeric_limits<mpfr::mpreal>::has_denorm_loss.
  template<>
    bool
    __has_denorm_loss<mpfr::mpreal>(mpfr::mpreal) _GLIBCXX_USE_NOEXCEPT
    { return std::numeric_limits<mpfr::mpreal>::has_denorm; }

  template<>
    mpfr::mpreal
    __infinity<mpfr::mpreal>(mpfr::mpreal) _GLIBCXX_USE_NOEXCEPT
    { return std::numeric_limits<mpfr::mpreal>::infinity(); }

  template<>
    mpfr::mpreal
    __quiet_NaN<mpfr::mpreal>(mpfr::mpreal) _GLIBCXX_USE_NOEXCEPT
    { return std::numeric_limits<mpfr::mpreal>::quiet_NaN(); }

  template<>
    mpfr::mpreal
    __signaling_NaN<mpfr::mpreal>(mpfr::mpreal) _GLIBCXX_USE_NOEXCEPT
    { return std::numeric_limits<mpfr::mpreal>::signaling_NaN(); }

  template<>
    mpfr::mpreal
    __denorm_min<mpfr::mpreal>(mpfr::mpreal) _GLIBCXX_USE_NOEXCEPT
    { return std::numeric_limits<mpfr::mpreal>::denorm_min(); }

  template<>
    bool
    __is_iec559<mpfr::mpreal>(mpfr::mpreal) _GLIBCXX_USE_NOEXCEPT
    { return std::numeric_limits<mpfr::mpreal>::is_iec559; }

  template<>
    bool
    __is_bounded<mpfr::mpreal>(mpfr::mpreal) _GLIBCXX_USE_NOEXCEPT
    { return std::numeric_limits<mpfr::mpreal>::is_bounded; }

  template<>
    bool
    __is_modulo<mpfr::mpreal>(mpfr::mpreal) _GLIBCXX_USE_NOEXCEPT
    { return std::numeric_limits<mpfr::mpreal>::is_modulo; }

  template<>
    bool
    __traps<mpfr::mpreal>(mpfr::mpreal) _GLIBCXX_USE_NOEXCEPT
    { return std::numeric_limits<mpfr::mpreal>::traps; }

  template<>
    bool
    __tinyness_before<mpfr::mpreal>(mpfr::mpreal) _GLIBCXX_USE_NOEXCEPT
    { return std::numeric_limits<mpfr::mpreal>::tinyness_before; }

  template<>
    std::float_round_style
    __round_style<mpfr::mpreal>(mpfr::mpreal) _GLIBCXX_USE_NOEXCEPT
    { return std::numeric_limits<mpfr::mpreal>::round_style(); }

  // Extra bits to help with numerics...
  // These depend math functions which aren't constexpr for mpfr::mpreal.
  // These are specializations of the functions in bits/numeric_limits.h

  template<>
    mpfr::mpreal
    __sqrt_max<mpfr::mpreal>(mpfr::mpreal __proto) _GLIBCXX_USE_NOEXCEPT
    { return std::sqrt(__max(__proto)); }

#ifdef NO_CBRT
  template<>
    mpfr::mpreal
    __cbrt_max<mpfr::mpreal>(mpfr::mpreal __proto) _GLIBCXX_USE_NOEXCEPT
    { return std::pow(__max(__proto), 1 / mpfr::mpreal{3}); }
#else
  template<>
    mpfr::mpreal
    __cbrt_max<mpfr::mpreal>(mpfr::mpreal __proto) _GLIBCXX_USE_NOEXCEPT
    { return std::cbrt(__max(__proto)); }
#endif

  template<>
    mpfr::mpreal
    __root_max(mpfr::mpreal __root) _GLIBCXX_USE_NOEXCEPT
    { return std::pow(__max(__root), 1 / __root); }

  template<>
    mpfr::mpreal
    __log_max<mpfr::mpreal>(mpfr::mpreal __proto) _GLIBCXX_USE_NOEXCEPT
    { return std::log(__max(__proto)); }

  template<>
    mpfr::mpreal
    __log10_max<mpfr::mpreal>(mpfr::mpreal __proto) _GLIBCXX_USE_NOEXCEPT
    { return std::log10(__max(__proto)); }


  template<>
    mpfr::mpreal
    __sqrt_min<mpfr::mpreal>(mpfr::mpreal __proto) _GLIBCXX_USE_NOEXCEPT
    { return std::sqrt(__min(__proto)); }

#ifdef NO_CBRT
  template<>
    mpfr::mpreal
    __cbrt_max<mpfr::mpreal>(mpfr::mpreal __proto) _GLIBCXX_USE_NOEXCEPT
    { return std::pow(__min(__proto), 1 / mpfr::mpreal{3}); }
#else
  template<>
    mpfr::mpreal
    __cbrt_min<mpfr::mpreal>(mpfr::mpreal __proto) _GLIBCXX_USE_NOEXCEPT
    { return std::cbrt(__min(__proto)); }
#endif

  template<>
    mpfr::mpreal
    __root_min(mpfr::mpreal __root) _GLIBCXX_USE_NOEXCEPT
    { return std::pow(__min(__root), 1 / __root); }

  template<>
    mpfr::mpreal
    __log_min<mpfr::mpreal>(mpfr::mpreal __proto) _GLIBCXX_USE_NOEXCEPT
    { return std::log(__min(__proto)); }

  template<>
    mpfr::mpreal
    __log10_min<mpfr::mpreal>(mpfr::mpreal __proto) _GLIBCXX_USE_NOEXCEPT
    { return std::log10(__min(__proto)); }

  template<>
    mpfr::mpreal
    __sqrt_eps<mpfr::mpreal>(mpfr::mpreal __proto) _GLIBCXX_USE_NOEXCEPT
    { return std::sqrt(__epsilon(__proto)); }

#ifdef NO_CBRT
  template<>
    mpfr::mpreal
    __cbrt_max<mpfr::mpreal>(mpfr::mpreal __proto) _GLIBCXX_USE_NOEXCEPT
    { return std::pow(__epsilon(__proto), 1 / mpfr::mpreal{3}); }
#else
  template<>
    mpfr::mpreal
    __cbrt_eps<mpfr::mpreal>(mpfr::mpreal __proto) _GLIBCXX_USE_NOEXCEPT
    { return std::cbrt(__epsilon(__proto)); }
#endif

  template<>
    mpfr::mpreal
    __root_eps(mpfr::mpreal __root) _GLIBCXX_USE_NOEXCEPT
    { return std::pow(__epsilon(__root), 1 / __root); }

  template<>
    mpfr::mpreal
    __log_eps<mpfr::mpreal>(mpfr::mpreal __proto) _GLIBCXX_USE_NOEXCEPT
    { return std::log(__epsilon(__proto)); }

  template<>
    mpfr::mpreal
    __log10_eps<mpfr::mpreal>(mpfr::mpreal __proto) _GLIBCXX_USE_NOEXCEPT
    { return std::log10(__epsilon(__proto)); }

_GLIBCXX_END_NAMESPACE_VERSION
} // namespace __gnu_cxx

#endif // _GLIBCXX_BITS_NUMERIC_LIMITS_MPREAL_H
