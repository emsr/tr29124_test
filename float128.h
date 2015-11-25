#ifndef FLOAT128_H
#define FLOAT128_H 1

#if !defined(__STRICT_ANSI__) && defined(_GLIBCXX_USE_FLOAT128)

#include <limits>
#include <iostream>
#include <iomanip> // For setw().
#include <sstream>
#include <quadmath.h>

// From <limits>
#define __glibcxx_max_digits10(T) \
  (2 + (T) * 643L / 2136)

namespace std
{

  inline __float128
  abs(__float128 __x) _GLIBCXX_USE_NOEXCEPT
  { return fabsq(__x); }

  inline __float128
  acos(__float128 __x) _GLIBCXX_USE_NOEXCEPT
  { return acosq(__x); }

  inline __float128
  asin(__float128 __x) _GLIBCXX_USE_NOEXCEPT
  { return asinq(__x); }

  inline __float128
  atan(__float128 __x) _GLIBCXX_USE_NOEXCEPT
  { return atanq(__x); }

  inline __float128
  atan2(__float128 __y, __float128 __x) _GLIBCXX_USE_NOEXCEPT
  { return atan2q(__y, __x); }

  inline __float128
  cbrt(__float128 __x) _GLIBCXX_USE_NOEXCEPT
  { return cbrtq(__x); }

  inline __float128
  ceil(__float128 __x) _GLIBCXX_USE_NOEXCEPT
  { return ceilq(__x); }

  inline __float128
  copysign(__float128 __x, __float128 __y) _GLIBCXX_USE_NOEXCEPT
  { return copysignq(__x, __y); }

  inline __float128
  cos(__float128 __x) _GLIBCXX_USE_NOEXCEPT
  { return cosq(__x); }

  inline __float128
  cosh(__float128 __x) _GLIBCXX_USE_NOEXCEPT
  { return coshq(__x); }

  inline __float128
  exp(__float128 __x) _GLIBCXX_USE_NOEXCEPT
  { return expq(__x); }

  inline __float128
  erf(__float128 __x) _GLIBCXX_USE_NOEXCEPT
  { return erfq(__x); }

  inline __float128
  erfc(__float128 __x) _GLIBCXX_USE_NOEXCEPT
  { return erfcq(__x); }

  inline __float128
  expm1(__float128 __x) _GLIBCXX_USE_NOEXCEPT
  { return expm1q(__x); }

  inline __float128
  fabs(__float128 __x) _GLIBCXX_USE_NOEXCEPT
  { return fabsq(__x); }

  inline __float128
  fdim(__float128 __x, __float128 __y) _GLIBCXX_USE_NOEXCEPT
  { return fdimq(__x, __y); }

  inline __float128
  floor(__float128 __x) _GLIBCXX_USE_NOEXCEPT
  { return floorq(__x); }

  inline __float128
  fma(__float128 __m, __float128 __x, __float128 __b) _GLIBCXX_USE_NOEXCEPT
  { return fmaq(__m, __x, __b); }

  inline __float128
  fmax(__float128 __x, __float128 __y) _GLIBCXX_USE_NOEXCEPT
  { return fmaxq(__x, __y); }

  inline __float128
  fmin(__float128 __x, __float128 __y) _GLIBCXX_USE_NOEXCEPT
  { return fminq(__x, __y); }

  inline __float128
  fmod(__float128 __x, __float128 __y) _GLIBCXX_USE_NOEXCEPT
  { return fmodq(__x, __y); }

  inline __float128
  frexp(__float128 __x, int* __exp) _GLIBCXX_USE_NOEXCEPT
  { return frexpq(__x, __exp); }

  inline __float128
  hypot(__float128 __x, __float128 __y) _GLIBCXX_USE_NOEXCEPT
  { return hypotq(__x, __y); }

  inline int
  isinf(__float128 __x) _GLIBCXX_USE_NOEXCEPT
  { return isinfq(__x); }

  inline int
  ilogb(__float128 __x) _GLIBCXX_USE_NOEXCEPT
  { return ilogbq(__x); }

  inline int
  isnan(__float128 __x) _GLIBCXX_USE_NOEXCEPT
  { return isnanq(__x); }

  inline __float128
  j0(__float128 __x) _GLIBCXX_USE_NOEXCEPT
  { return j0q(__x); }

  inline __float128
  j1(__float128 __x) _GLIBCXX_USE_NOEXCEPT
  { return j1q(__x); }

  inline __float128
  jn(int __n, __float128 __x) _GLIBCXX_USE_NOEXCEPT
  { return jnq(__n, __x); }

  inline __float128
  ldexp(__float128 __x, int __exp) _GLIBCXX_USE_NOEXCEPT
  { return ldexpq(__x, __exp); }

  inline __float128
  lgamma(__float128 __x) _GLIBCXX_USE_NOEXCEPT
  { return lgammaq(__x); }

  inline long long int
  llrint(__float128 __x) _GLIBCXX_USE_NOEXCEPT
  { return llrintq(__x); }

  inline long long int
  llround(__float128 __x) _GLIBCXX_USE_NOEXCEPT
  { return llroundq(__x); }

  inline __float128
  logb(__float128 __x) _GLIBCXX_USE_NOEXCEPT
  { return logbq(__x); }

  inline __float128
  log(__float128 __x) _GLIBCXX_USE_NOEXCEPT
  { return logq(__x); }

  inline __float128
  log10(__float128 __x) _GLIBCXX_USE_NOEXCEPT
  { return log10q(__x); }

  inline __float128
  log2(__float128 __x) _GLIBCXX_USE_NOEXCEPT
  { return log2q(__x); }

  inline __float128
  log1p(__float128 __x) _GLIBCXX_USE_NOEXCEPT
  { return log1pq(__x); }

  inline long int
  lrint(__float128 __x) _GLIBCXX_USE_NOEXCEPT
  { return lrintq(__x); }

  inline long int
  lround(__float128 __x) _GLIBCXX_USE_NOEXCEPT
  { return lroundq(__x); }

  inline __float128
  modf(__float128 __x, __float128* __iptr) _GLIBCXX_USE_NOEXCEPT
  { return modfq(__x, __iptr); }

  inline __float128
  nan(const char * __msg) _GLIBCXX_USE_NOEXCEPT
  { return nanq(__msg); }

  inline __float128
  nearbyint(__float128 __x) _GLIBCXX_USE_NOEXCEPT
  { return nearbyintq(__x); }

  inline __float128
  nextafter(__float128 __x, __float128 __y) _GLIBCXX_USE_NOEXCEPT
  { return nextafterq(__x, __y); }

  //inline __float128
  //powi(__float128 __x, int __n) _GLIBCXX_USE_NOEXCEPT
  //{ return powiq(__x, __n); }

  inline __float128
  remainder(__float128 __x, __float128 __y) _GLIBCXX_USE_NOEXCEPT
  { return remainderq(__x, __y); }

  inline __float128
  remquo(__float128 __x, __float128 __y, int* __n) _GLIBCXX_USE_NOEXCEPT
  { return remquoq(__x, __y, __n); }

  inline __float128
  rint(__float128 __x) _GLIBCXX_USE_NOEXCEPT
  { return rintq(__x); }

  inline __float128
  round(__float128 __x) _GLIBCXX_USE_NOEXCEPT
  { return roundq(__x); }

  inline __float128
  scalbln(__float128 __x, long int __n) _GLIBCXX_USE_NOEXCEPT
  { return scalblnq(__x, __n); }

  inline __float128
  scalbn(__float128 __x, int __n) _GLIBCXX_USE_NOEXCEPT
  { return scalbnq(__x, __n); }

  inline int
  signbitq(__float128 __x) _GLIBCXX_USE_NOEXCEPT
  { return signbitq(__x); }

  inline void
  sincos(__float128 __x, __float128 * __sin, __float128 * __cos)
  _GLIBCXX_USE_NOEXCEPT
  { return sincosq(__x, __sin, __cos); }

  inline __float128
  sin(__float128 __x) _GLIBCXX_USE_NOEXCEPT
  { return sinq(__x); }

  inline __float128
  sinh(__float128 __x) _GLIBCXX_USE_NOEXCEPT
  { return sinhq(__x); }

  inline __float128
  sqrt(__float128 __x) _GLIBCXX_USE_NOEXCEPT
  { return sqrtq(__x); }

  inline __float128
  tan(__float128 __x) _GLIBCXX_USE_NOEXCEPT
  { return tanq(__x); }

  inline __float128
  tanh(__float128 __x) _GLIBCXX_USE_NOEXCEPT
  { return tanhq(__x); }

  inline __float128
  tgamma(__float128 __x) _GLIBCXX_USE_NOEXCEPT
  { return tgammaq(__x); }

  inline __float128
  trunc(__float128 __x) _GLIBCXX_USE_NOEXCEPT
  { return truncq(__x); }

  inline __float128
  y0(__float128 __x) _GLIBCXX_USE_NOEXCEPT
  { return y0q(__x); }

  inline __float128
  y1(__float128 __x) _GLIBCXX_USE_NOEXCEPT
  { return y1q(__x); }

  inline __float128
  yn(int __n, __float128 __x) _GLIBCXX_USE_NOEXCEPT
  { return ynq(__n, __x); }

  inline std::ostream&
  operator<<(std::ostream& __os, const __float128& __x)
  {
    std::ostringstream __fmt;
    __fmt << "%." << __os.precision() << "Qg";
    constexpr int __strlen = 1000;
    char __str[__strlen];
    quadmath_snprintf(__str, __strlen, __fmt.str().c_str(), __x) ;
    __os.precision(30);
    __os << __str;
    return __os;
  }

  inline std::istream&
  operator>>(std::istream& __is, __float128& __x)
  {
    constexpr int __strlen = 1000;
    char __str[__strlen];
    __is >> std::setw(__strlen) >> __str;
    __x = strtoflt128(__str, 0);
    return __is;
  }

  /// numeric_limits<__float128> specialization.
  template<>
    struct numeric_limits<__float128>
    {
      static _GLIBCXX_USE_CONSTEXPR bool is_specialized = true;

      static _GLIBCXX_CONSTEXPR __float128 
      min() _GLIBCXX_USE_NOEXCEPT { return FLT128_MIN; }

      static _GLIBCXX_CONSTEXPR __float128 
      max() _GLIBCXX_USE_NOEXCEPT { return FLT128_MAX; }

#if __cplusplus >= 201103L
      static constexpr __float128 
      lowest() _GLIBCXX_USE_NOEXCEPT { return -FLT128_MAX; }
#endif

      static _GLIBCXX_USE_CONSTEXPR int digits = FLT128_MANT_DIG;
      static _GLIBCXX_USE_CONSTEXPR int digits10 = FLT128_DIG;
#if __cplusplus >= 201103L
      static _GLIBCXX_USE_CONSTEXPR int max_digits10
	 = __glibcxx_max_digits10 (FLT128_MANT_DIG);
#endif
      static _GLIBCXX_USE_CONSTEXPR bool is_signed = true;
      static _GLIBCXX_USE_CONSTEXPR bool is_integer = false;
      static _GLIBCXX_USE_CONSTEXPR bool is_exact = false;
      static _GLIBCXX_USE_CONSTEXPR int radix = __FLT_RADIX__;

      static _GLIBCXX_CONSTEXPR __float128 
      epsilon() _GLIBCXX_USE_NOEXCEPT { return FLT128_EPSILON; }

      static _GLIBCXX_CONSTEXPR __float128 
      round_error() _GLIBCXX_USE_NOEXCEPT { return 0.5Q; }

      static _GLIBCXX_USE_CONSTEXPR int min_exponent = FLT128_MIN_EXP;
      static _GLIBCXX_USE_CONSTEXPR int min_exponent10 = FLT128_MIN_10_EXP;
      static _GLIBCXX_USE_CONSTEXPR int max_exponent = FLT128_MAX_EXP;
      static _GLIBCXX_USE_CONSTEXPR int max_exponent10 = FLT128_MAX_10_EXP;

      static _GLIBCXX_USE_CONSTEXPR bool has_infinity = true;
      static _GLIBCXX_USE_CONSTEXPR bool has_quiet_NaN = true;
      static _GLIBCXX_USE_CONSTEXPR bool has_signaling_NaN = has_quiet_NaN;
      static _GLIBCXX_USE_CONSTEXPR float_denorm_style has_denorm
	= denorm_present;
      static _GLIBCXX_USE_CONSTEXPR bool has_denorm_loss = true;

      static _GLIBCXX_CONSTEXPR __float128
      infinity() _GLIBCXX_USE_NOEXCEPT { return __builtin_huge_valq(); }

      static /*_GLIBCXX_CONSTEXPR*/ __float128 
      quiet_NaN() _GLIBCXX_USE_NOEXCEPT { return nanq(""); }

      //static _GLIBCXX_CONSTEXPR __float128 
      //signaling_NaN() _GLIBCXX_USE_NOEXCEPT { return __builtin_nansq(""); }

      static _GLIBCXX_CONSTEXPR __float128 
      denorm_min() _GLIBCXX_USE_NOEXCEPT { return FLT128_DENORM_MIN; }

      static _GLIBCXX_USE_CONSTEXPR bool is_iec559
	= has_infinity && has_quiet_NaN && has_denorm == denorm_present;
      static _GLIBCXX_USE_CONSTEXPR bool is_bounded = true;
      static _GLIBCXX_USE_CONSTEXPR bool is_modulo = false;

      static _GLIBCXX_USE_CONSTEXPR bool traps = false;//???
      static _GLIBCXX_USE_CONSTEXPR bool tinyness_before = 
					 false;//???
      static _GLIBCXX_USE_CONSTEXPR float_round_style round_style = 
						      round_to_nearest;
    };

} // namespace std

// From <limits>
#undef __glibcxx_max_digits10

#endif // __STRICT_ANSI__ && _GLIBCXX_USE_FLOAT128

#endif // FLOAT128_H
