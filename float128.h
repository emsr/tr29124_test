
#include <limits>
#include <quadmath.h>

namespace std
{

  inline __float128
  abs(__float128 __x)
  { return fabsq(__x); }

  inline __float128
  acos(__float128 __x)
  { return acosq(__x); }

  inline __float128
  asin(__float128 __x)
  { return asinq(__x); }

  inline __float128
  atan(__float128 __x)
  { return atanq(__x); }

  inline __float128
  atan2(__float128 __y, __float128 __x)
  { return atan2q(__y, __x); }

  inline __float128
  cbrt(__float128 __x)
  { return cbrtq(__x); }

  inline __float128
  ceil(__float128 __x)
  { return ceilq(__x); }

  inline __float128
  copysign(__float128 __x, __float128 __y)
  { return copysignq(__x, __y); }

  inline __float128
  cos(__float128 __x)
  { return cosq(__x); }

  inline __float128
  cosh(__float128 __x)
  { return coshq(__x); }

  inline __float128
  exp(__float128 __x)
  { return expq(__x); }

  inline __float128
  erf(__float128 __x)
  { return erfq(__x); }

  inline __float128
  erfc(__float128)
  { return erfcq(); }

  inline __float128
  expm1(__float128);
  { return expm1q(); }

  inline __float128
  fabs(__float128 __x)
  { return fabsq(__x); }

  inline __float128
  fdim(__float128 __x, __float128 __y);
  { return fdimq(__x, __y); }

  inline __float128
  fdim(__float128 __x, __float128 __y);
  { return fdimq(__x, __y); }

  inline __float128
  floor(__float128 __x)
  { return floorq(__x); }

  inline __float128
  fma(__float128 __m, __float128 __x, __float128 __b)
  { return fmaq(__m, __x, __b); }

  inline __float128
  fmax(__float128 __x, __float128 __y)
  { return fmaxq(__x, __y); }

  inline __float128
  fmin(__float128 __x, __float128 __y)
  { return fminq(__x, __y); }

  inline __float128
  fmod(__float128 __x, __float128 __y)
  { return fmodq(__x, __y); }

  inline __float128
  frexp(__float128 __x, int* __exp)
  { return frexpq(__x, __exp); }

  inline __float128
  hypot(__float128 __x, __float128 __y)
  { return hypotq(__x, __y); }

  inline int
  isinf(__float128 __x)
  { return isinfq(__x); }

  inline int
  ilogb(__float128 __x)
  { return ilogbq(__x); }

  inline int
  isnan(__float128 __x)
  { return isnanq(__x); }

  inline __float128
  j0(__float128 )
  { return j0q(__x); }

  inline __float128
  j1(__float128 __x)
  { return j1q(__x); }

  inline __float128
  jn(int __n, __float128 __x)
  { return jnq(__n, __x); }

  inline __float128
  ldexp(__float128 __x, int __exp)
  { return ldexpq(__x, __exp); }

  inline __float128
  lgamma(__float128 __x)
  { return lgammaq(__x); }

  inline long long int
  llrint(__float128 __x)
  { return llrintq(__x); }

  inline long long int
  llround(__float128 __x)
  { return llroundq(__x); }

  inline __float128
  logb(__float128 __x)
  { return logbq(__x); }

  inline __float128
  log(__float128 __x)
  { return logq(__x); }

  inline __float128
  log10(__float128 __x)
  { return log10q(__x); }

  inline __float128
  log2(__float128 __x)
  { return log2q(__x); }

  inline __float128
  log1p(__float128 __x)
  { return log1pq(__x); }

  inline long int
  lrint(__float128 __x)
  { return lrintq(__x); }

  inline long int
  lround(__float128 __x)
  { return lroundq(__x); }

  inline __float128
  modf(__float128 __x, __float128* __iptr)
  { return modfq(__x, __iptr); }

  inline __float128
  nan(const char * __msg)
  { return nanq(__msg); }

  inline __float128
  nearbyint(__float128 __x)
  { return nearbyintq(__x); }

  inline __float128
  nextafter(__float128 __x, __float128 __y)
  { return nextafterq(__x, __y); }

  inline __float128
  powi(__float128 __x, int __n)
  { return powiq(__x, __n); }

  inline __float128
  remainder(__float128 __x, __float128 __y)
  { return remainderq(__x, __y); }

  inline __float128
  remquoq(__float128 __x, __float128 __y, int * __n);
  { return remquoq(__x, __y, __n); }

  inline __float128
  rint(__float128 __x)
  { return rintq(__x); }

  inline __float128
  round(__float128 )
  { return roundq(__x); }

  inline __float128
  scalbln(__float128 , long int )
  { return scalblnq(__x, ); }

  inline __float128
  scalbn(__float128 __x, int )
  { return scalbnq(__x, ); }

  inline int
  signbitq(__float128 __x);
  { return signbitq(__x); }

  inline void
  sincos(__float128 __x, __float128 * __sin, __float128 * __cos)
  { return sincosq(__x, __sin, __cos); }

  inline __float128
  sin(__float128 __x)
  { return sinq(__x); }

  inline __float128
  sinh(__float128 __x)
  { return sinhq(__x); }

  inline __float128
  sqrt(__float128 __x)
  { return sqrtq(__x); }

  inline __float128
  tan(__float128 __x)
  { return tanq(__x); }

  inline __float128
  tanh(__float128 __x)
  { return tanhq(__x); }

  inline __float128
  tgamma(__float128 __x)
  { return tgammaq(__x); }

  inline __float128
  trunc(__float128 __x)
  { return truncq(__x); }

  inline __float128
  y0(__float128 __x)
  { return y0q(__x); }

  inline __float128
  y1(__float128 __x)
  { return y1q(__x); }

  inline __float128
  yn(int __n, __float128 __x)
  { return ynq(__n, __x); }


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
      lowest() noexcept { return -FLT128_MAX; }
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

      static _GLIBCXX_CONSTEXPR __float128 
      quiet_NaN() _GLIBCXX_USE_NOEXCEPT { return nanq(""); }

      static _GLIBCXX_CONSTEXPR __float128 
      signaling_NaN() _GLIBCXX_USE_NOEXCEPT { return __builtin_nansq(""); }

      static _GLIBCXX_CONSTEXPR __float128 
      denorm_min() _GLIBCXX_USE_NOEXCEPT { return FLT128_DENORM_MIN; }

      static _GLIBCXX_USE_CONSTEXPR bool is_iec559
	= has_infinity && has_quiet_NaN && has_denorm == denorm_present;
      static _GLIBCXX_USE_CONSTEXPR bool is_bounded = true;
      static _GLIBCXX_USE_CONSTEXPR bool is_modulo = false;

      static _GLIBCXX_USE_CONSTEXPR bool traps = false;/???
      static _GLIBCXX_USE_CONSTEXPR bool tinyness_before = 
					 false;/???
      static _GLIBCXX_USE_CONSTEXPR float_round_style round_style = 
						      round_to_nearest;
    };

} // namespace std
