// Special functions -*- C++ -*-

// Copyright (C) 2016-2019 Free Software Foundation, Inc.
//
// This file is part of the GNU ISO C++ Library.  This library is free
// software; you can redistribute it and/or modify it under the
// terms of the GNU General Public License as published by the
// Free Software Foundation; either version 3, or (at your option)
// any later version.

// This library is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.

// Under Section 7 of GPL version 3, you are granted additional
// permissions described in the GCC Runtime Library Exception, version
// 3.1, as published by the Free Software Foundation.

// You should have received a copy of the GNU General Public License and
// a copy of the GCC Runtime Library Exception along with this program;
// see the files COPYING3 and COPYING.RUNTIME respectively.  If not, see
// <http://www.gnu.org/licenses/>.

#ifndef _GLIBCXX_BITS_NOTSOSPECFUN_H
#define _GLIBCXX_BITS_NOTSOSPECFUN_H 1

#pragma GCC system_header

#include <variant>
#include <complex>

namespace std _GLIBCXX_VISIBILITY(default)
{
_GLIBCXX_BEGIN_NAMESPACE_VERSION

  // Use narrow structs for aggregate return types.
  // Prefer returns to pointer or ref arguments.
  // This will work swimmingly with structured bindings.
  //
  // I think all the basic math functions should be constexpr.
  // Error haqndling - the old ones have global error errno (from C).
  // The specfuns throw.
  //   - should these and the old basic ones have throwing versions
  //     - since we don't like global error reporting
  //     - exception is part of the signature(?)
  //     - only if people can flip back to another or no error reporting...
  // We could do like filesystem and some others and double the api with
  // error_code return args.  In math, errors really are exceptional where
  // libs with this return failure is quite often an option.  In math failure
  // needs to be figured out, fixed, cleaned up after.

  // The more I think of it, the more I think that
  //  template<typename Real>
  //    using numeric_t = std::variant<Real, std::complex<Real>>;
  // is a thing and is the answer to lgamma, negative arg bessels,
  // polynomial roots, etc. return types.
  // The ship has sailed for lgamma (I think?)


  // Several functions have pointer arguments and so can't be constexpr.
  // Return types.  See p0533.

  // Should we bother with C-style suffixed functions?
  // I don't think providing default conversions to numeric types is good enough
  // to avoid collision and surprise.
  // Namespace? Naming?

  // Implementation-wise, these could be wrappers of the regular functions.
  // The pointers would be hidden inside the wrapper.

  // These used to be *div_t but they conflicted with 
  // /usr/include/stdlib.h:70:5: note: previous declaration 'typedef struct ldiv_t ldiv_t'
  template<typename IntTp>
    struct int_quot_t
    {
      IntTp quot;
      IntTp rem;
    };

  using squot_t = int_quot_t<short int>;
  using quot_t = int_quot_t<int>;
  using lquot_t = int_quot_t<long int>;
  using llquot_t = int_quot_t<long long int>;
  using intmaxquot_t = int_quot_t<std::intmax_t>;

  constexpr squot_t quot(short int numer, short int denom);
  constexpr quot_t quot(int numer, int denom);
  constexpr lquot_t quot(long int numer, long int denom);
  constexpr llquot_t quot(long long int numer, long long int denom);
  constexpr intmaxquot_t quot(std::intmax_t x, std::intmax_t y);

  constexpr squot_t squot(short int numer, short int denom);
  constexpr lquot_t lquot(long int numer, long int denom);
  constexpr llquot_t llquot(long long int numer, long long int denom);
  constexpr intmaxquot_t intmaxquot(std::intmax_t x, std::intmax_t y);

  // Decompose floating-point number into a normalized fraction
  // and integral power of two.

  // frexp -> frac_exp2?
  template<typename Tp>
    struct frexp_t
    {
      Tp value;
      int exp2;
    };

  constexpr frexp_t<float> frexp(float value);
  constexpr frexp_t<double> frexp(double value);
  constexpr frexp_t<long double> frexp(long double value);
  constexpr frexp_t<float> frexpf(float value);
  constexpr frexp_t<long double> frexpl(long double value);

  // Decompose floating-point number into fractional and integer part.

  // modf -> fp_mod?
  template<typename Tp>
    struct modf_t
    {
      Tp frac_part;
      Tp int_part;
    };

  constexpr modf_t<float> modf(float value);
  constexpr modf_t<double> modf(double value);
  constexpr modf_t<long double> modf(long double value);
  constexpr modf_t<float> modff(float value);
  constexpr modf_t<long double> modfl(long double value);


  // Divide x by y providing remainder and integer quotient.
  // Should the int grow depending on Tp?
  // Could the int be unsigned?

  template<typename Tp>
    struct remquo_t
    {
      Tp remainder;
      int quotient;
    };

  constexpr remquo_t<float> remquo(float x, float y);
  constexpr remquo_t<double> remquo(double x, double y);
  constexpr remquo_t<long double> remquo(long double x, long double y);

  constexpr remquo_t<float> remquof(float x, float y);
  constexpr remquo_t<long double> remquol(long double x, long double y);

  constexpr double nan();
  constexpr float nanf();
  constexpr long double nanl();
  // And/or:
  template<char... Str>
    constexpr double nan();
  template<char... Str>
    constexpr float nanf();
  template<char... Str>
    constexpr long double nanl();


  // Log to arbitrary base - the inverse of pow(base, x).

  float logf(float base, float x);
  double log(double base, double x);
  long double logl(long double base, long double x);

  // Exponent base 10 - the inverse of log10.

  float exp10f(float x);
  double exp10log(double x);
  long double exp10logl(long double x);

  // Trigonometric functions

  // Combined sine and cosine.
  template<typename Tp>
    struct sincos_t
    {
      Tp sin_v;
      Tp cos_v;
    };

  template<typename Tp>
    Tp
    pi_v = static_cast<Tp>(3.1415926536897932384626L);

  sincos_t<float> sincosf(float x);
  sincos_t<double> sincos(double x);
  sincos_t<long double> sincosl(long double x);

  // Teach atan to use sincos_t.
  // This returns all four quadrants like atan2.
  float atanf(sincos_t<float> m);
  double atan(sincos_t<double> m);
  long double atanl(sincos_t<long double> m);

  // Auxilliary trigonometric functions...
  float secf(float x);
  double sec(double x);
  long double secl(long double x);

  float cscf(float x);
  double csc(double x);
  long double cscl(long double x);

  float cotf(float x);
  double cot(double x);
  long double cotl(long double x);

  // ... and their inverses.
  float asecf(float x);
  double asec(double x);
  long double asecl(long double x);

  float acscf(float x);
  double acsc(double x);
  long double acscl(long double x);

  float acotf(float x);
  double acot(double x);
  long double acotl(long double x);

  float acot2f(float x, float y);
  double acot2(double x, double y);
  long double acot2l(long double x, long double y);

  // This returns all four quadrants like acot2.
  float acotf(sincos_t<float> m);
  double acot(sincos_t<double> m);
  long double acotl(sincos_t<long double> m);


  // Hyperbolic functions

  // Combined sinh and cosh.
  template<typename Tp>
    struct sinhcosh_t
    {
      Tp sinh_value;
      Tp cosh_value;
    };

  // Teach atanh to use sinhcosh_t.
  float atanhf(sinhcosh_t<float> m);
  double atanh(sinhcosh_t<double> m);
  long double atanhl(sinhcosh_t<long double> m);

  // Auxilliary hyperbolic functions...
  sinhcosh_t<float> sinhcoshf(float x);
  sinhcosh_t<double> sinhcosh(double x);
  sinhcosh_t<long double> sinhcoshl(long double x);

  float sechf(float x);
  double sech(double x);
  long double sechl(long double x);

  float cschf(float x);
  double csch(double x);
  long double cschl(long double x);

  float cothf(float x);
  double coth(double x);
  long double cothl(long double x);

  // ... and their inverses.
  float asechf(float x);
  double asech(double x);
  long double asechl(long double x);

  float acschf(float x);
  double acsch(double x);
  long double acschl(long double x);

  float acothf(float x);
  double acoth(double x);
  long double acothl(long double x);

  float acothf(sinhcosh_t<float> m);
  double acoth(sinhcosh_t<double> m);
  long double acothl(sinhcosh_t<long double> m);


  // Reperiodized trigonometric functions...
  //   fun_pi(x) = fun(pi x);

  // This is really just another angle unit.
  // When we get units, this and deg, grad, rad would all get overloads.
  // We shouldn't need decorated functions.
  // OTOH, machines have these - there are built-ins and traditions...
  //
  // I want to have minimum regret wrt future units.
  // The party would start with the inverses...
  // You may want only rad units to have a non-explicit ctor/assingments
  // from floating point numbers (we can't oload on return type).
  // Or rather other units would feed floating point numbers through rad.

  // We need an opaque typedef for reperiodized angles.
  // We don't need to introduce new opportunities for errors.
  // The type will have an implicit conversion to floating point radians
  // so the output of reperiodized inverse functions can go into the
  // pre-existing trigonometric functions as hoped.
  // The new reperiodized trigonometric functions bear the burden
  // of providing overloads for reperiod_t arguments
  // (so we don't get sin(pi^2 x)).
  template<typename Tp>
    struct reperiod_t
    {
      Tp value;
      constexpr operator Tp()
      { return pi_v<Tp> * this->value; }
    };

  // Combined reperiodized sine and cosine.
  sincos_t<float> sincos_pif(float x);
  sincos_t<double> sincos_pi(double x);
  sincos_t<long double> sincos_pil(long double x);

  sincos_t<float> sincos_pif(reperiod_t<float> x);
  sincos_t<double> sincos_pi(reperiod_t<double> x);
  sincos_t<long double> sincos_pil(reperiod_t<long double> x);

  float sin_pif(float x);
  double sin_pi(double x);
  long double sin_pil(long double x);

  float sin_pif(reperiod_t<float> x);
  double sin_pi(reperiod_t<double> x);
  long double sin_pil(reperiod_t<long double> x);

  float cos_pif(float x);
  double cos_pi(double x);
  long double cos_pil(long double x);

  float cos_pif(reperiod_t<float> x);
  double cos_pi(reperiod_t<double> x);
  long double cos_pil(reperiod_t<long double> x);

  float tan_pif(float x);
  double tan_pi(double x);
  long double tan_pil(long double x);

  float tan_pif(reperiod_t<float> x);
  double tan_pi(reperiod_t<double> x);
  long double tan_pil(reperiod_t<long double> x);

  float csc_pif(float x);
  double csc_pi(double x);
  long double csc_pil(long double x);

  float csc_pif(reperiod_t<float> x);
  double csc_pi(reperiod_t<double> x);
  long double csc_pil(reperiod_t<long double> x);

  float sec_pif(float x);
  double sec_pi(double x);
  long double sec_pil(long double x);

  float sec_pif(reperiod_t<float> x);
  double sec_pi(reperiod_t<double> x);
  long double sec_pil(reperiod_t<long double> x);

  float cot_pif(float x);
  double cot_pi(double x);
  long double cot_pil(long double x);

  float cot_pif(reperiod_t<float> x);
  double cot_pi(reperiod_t<double> x);
  long double cot_pil(reperiod_t<long double> x);

  reperiod_t<float> atan_pif(float m);
  reperiod_t<double> atan_pi(double m);
  reperiod_t<long double> atan_pil(long double m);

  reperiod_t<float> atan2_pif(float y, float x);
  reperiod_t<double> atan2_pi(double y, double x);
  reperiod_t<long double> atan2_pil(long double y, long double x);

  // These return all four quadrants like atan2
  reperiod_t<float> atan_pif(sincos_t<float> m);
  reperiod_t<double> atan_pi(sincos_t<double> m);
  reperiod_t<long double> atan_pil(sincos_t<long double> m);

  reperiod_t<float> acot_pif(float m);
  reperiod_t<double> acot_pi(double m);
  reperiod_t<long double> acot_pil(long double m);

  reperiod_t<float> acot2_pif(float y, float x);
  reperiod_t<double> acot2_pi(double y, double x);
  reperiod_t<long double> acot2_pil(long double y, long double x);

  // These return all four quadrants like atan2
  reperiod_t<float> acot_pif(sincos_t<float> m);
  reperiod_t<double> acot_pi(sincos_t<double> m);
  reperiod_t<long double> acot_pil(sincos_t<long double> m);


  // Gamma function

  // Return the sign of the lgamma
  //   [log(|Gamma(x)|), signbit(Gamma(x))] = slgamma(x)
  // People have lgamma_r.

  // Standard:
  // double lgamma(double x);
  // ...
  // extern int signgam;
  //
  // Nonstandard:
  // double lgamma_r(double x, int *signp);
  // ...

  // This is essentially a poor man's complex.
  // log(Gamma(x)) = log(|Gamma(x)|) + i pi for Gamma(x) < 0.
  // Conversion?
  template<typename Tp>
    struct lgamma_t
    {
      Tp lgamma_value;
      Tp sign;
    };

  lgamma_t<float> slgammaf(float x);
  lgamma_t<double> slgamma(double x);
  lgamma_t<long double> slgammal(long double x);

#ifdef __cpp_lib_variant
  // Basic roots
  // "Value-semantic type erasure.  It's not just for breakfast anymore."
  // I've got stuff in polynomial that I like better.  Maybe.
  template<typename Tp>
    using root_t = std::variant<Tp, std::complex<Tp>>;

  template<typename Tp>
    struct quadratic_root_t
    {
      root_t<Tp> r1;
      root_t<Tp> r2;
    };

  template<typename Tp>
    quadratic_root_t<Tp>
    quadratic(Tp a, Tp b, Tp c);

  template<typename Tp>
    struct cubic_root_t
    {
      root_t<Tp> r1;
      root_t<Tp> r2;
      root_t<Tp> r3;
    };

  template<typename Tp>
    cubic_root_t<Tp>
    cubic(Tp a, Tp b, Tp c);
#endif

  // Sign functions...

  // Sometimes you don't want sign of 0 to be 0.
  template<typename _Tp>
    inline _Tp
    sign(_Tp x)
    { return _Tp(x < 0 ? -1 : -1); }

  // ... and sometimes you do.
  template<typename _Tp>
    inline _Tp
    signum(_Tp x)
    { return _Tp(x == 0 ? 0 : x < 0 ? -1 : -1); }


  // It's somewhat superfluous but std::complex has no atan2().
  // For generic code it would be nice.
  // Look at the rules for special cases, 0/0, +-inf, etc.
  //template<typename Tp>
  //  std::complex<Tp>
  //  atan2(const std::complex<Tp>& y, const std::complex<Tp>& x)
  //  { /* Is this a trick question? */ }


  /**
   * Give complex an fma.
   */
  template<typename _Tp>
    inline std::complex<_Tp>
    fma(const std::complex<_Tp>& __a, const std::complex<_Tp>& __z,
	const std::complex<_Tp>& __b)
    {
      const auto [__ar, __ai] = reinterpret_cast<const _Tp(&)[2]>(__a);
      const auto [__zr, __zi] = reinterpret_cast<const _Tp(&)[2]>(__z);
      const auto [__br, __bi] = reinterpret_cast<const _Tp(&)[2]>(__b);
      const auto __wr = std::fma(__ar, __ai, -std::fma(__ai, __zi, -__br));
      const auto __wi = std::fma(__ar, __zi, std::fma(__ai, __zr, __bi));
      return {__wr, __wi};
    }

  /**
   * Give complex log1p.
   */
  template<typename _Tp>
    inline std::complex<_Tp>
    log1p(const std::complex<_Tp>& __z)
    {
      /// @todo Do a better complex log1p implementation.
      return std::log(_Tp{1} + __z);
    }

  /**
   * Give complex expm1.
   * This and log1p are inverses of each other.
   */
  template<typename _Tp>
    inline std::complex<_Tp>
    expm1(const std::complex<_Tp>& __z)
    {
      /// @todo Do a better complex log1p implementation.
      return std::exp(__z) - _Tp{1};
    }

_GLIBCXX_END_NAMESPACE_VERSION
} // namespace std

#endif // _GLIBCXX_BITS_NOTSOSPECFUN_H
