
#include <algorithm>

#include <boost/python.hpp>

#include <functions/functions.h>
#include <interop/pyd/e_float_pyd_adapt.h>

namespace e_float_pyd
{
  namespace pyd
  {
    typedef ::e_float    e_float;
    typedef ::ef_complex ef_complex;

    struct ef
    {
      // <functions/constants/constants.h>
      static const pyd::e_float& zero         (void) { return ::ef::zero         (); }
      static const pyd::e_float& one          (void) { return ::ef::one          (); }
      static const pyd::e_float& half         (void) { return ::ef::half         (); }
      static const pyd::e_float& value_min    (void) { return ::ef::value_min    (); }
      static const pyd::e_float& value_max    (void) { return ::ef::value_max    (); }
      static const pyd::e_float& value_eps    (void) { return ::ef::value_eps    (); }
      static const pyd::e_float& value_inf    (void) { return ::ef::value_inf    (); }
      static const pyd::e_float& value_nan    (void) { return ::ef::value_nan    (); }
      static const pyd::e_float& two          (void) { return ::ef::two          (); }
      static const pyd::e_float& three        (void) { return ::ef::three        (); }
      static const pyd::e_float& four         (void) { return ::ef::four         (); }
      static const pyd::e_float& five         (void) { return ::ef::five         (); }
      static const pyd::e_float& six          (void) { return ::ef::six          (); }
      static const pyd::e_float& seven        (void) { return ::ef::seven        (); }
      static const pyd::e_float& eight        (void) { return ::ef::eight        (); }
      static const pyd::e_float& nine         (void) { return ::ef::nine         (); }
      static const pyd::e_float& ten          (void) { return ::ef::ten          (); }
      static const pyd::e_float& twenty       (void) { return ::ef::twenty       (); }
      static const pyd::e_float& thirty       (void) { return ::ef::thirty       (); }
      static const pyd::e_float& forty        (void) { return ::ef::forty        (); }
      static const pyd::e_float& fifty        (void) { return ::ef::fifty        (); }
      static const pyd::e_float& hundred      (void) { return ::ef::hundred      (); }
      static const pyd::e_float& two_hundred  (void) { return ::ef::two_hundred  (); }
      static const pyd::e_float& three_hundred(void) { return ::ef::three_hundred(); }
      static const pyd::e_float& four_hundred (void) { return ::ef::four_hundred (); }
      static const pyd::e_float& five_hundred (void) { return ::ef::five_hundred (); }
      static const pyd::e_float& thousand     (void) { return ::ef::thousand     (); }
      static const pyd::e_float& two_k        (void) { return ::ef::two_k        (); }
      static const pyd::e_float& three_k      (void) { return ::ef::three_k      (); }
      static const pyd::e_float& four_k       (void) { return ::ef::four_k       (); }
      static const pyd::e_float& five_k       (void) { return ::ef::five_k       (); }
      static const pyd::e_float& ten_k        (void) { return ::ef::ten_k        (); }
      static const pyd::e_float& twenty_k     (void) { return ::ef::twenty_k     (); }
      static const pyd::e_float& thirty_k     (void) { return ::ef::thirty_k     (); }
      static const pyd::e_float& forty_k      (void) { return ::ef::forty_k      (); }
      static const pyd::e_float& fifty_k      (void) { return ::ef::fifty_k      (); }
      static const pyd::e_float& hundred_k    (void) { return ::ef::hundred_k    (); }
      static const pyd::e_float& million      (void) { return ::ef::million      (); }
      static const pyd::e_float& ten_M        (void) { return ::ef::ten_M        (); }
      static const pyd::e_float& hundred_M    (void) { return ::ef::hundred_M    (); }
      static const pyd::e_float& billion      (void) { return ::ef::billion      (); }
      static const pyd::e_float& trillion     (void) { return ::ef::trillion     (); }
      static const pyd::e_float& googol       (void) { return ::ef::googol       (); }
      static const pyd::e_float& int32max     (void) { return ::ef::int32max     (); }
      static const pyd::e_float& int32min     (void) { return ::ef::int32min     (); }
      static const pyd::e_float& int64max     (void) { return ::ef::int64max     (); }
      static const pyd::e_float& int64min     (void) { return ::ef::int64min     (); }
      static const pyd::e_float& one_minus    (void) { return ::ef::one_minus    (); }
      static const pyd::e_float& tenth        (void) { return ::ef::tenth        (); }
      static const pyd::e_float& eighth       (void) { return ::ef::eighth       (); }
      static const pyd::e_float& fifth        (void) { return ::ef::fifth        (); }
      static const pyd::e_float& quarter      (void) { return ::ef::quarter      (); }
      static const pyd::e_float& third        (void) { return ::ef::third        (); }
      static const pyd::e_float& two_third    (void) { return ::ef::two_third    (); }
      static const pyd::e_float& four_third   (void) { return ::ef::four_third   (); }
      static const pyd::e_float& three_half   (void) { return ::ef::three_half   (); }
      static const pyd::e_float& sqrt2        (void) { return ::ef::sqrt2        (); }
      static const pyd::e_float& sqrt3        (void) { return ::ef::sqrt3        (); }
      static const pyd::e_float& pi           (void) { return ::ef::pi           (); }
      static const pyd::e_float& pi_half      (void) { return ::ef::pi_half      (); }
      static const pyd::e_float& pi_quarter   (void) { return ::ef::pi_quarter   (); }
      static const pyd::e_float& pi_squared   (void) { return ::ef::pi_squared   (); }
      static const pyd::e_float& two_pi       (void) { return ::ef::two_pi       (); }
      static const pyd::e_float& sqrt_pi      (void) { return ::ef::sqrt_pi      (); }
      static const pyd::e_float& degree       (void) { return ::ef::degree       (); }
      static const pyd::e_float& exp1         (void) { return ::ef::exp1         (); }
      static const pyd::e_float& ln2          (void) { return ::ef::ln2          (); }
      static const pyd::e_float& ln3          (void) { return ::ef::ln3          (); }
      static const pyd::e_float& ln10         (void) { return ::ef::ln10         (); }
      static const pyd::e_float& log10_2      (void) { return ::ef::log10_2      (); }
      static const pyd::e_float& golden_ratio (void) { return ::ef::golden_ratio (); }
      static const pyd::e_float& euler_gamma  (void) { return ::ef::euler_gamma  (); }
      static const pyd::e_float& catalan      (void) { return ::ef::catalan      (); }
      static const pyd::e_float& khinchin     (void) { return ::ef::khinchin     (); }
      static const pyd::e_float& glaisher     (void) { return ::ef::glaisher     (); }

      // <functions/elementary/elementary_math.h>
      static INT32 max_iteration(void) { return ::ef::max_iteration(); }

      static e_float floor(const pyd::e_float& x) { return ::ef::floor(x); }
      static e_float ceil (const pyd::e_float& x) { return ::ef::ceil (x); }
      static INT32   sgn  (const pyd::e_float& x) { return ::ef::sgn  (x); }

      static bool isnan   (const pyd::e_float& x)    { return ::ef::isnan(x); }
      static bool isnan   (const pyd::ef_complex& z) { return ::ef::isnan(z); }
      static bool isfinite(const pyd::e_float& x)    { return ::ef::isfinite(x); }
      static bool isfinite(const pyd::ef_complex& z) { return ::ef::isfinite(z); }
      static bool isinf   (const pyd::e_float& x)    { return ::ef::isinf(x); }
      static bool isinf   (const pyd::ef_complex& z) { return ::ef::isinf(z); }
      static bool isneg   (const pyd::e_float& x)    { return ::ef::isneg(x); }
      static bool isneg   (const pyd::ef_complex& z) { return ::ef::isneg(z); }
      static bool ispos   (const pyd::e_float& x)    { return ::ef::ispos(x); }
      static bool ispos   (const pyd::ef_complex& z) { return ::ef::ispos(z); }
      static bool isint   (const pyd::e_float& x)    { return ::ef::isint(x); }
      static bool isint   (const pyd::ef_complex& z) { return ::ef::isint(z); }
      static bool isone   (const pyd::e_float& x)    { return ::ef::isone(x); }
      static bool isone   (const pyd::ef_complex& z) { return ::ef::isone(z); }
      static bool iszero  (const pyd::e_float& x)    { return ::ef::iszero(x); }
      static bool iszero  (const pyd::ef_complex& z) { return ::ef::iszero(z); }

      static INT64 tol(void) { return ::ef::tol(); }

      static e_float fabs(const pyd::e_float& x)   { return ::ef::fabs(x); }
      static e_float abs (const pyd::e_float& x)   { return ::ef::fabs(x); }
      static e_float real(const pyd::e_float& x)   { return x; }
      static e_float imag(const pyd::e_float& x)   { static_cast<void>(x); return ::ef::zero(); }

      static e_float integer_part(const pyd::e_float& x) { return ::ef::integer_part(x); }
      static e_float decimal_part(const pyd::e_float& x) { return ::ef::decimal_part(x); }

      static void to_parts(const pyd::e_float& x, double& mantissa, INT64& exponent) { ::ef::to_parts(x, mantissa, exponent); }

      static double to_double(const pyd::e_float& x)    { return ::ef::to_double(x); }
      static double to_double(const pyd::ef_complex& z) { return ::ef::to_double(z); }

      static INT64 to_int64(const pyd::e_float& x)      { return ::ef::to_int64(x); }
      static INT64 to_int64(const pyd::ef_complex& z)   { return ::ef::to_int64(z); }
      static INT32 to_int32(const pyd::e_float& x)      { return ::ef::to_int32(x); }
      static INT32 to_int32(const pyd::ef_complex& z)   { return ::ef::to_int32(z); }

      static bool small_arg(const pyd::e_float& x)      { return ::ef::small_arg(x); }
      static bool small_arg(const pyd::ef_complex& z)   { return ::ef::small_arg(z); }
      static bool large_arg(const pyd::e_float& x)      { return ::ef::large_arg(x); }
      static bool large_arg(const pyd::ef_complex& z)   { return ::ef::large_arg(z); }
      static bool near_one (const pyd::e_float& x)      { return ::ef::near_one(x); }
      static bool near_one (const pyd::ef_complex& z)   { return ::ef::near_one(z); }
      static bool near_int (const pyd::e_float& x)      { return ::ef::near_int(x); }
      static bool near_int (const pyd::ef_complex& z)   { return ::ef::near_int(z); }

      // <functions/elementary/elementary_trig.h>
      static void         sincos(const pyd::e_float& x, pyd::e_float* s, pyd::e_float* c) { ::ef::sincos(x, s, c); }
      static pyd::e_float sin   (const pyd::e_float& x) { return ::ef::sin (x); }
      static pyd::e_float cos   (const pyd::e_float& x) { return ::ef::cos (x); }
      static pyd::e_float tan   (const pyd::e_float& x) { return ::ef::tan (x); }
      static pyd::e_float csc   (const pyd::e_float& x) { return ::ef::csc (x); }
      static pyd::e_float sec   (const pyd::e_float& x) { return ::ef::sec (x); }
      static pyd::e_float cot   (const pyd::e_float& x) { return ::ef::cot (x); }
      static pyd::e_float asin  (const pyd::e_float& x) { return ::ef::asin(x); }
      static pyd::e_float acos  (const pyd::e_float& x) { return ::ef::acos(x); }
      static pyd::e_float atan  (const pyd::e_float& x) { return ::ef::atan(x); }
      static pyd::e_float atan2 (const pyd::e_float& y, const pyd::e_float& x) { return ::ef::atan2(y, x); }

      // <functions/elementary/elementary_trans.h>
      static pyd::e_float pow2    (const INT64 p)                                           { return ::ef::pow2(p); }
      static pyd::e_float pown    (const pyd::e_float& x, const INT64 p)                    { return ::ef::pown(x, p); }
      static pyd::e_float inv     (const pyd::e_float& x)                                   { return ::ef::inv(x); }
      static pyd::e_float sqrt    (const pyd::e_float& x)                                   { return ::ef::sqrt(x); }
      static pyd::e_float cbrt    (const pyd::e_float& x)                                   { return ::ef::cbrt(x); }
      static pyd::e_float rootn   (const pyd::e_float& x, const INT32 p)                    { return ::ef::rootn(x, p); }
      static pyd::e_float exp     (const pyd::e_float& x)                                   { return ::ef::exp(x); }
      static pyd::e_float log     (const pyd::e_float& x)                                   { return ::ef::log(x); }
      static pyd::e_float log10   (const pyd::e_float& x)                                   { return ::ef::log10(x); }
      static pyd::e_float loga    (const pyd::e_float& a, const pyd::e_float& x)            { return ::ef::loga(a, x); }
      static pyd::e_float log1p   (const pyd::e_float& x)                                   { return ::ef::log1p(x); }
      static pyd::e_float log1p1m2(const pyd::e_float& x)                                   { return ::ef::log1p1m2(x); }
      static pyd::e_float pow     (const pyd::e_float& x, const pyd::e_float& a)            { return ::ef::pow(x, a); }
      static void         sinhcosh(const pyd::e_float& x, pyd::e_float* s, pyd::e_float* c) { ::ef::sinhcosh(x, s, c); }
      static pyd::e_float sinh    (const pyd::e_float& x)                                   { return ::ef::sinh(x); }
      static pyd::e_float cosh    (const pyd::e_float& x)                                   { return ::ef::cosh(x); }
      static pyd::e_float tanh    (const pyd::e_float& x)                                   { return ::ef::tanh(x); }
      static pyd::e_float asinh   (const pyd::e_float& x)                                   { return ::ef::asinh(x); }
      static pyd::e_float acosh   (const pyd::e_float& x)                                   { return ::ef::acosh(x); }
      static pyd::e_float atanh   (const pyd::e_float& x)                                   { return ::ef::atanh(x); }

      // <functions/bessel/airy.h>
      static pyd::e_float airy_a      (const pyd::e_float& x) { return ::ef::airy_a(x); }
      static pyd::e_float airy_a_prime(const pyd::e_float& x) { return ::ef::airy_a_prime(x); }
      static pyd::e_float airy_b      (const pyd::e_float& x) { return ::ef::airy_b(x); }
      static pyd::e_float airy_b_prime(const pyd::e_float& x) { return ::ef::airy_b_prime(x); }

      static boost::python::list airy_a_zero(const UINT32 k) { return pyd::adapt::std_deque_to_pylist(::ef::airy_a_zero(k)); }
      static boost::python::list airy_b_zero(const UINT32 k) { return pyd::adapt::std_deque_to_pylist(::ef::airy_b_zero(k)); }

      // <functions/bessel/bessel.h>
      static pyd::e_float cyl_bessel_i(const pyd::e_float& v, const pyd::e_float& x) { return ::ef::cyl_bessel_i(v, x); }
      static pyd::e_float cyl_bessel_i(const INT32 n,         const pyd::e_float& x) { return ::ef::cyl_bessel_i(n, x); }
      static pyd::e_float cyl_bessel_j(const pyd::e_float& v, const pyd::e_float& x) { return ::ef::cyl_bessel_j(v, x); }
      static pyd::e_float cyl_bessel_j(const INT32 n,         const pyd::e_float& x) { return ::ef::cyl_bessel_j(n, x); }
      static pyd::e_float cyl_bessel_k(const pyd::e_float& v, const pyd::e_float& x) { return ::ef::cyl_bessel_k(v, x); }
      static pyd::e_float cyl_bessel_k(const INT32 n,         const pyd::e_float& x) { return ::ef::cyl_bessel_k(n, x); }
      static pyd::e_float cyl_bessel_y(const pyd::e_float& v, const pyd::e_float& x) { return ::ef::cyl_bessel_y(v, x); }
      static pyd::e_float cyl_bessel_y(const INT32 n,         const pyd::e_float& x) { return ::ef::cyl_bessel_y(n, x); }

      static pyd::e_float cyl_bessel_i_prime(const pyd::e_float& v, const pyd::e_float& x) { return ::ef::cyl_bessel_i_prime(v, x); }
      static pyd::e_float cyl_bessel_i_prime(const INT32 n,         const pyd::e_float& x) { return ::ef::cyl_bessel_i_prime(n, x); }
      static pyd::e_float cyl_bessel_j_prime(const pyd::e_float& v, const pyd::e_float& x) { return ::ef::cyl_bessel_j_prime(v, x); }
      static pyd::e_float cyl_bessel_j_prime(const INT32 n,         const pyd::e_float& x) { return ::ef::cyl_bessel_j_prime(n, x); }
      static pyd::e_float cyl_bessel_k_prime(const pyd::e_float& v, const pyd::e_float& x) { return ::ef::cyl_bessel_k_prime(v, x); }
      static pyd::e_float cyl_bessel_k_prime(const INT32 n,         const pyd::e_float& x) { return ::ef::cyl_bessel_k_prime(n, x); }
      static pyd::e_float cyl_bessel_y_prime(const pyd::e_float& v, const pyd::e_float& x) { return ::ef::cyl_bessel_y_prime(v, x); }
      static pyd::e_float cyl_bessel_y_prime(const INT32 n,         const pyd::e_float& x) { return ::ef::cyl_bessel_y_prime(n, x); }

      static boost::python::list cyl_bessel_j_zero(const pyd::e_float& v, const UINT32 k) { return pyd::adapt::std_deque_to_pylist(::ef::cyl_bessel_j_zero(v, k)); }
      static boost::python::list cyl_bessel_j_zero(const INT32 n, const UINT32 k)         { return pyd::adapt::std_deque_to_pylist(::ef::cyl_bessel_j_zero(n, k)); }

      // <functions/elliptic/elliptic.h>
      static pyd::e_float comp_ellint_1(const pyd::e_float& m)                          { return ::ef::comp_ellint_1(m); }
      static pyd::e_float      ellint_1(const pyd::e_float& m, const pyd::e_float& phi) { return ::ef::ellint_1(m, phi); }
      static pyd::e_float comp_ellint_2(const pyd::e_float& m)                          { return ::ef::comp_ellint_2(m); }
      static pyd::e_float      ellint_2(const pyd::e_float& m, const pyd::e_float& phi) { return ::ef::ellint_2(m, phi); }

      // <functions/gamma/gamma.h>
      static pyd::e_float gamma               (const pyd::e_float& x)                                                 { return ::ef::gamma(x); } 
      static pyd::e_float gamma_near_n        (const INT32 n, const pyd::e_float& x)                                  { return ::ef::gamma_near_n(n, x);}
      static pyd::e_float incomplete_gamma    (const pyd::e_float& a, const pyd::e_float& x)                          { return ::ef::incomplete_gamma(a, x); }
      static pyd::e_float gen_incomplete_gamma(const pyd::e_float& a, const pyd::e_float& x0, const pyd::e_float& x1) { return ::ef::gen_incomplete_gamma(a, x0, x1); }
      static pyd::e_float beta                (const pyd::e_float& a, const pyd::e_float& b)                          { return ::ef::beta(a, b); }
      static pyd::e_float incomplete_beta     (const pyd::e_float& x, const pyd::e_float& a, const pyd::e_float& b)   { return ::ef::incomplete_beta(x, a, b); }
      static pyd::e_float factorial           (const UINT32 n)                                                        { return ::ef::factorial(n); }
      static pyd::e_float factorial2          (const  INT32 n)                                                        { return ::ef::factorial2(n); }
      static pyd::e_float binomial            (const UINT32 n, const UINT32 k)                                        { return ::ef::binomial(n, k); }
      static pyd::e_float binomial            (const UINT32 n, const pyd::e_float& y)                                 { return ::ef::binomial(n, y); }
      static pyd::e_float binomial            (const pyd::e_float& x, const UINT32 k)                                 { return ::ef::binomial(x, k); }
      static pyd::e_float binomial            (const pyd::e_float& x, const pyd::e_float& y)                          { return ::ef::binomial(x, y); }
      static pyd::e_float pochhammer          (const pyd::e_float& x, const UINT32 n)                                 { return ::ef::pochhammer(x, n); }
      static pyd::e_float pochhammer          (const pyd::e_float& x, const pyd::e_float& a)                          { return ::ef::pochhammer(x, a); }

      // <functions/gamma/poly_gamma.h>
      static pyd::e_float poly_gamma(const pyd::e_float& x)                { return ::ef::poly_gamma(x); }
      static pyd::e_float poly_gamma(const INT32 n, const pyd::e_float& x) { return ::ef::poly_gamma(n, x); }

      // <functions/hypergeometric/hermite.h>
      // and
      // <functions/polynomials/polynomials.h>
      static pyd::e_float hermite(const pyd::e_float& v, const pyd::e_float& x)                          { return ::ef::hermite(v, x); }
      static pyd::e_float hermite(const INT32 n,         const pyd::e_float& x)                          { return ::ef::hermite(n, x); }
      static pyd::e_float hermite(const UINT32 n,        const pyd::e_float& x, boost::python::list& lp) { std::vector<pyd::e_float> vp; const pyd::e_float y = ::ef::hermite(n, x, &vp); lp = pyd::adapt::std_vector_to_pylist(vp); return y; }

      // <functions/hypergeometric/hypergeometric.h>
      static pyd::e_float hyperg_0f0    (const pyd::e_float& x)                                                                      { return ::ef::hyperg_0f0(x); }
      static pyd::e_float hyperg_0f1    (const pyd::e_float& b, const pyd::e_float& x)                                               { return ::ef::hyperg_0f1(b, x); }
      static pyd::e_float hyperg_0f1_reg(const pyd::e_float& a, const pyd::e_float& x)                                               { return ::ef::hyperg_0f1_reg(a, x); }
      static pyd::e_float hyperg_1f0    (const pyd::e_float& a, const pyd::e_float& x)                                               { return ::ef::hyperg_1f0(a, x); }
      static pyd::e_float hyperg_1f1    (const pyd::e_float& a, const pyd::e_float& b, const pyd::e_float& x)                        { return ::ef::hyperg_1f1(a, b, x); }
      static pyd::e_float hyperg_1f1_reg(const pyd::e_float& a, const pyd::e_float& b, const pyd::e_float& x)                        { return ::ef::hyperg_1f1_reg(a, b, x); }
      static pyd::e_float hyperg_2f0    (const pyd::e_float& a, const pyd::e_float& b, const pyd::e_float& x)                        { return ::ef::hyperg_2f0(a, b, x); }
      static pyd::e_float hyperg_2f1    (const pyd::e_float& a, const pyd::e_float& b, const pyd::e_float& c, const pyd::e_float& x) { return ::ef::hyperg_2f1(a, b, c, x); }
      static pyd::e_float hyperg_2f1_reg(const pyd::e_float& a, const pyd::e_float& b, const pyd::e_float& c, const pyd::e_float& x) { return ::ef::hyperg_2f1_reg(a, b, c, x); }
      static pyd::e_float hyperg_pfq    (const boost::python::list& a, const boost::python::list& b, const pyd::e_float& x)          { return ::ef::hyperg_pfq(pyd::adapt::pylist_to_std_deque<pyd::e_float>(a), pyd::adapt::pylist_to_std_deque<pyd::e_float>(b), x); }
      static pyd::e_float conf_hyperg   (const pyd::e_float& a, const pyd::e_float& c, const pyd::e_float& x)                        { return ::ef::hyperg_1f1(a, c, x); }
      static pyd::e_float      hyperg   (const pyd::e_float& a, const pyd::e_float& b, const pyd::e_float& c, const pyd::e_float& x) { return ::ef::hyperg_2f1(a, b, c, x); }

      // <functions/hypergeometric/laguerre.h>
      // and
      // <functions/polynomials/polynomials.h>
      static pyd::e_float laguerre(const pyd::e_float& v, const pyd::e_float& x)                          { return ::ef::laguerre(v, x); }
      static pyd::e_float laguerre(const pyd::e_float& v, const pyd::e_float& L, const pyd::e_float& x)   { return ::ef::laguerre(v, L, x); }
      static pyd::e_float laguerre(const INT32 n,         const pyd::e_float& L, const pyd::e_float& x)   { return ::ef::laguerre(n, L, x); }
      static pyd::e_float laguerre(const INT32 n,         const INT32 m,         const pyd::e_float& x)   { return ::ef::laguerre(n, m, x); }
      static pyd::e_float laguerre(const INT32 n,         const pyd::e_float& x)                          { return ::ef::laguerre(n, x); }
      static pyd::e_float laguerre(const UINT32 n,        const pyd::e_float& x, boost::python::list& lp) { std::vector<pyd::e_float> vp; const pyd::e_float y = ::ef::laguerre(n, x, &vp); lp = pyd::adapt::std_vector_to_pylist(vp); return y; }

      // <functions/hypergeometric/legendre.h>
      // and
      // <functions/polynomials/polynomials.h>
      static pyd::e_float legendre_p(const pyd::e_float& v, const pyd::e_float& x)                        { return ::ef::legendre_p(v, x); }
      static pyd::e_float legendre_p(const pyd::e_float& v, const pyd::e_float& u, const pyd::e_float& x) { return ::ef::legendre_p(v, u, x); }
      static pyd::e_float legendre_p(const pyd::e_float& v, const INT32   m,       const pyd::e_float& x) { return ::ef::legendre_p(v, m, x); }
      static pyd::e_float legendre_p(const INT32 n,         const pyd::e_float& u, const pyd::e_float& x) { return ::ef::legendre_p(n, u, x); }
      static pyd::e_float legendre_p(const INT32 n,         const INT32   m,       const pyd::e_float& x) { return ::ef::legendre_p(n, m, x); }
      static pyd::e_float legendre_p(const INT32 n,         const pyd::e_float& x)                        { return ::ef::legendre_p(n, x); }
      static pyd::e_float legendre_q(const pyd::e_float& v, const pyd::e_float& x)                        { return ::ef::legendre_q(v, x); }
      static pyd::e_float legendre_q(const pyd::e_float& v, const pyd::e_float& u, const pyd::e_float& x) { return ::ef::legendre_q(v, u, x); }
      static pyd::e_float legendre_q(const pyd::e_float& v, const INT32   m,       const pyd::e_float& x) { return ::ef::legendre_q(v, m, x); }
      static pyd::e_float legendre_q(const INT32 n,         const pyd::e_float& u, const pyd::e_float& x) { return ::ef::legendre_q(n, u, x); }
      static pyd::e_float legendre_q(const INT32 n,         const INT32   m,       const pyd::e_float& x) { return ::ef::legendre_q(n, m, x); }
      static pyd::e_float legendre_q(const INT32 n,         const pyd::e_float& x)                        { return ::ef::legendre_q(n, x); }

      // <functions/hypergeometric/parabolic_cylinder_d.h>
      static pyd::e_float weber_d(const pyd::e_float& v, const pyd::e_float& x) { return ::ef::weber_d(v, x); }

      // <functions/integer/integer.h>
      static pyd::e_float euler    (const UINT32 n)                 { return ::ef::euler(n); }
      static pyd::e_float bernoulli(const UINT32 n)                 { return ::ef::bernoulli(n); }
      static pyd::e_float stirling2(const UINT32 n, const UINT32 k) { return ::ef::stirling2(n, k); }

      // <functions/integer/prime.h>
      static boost::python::list prime(const UINT32 n) { std::deque<UINT32> pd; ::ef::prime(n, pd); return pyd::adapt::std_deque_to_pylist(pd); }

      // <functions/polynomials/polynomials.h>
      static pyd::e_float chebyshev_t(const INT32 n,  const pyd::e_float& x)                          { return ::ef::chebyshev_t(n, x); }
      static pyd::e_float chebyshev_t(const UINT32 n, const pyd::e_float& x, boost::python::list& lp) { std::vector<pyd::e_float> vp; const pyd::e_float y = ::ef::chebyshev_t(n, x, &vp); lp = pyd::adapt::std_vector_to_pylist(vp); return y; }
      static pyd::e_float chebyshev_u(const INT32 n,  const pyd::e_float& x)                          { return ::ef::chebyshev_u(n, x); }
      static pyd::e_float chebyshev_u(const UINT32 n, const pyd::e_float& x, boost::python::list& lp) { std::vector<pyd::e_float> vp; const pyd::e_float y = ::ef::chebyshev_u(n, x, &vp); lp = pyd::adapt::std_vector_to_pylist(vp); return y; }

      // <functions/zeta/polylog.h>
      static pyd::e_float poly_logarithm(const INT32 n, const pyd::e_float& x) { return ::ef::poly_logarithm(n, x); }

      // <functions/zeta/zeta.h>
      static pyd::e_float riemann_zeta(const INT32 n)                                { return ::ef::riemann_zeta(n); }
      static pyd::e_float riemann_zeta(const pyd::e_float& s)                        { return ::ef::riemann_zeta(s); }
      static pyd::e_float hurwitz_zeta(const pyd::e_float& s, const INT32 n)         { return ::ef::hurwitz_zeta(s, n); }
      static pyd::e_float hurwitz_zeta(const pyd::e_float& s, const pyd::e_float& a) { return ::ef::hurwitz_zeta(s, a); }
    };

    struct efz
    {
      // <functions/elementary/elementary_complex.h>
      static pyd::ef_complex norm    (const pyd::ef_complex& z)                                         { return z.norm(); }
      static pyd::ef_complex abs     (const pyd::ef_complex& z)                                         { return ::efz::abs(z); }
      static pyd::ef_complex arg     (const pyd::ef_complex& z)                                         { return ::efz::arg(z); }
      static pyd::ef_complex real    (const pyd::ef_complex& z)                                         { return z.real(); }
      static pyd::ef_complex imag    (const pyd::ef_complex& z)                                         { return z.imag(); }
      static pyd::ef_complex conj    (const pyd::ef_complex& z)                                         { return ::efz::conj(z); }
      static pyd::ef_complex iz      (const pyd::ef_complex& z)                                         { return ::efz::iz(z); }
      static pyd::ef_complex polar   (const pyd::e_float& mod, const pyd::e_float& arg)                 { return ::efz::polar(mod, arg); }
      static pyd::ef_complex sin     (const pyd::ef_complex& z)                                         { return ::efz::sin(z); }
      static pyd::ef_complex cos     (const pyd::ef_complex& z)                                         { return ::efz::cos(z); }
      static pyd::ef_complex tan     (const pyd::ef_complex& z)                                         { return ::efz::tan(z); }
      static void            sincos  (const pyd::ef_complex& z, pyd::ef_complex* s, pyd::ef_complex* c) { ::efz::sincos(z, s, c); }
      static pyd::ef_complex csc     (const pyd::ef_complex& z)                                         { return ::efz::csc  (z); }
      static pyd::ef_complex sec     (const pyd::ef_complex& z)                                         { return ::efz::sec  (z); }
      static pyd::ef_complex cot     (const pyd::ef_complex& z)                                         { return ::efz::cot  (z); }
      static pyd::ef_complex asin    (const pyd::ef_complex& z)                                         { return ::efz::asin (z); }
      static pyd::ef_complex acos    (const pyd::ef_complex& z)                                         { return ::efz::acos (z); }
      static pyd::ef_complex atan    (const pyd::ef_complex& z)                                         { return ::efz::atan (z); }
      static pyd::ef_complex inv     (const pyd::ef_complex& z)                                         { return ::efz::inv  (z); }
      static pyd::ef_complex sqrt    (const pyd::ef_complex& z)                                         { return ::efz::sqrt (z); }
      static pyd::ef_complex exp     (const pyd::ef_complex& z)                                         { return ::efz::exp  (z); }
      static pyd::ef_complex log     (const pyd::ef_complex& z)                                         { return ::efz::log  (z); }
      static pyd::ef_complex log10   (const pyd::ef_complex& z)                                         { return ::efz::log10(z); }
      static pyd::ef_complex loga    (const pyd::ef_complex& a, const pyd::ef_complex& z)               { return ::efz::loga(a, z); }
      static pyd::ef_complex pown    (const pyd::ef_complex& z, const INT64 p)                          { return ::efz::pown(z, p); }
      static pyd::ef_complex pow     (const pyd::ef_complex& z, const pyd::ef_complex& a)               { return ::efz::pow(z, a); }
      static pyd::ef_complex rootn   (const pyd::ef_complex& z, const INT32 p)                          { return ::efz::rootn(z, p); }
      static pyd::ef_complex sinh    (const pyd::ef_complex& z)                                         { return ::efz::sinh(z); }
      static pyd::ef_complex cosh    (const pyd::ef_complex& z)                                         { return ::efz::cosh(z); }
      static pyd::ef_complex tanh    (const pyd::ef_complex& z)                                         { return ::efz::tanh(z); }
      static void            sinhcosh(const pyd::ef_complex& z, pyd::ef_complex* s, pyd::ef_complex* c) { ::efz::sinhcosh(z, s, c); }
      static pyd::ef_complex asinh   (const pyd::ef_complex& z)                                         { return ::efz::asinh(z); }
      static pyd::ef_complex acosh   (const pyd::ef_complex& z)                                         { return ::efz::acosh(z); }
      static pyd::ef_complex atanh   (const pyd::ef_complex& z)                                         { return ::efz::atanh(z); }

      // <functions/gamma/gamma.h>
      static pyd::ef_complex gamma     (const pyd::ef_complex& z)                           { return ::efz::gamma(z); }
      static pyd::ef_complex beta      (const pyd::ef_complex& a, const pyd::ef_complex& b) { return ::efz::beta(a, b); }
      static pyd::ef_complex pochhammer(const pyd::ef_complex& z, const UINT32 n)           { return ::efz::pochhammer(z, n); }
      static pyd::ef_complex pochhammer(const pyd::ef_complex& z, const pyd::ef_complex& a) { return ::efz::pochhammer(z, a); }

      // <functions/polynomials/polynomials.h>
      static pyd::ef_complex chebyshev_t(const INT32 n, const pyd::ef_complex& z) { return ::efz::chebyshev_t(n, z); }
      static pyd::ef_complex chebyshev_u(const INT32 n, const pyd::ef_complex& z) { return ::efz::chebyshev_u(n, z); }
      static pyd::ef_complex hermite    (const INT32 n, const pyd::ef_complex& z) { return ::efz::hermite    (n, z); }
      static pyd::ef_complex laguerre   (const INT32 n, const pyd::ef_complex& z) { return ::efz::laguerre   (n, z); }
      static pyd::ef_complex chebyshev_t(const UINT32 n, const pyd::ef_complex& z, boost::python::list& lp) { std::vector<pyd::ef_complex> vp; const pyd::ef_complex y = ::efz::chebyshev_t(n, z, &vp); lp = pyd::adapt::std_vector_to_pylist(vp); return y; }
      static pyd::ef_complex chebyshev_u(const UINT32 n, const pyd::ef_complex& z, boost::python::list& lp) { std::vector<pyd::ef_complex> vp; const pyd::ef_complex y = ::efz::chebyshev_u(n, z, &vp); lp = pyd::adapt::std_vector_to_pylist(vp); return y; }
      static pyd::ef_complex hermite    (const UINT32 n, const pyd::ef_complex& z, boost::python::list& lp) { std::vector<pyd::ef_complex> vp; const pyd::ef_complex y = ::efz::hermite    (n, z, &vp); lp = pyd::adapt::std_vector_to_pylist(vp); return y; }
      static pyd::ef_complex laguerre   (const UINT32 n, const pyd::ef_complex& z, boost::python::list& lp) { std::vector<pyd::ef_complex> vp; const pyd::ef_complex y = ::efz::laguerre   (n, z, &vp); lp = pyd::adapt::std_vector_to_pylist(vp); return y; }

      // <functions/zeta/zeta.h>
      static pyd::ef_complex riemann_zeta(const pyd::ef_complex& s)                           { return ::efz::riemann_zeta(s); }
      static pyd::ef_complex hurwitz_zeta(const pyd::ef_complex& s, const INT32 n)            { return ::efz::hurwitz_zeta(s, n); }
      static pyd::ef_complex hurwitz_zeta(const pyd::ef_complex& s, const pyd::ef_complex& a) { return ::efz::hurwitz_zeta(s, a); }
    };
  }
}

using namespace e_float_pyd;

BOOST_PYTHON_MODULE(e_float_pyd)
{
  boost::python::class_<pyd::e_float>("e_float", "e_float type")
    .def(boost::python::init<const pyd::e_float&>())
    .def(boost::python::init<INT32>())
    .def(boost::python::init<UINT32>())
    .def(boost::python::init<INT64>())
    .def(boost::python::init<UINT64>())
    .def(boost::python::init<double>())
    .def(boost::python::init<const std::string&>())
    .def("__pos__",   &pyd::e_float::get_pos)
    .def("__neg__",   &pyd::e_float::get_neg)
    .def("__abs__",   &pyd::e_float::get_abs)
    .def("__int__",   &pyd::e_float::get_int)
    .def("__float__", &pyd::e_float::get_flt)
    .def("__str__",   &pyd::e_float::get_str)
    .def(boost::python::self += boost::python::self)
    .def(boost::python::self -= boost::python::self)
    .def(boost::python::self *= boost::python::self)
    .def(boost::python::self /= boost::python::self)
    .def(boost::python::self + pyd::e_float())
    .def(boost::python::self - pyd::e_float())
    .def(boost::python::self * pyd::e_float())
    .def(boost::python::self / pyd::e_float())
    .def(pyd::e_float() + boost::python::self)
    .def(pyd::e_float() - boost::python::self)
    .def(pyd::e_float() * boost::python::self)
    .def(pyd::e_float() / boost::python::self)
    .def(boost::python::self += INT32())
    .def(boost::python::self -= INT32())
    .def(boost::python::self *= INT32())
    .def(boost::python::self /= INT32())
    .def(boost::python::self + INT32())
    .def(boost::python::self - INT32())
    .def(boost::python::self * INT32())
    .def(boost::python::self / INT32())
    .def(INT32() + boost::python::self)
    .def(INT32() - boost::python::self)
    .def(INT32() * boost::python::self)
    .def(INT32() / boost::python::self)
    .def(boost::python::self <  boost::python::self)
    .def(boost::python::self <= boost::python::self)
    .def(boost::python::self == boost::python::self)
    .def(boost::python::self != boost::python::self)
    .def(boost::python::self >= boost::python::self)
    .def(boost::python::self >  boost::python::self)
    .def(INT32() <  boost::python::self)
    .def(INT32() <= boost::python::self)
    .def(INT32() == boost::python::self)
    .def(INT32() != boost::python::self)
    .def(INT32() >= boost::python::self)
    .def(INT32() >  boost::python::self)
    .def(boost::python::self <  INT32())
    .def(boost::python::self <= INT32())
    .def(boost::python::self == INT32())
    .def(boost::python::self != INT32())
    .def(boost::python::self >= INT32())
    .def(boost::python::self >  INT32())
    .def("isnan",    &pyd::e_float::isnan)
    .def("isinf",    &pyd::e_float::isinf)
    .def("isfinite", &pyd::e_float::isfinite)
    .def("iszero",   &pyd::e_float::iszero)
    .def("isone",    &pyd::e_float::isone)
    .def("isint",    &pyd::e_float::isint)
    .def("isneg",    &pyd::e_float::isneg)
    .def("order",    &pyd::e_float::order)
    ;

  boost::python::class_<pyd::ef_complex>("ef_complex", "ef_complex type")
    .def(boost::python::init<const pyd::ef_complex&>())
    .def(boost::python::init<const pyd::e_float&>())
    .def(boost::python::init<const pyd::e_float&, const pyd::e_float&>())
    .def(boost::python::init<INT32>())
    .def(boost::python::init<UINT32>())
    .def(boost::python::init<INT64>())
    .def(boost::python::init<UINT64>())
    .def(boost::python::init<double>())
    .def("__pos__",   &pyd::ef_complex::get_pos)
    .def("__neg__",   &pyd::ef_complex::get_neg)
    .def("__abs__",   &pyd::ef_complex::get_abs)
    .def("__int__",   &pyd::ef_complex::get_int)
    .def("__float__", &pyd::ef_complex::get_flt)
    .def("__str__",   &pyd::ef_complex::get_str)
    .def(boost::python::self += boost::python::self)
    .def(boost::python::self -= boost::python::self)
    .def(boost::python::self *= boost::python::self)
    .def(boost::python::self /= boost::python::self)
    .def(boost::python::self + pyd::ef_complex())
    .def(boost::python::self - pyd::ef_complex())
    .def(boost::python::self * pyd::ef_complex())
    .def(boost::python::self / pyd::ef_complex())
    .def(pyd::ef_complex() + boost::python::self)
    .def(pyd::ef_complex() - boost::python::self)
    .def(pyd::ef_complex() * boost::python::self)
    .def(pyd::ef_complex() / boost::python::self)
    .def(boost::python::self + pyd::e_float())
    .def(boost::python::self - pyd::e_float())
    .def(boost::python::self * pyd::e_float())
    .def(boost::python::self / pyd::e_float())
    .def(pyd::e_float() + boost::python::self)
    .def(pyd::e_float() - boost::python::self)
    .def(pyd::e_float() * boost::python::self)
    .def(pyd::e_float() / boost::python::self)
    .def(boost::python::self += INT32())
    .def(boost::python::self -= INT32())
    .def(boost::python::self *= INT32())
    .def(boost::python::self /= INT32())
    .def(boost::python::self + INT32())
    .def(boost::python::self - INT32())
    .def(boost::python::self * INT32())
    .def(boost::python::self / INT32())
    .def(INT32() + boost::python::self)
    .def(INT32() - boost::python::self)
    .def(INT32() * boost::python::self)
    .def(INT32() / boost::python::self)
    .def(boost::python::self == boost::python::self)
    .def(boost::python::self != boost::python::self)
    .def(boost::python::self == pyd::e_float())
    .def(boost::python::self != pyd::e_float())
    .def(pyd::e_float() == boost::python::self)
    .def(pyd::e_float() != boost::python::self)
    .def(INT32() == boost::python::self)
    .def(INT32() != boost::python::self)
    .def(boost::python::self == INT32())
    .def(boost::python::self != INT32())
    .def("isnan",    &pyd::ef_complex::isnan)
    .def("isinf",    &pyd::ef_complex::isinf)
    .def("isfinite", &pyd::ef_complex::isfinite)
    .def("iszero",   &pyd::ef_complex::iszero)
    .def("isone",    &pyd::ef_complex::isone)
    .def("isint",    &pyd::ef_complex::isint)
    .def("isneg",    &pyd::ef_complex::isneg)
    ;

  boost::python::class_<e_float_pyd::pyd::ef>("ef", "ef interface structure")
    // <functions/constants/constants.h>
    .def("zero",          &pyd::ef::zero         , boost::python::return_value_policy<boost::python::copy_const_reference>()) .staticmethod("zero")
    .def("one",           &pyd::ef::one          , boost::python::return_value_policy<boost::python::copy_const_reference>()) .staticmethod("one")
    .def("half",          &pyd::ef::half         , boost::python::return_value_policy<boost::python::copy_const_reference>()) .staticmethod("half")
    .def("value_min",     &pyd::ef::value_min    , boost::python::return_value_policy<boost::python::copy_const_reference>()) .staticmethod("value_min")
    .def("value_max",     &pyd::ef::value_max    , boost::python::return_value_policy<boost::python::copy_const_reference>()) .staticmethod("value_max")
    .def("value_eps",     &pyd::ef::value_eps    , boost::python::return_value_policy<boost::python::copy_const_reference>()) .staticmethod("value_eps")
    .def("value_inf",     &pyd::ef::value_inf    , boost::python::return_value_policy<boost::python::copy_const_reference>()) .staticmethod("value_inf")
    .def("value_nan",     &pyd::ef::value_nan    , boost::python::return_value_policy<boost::python::copy_const_reference>()) .staticmethod("value_nan")
    .def("two",           &pyd::ef::two          , boost::python::return_value_policy<boost::python::copy_const_reference>()) .staticmethod("two")
    .def("three",         &pyd::ef::three        , boost::python::return_value_policy<boost::python::copy_const_reference>()) .staticmethod("three")
    .def("four",          &pyd::ef::four         , boost::python::return_value_policy<boost::python::copy_const_reference>()) .staticmethod("four")
    .def("five",          &pyd::ef::five         , boost::python::return_value_policy<boost::python::copy_const_reference>()) .staticmethod("five")
    .def("six",           &pyd::ef::six          , boost::python::return_value_policy<boost::python::copy_const_reference>()) .staticmethod("six")
    .def("seven",         &pyd::ef::seven        , boost::python::return_value_policy<boost::python::copy_const_reference>()) .staticmethod("seven")
    .def("eight",         &pyd::ef::eight        , boost::python::return_value_policy<boost::python::copy_const_reference>()) .staticmethod("eight")
    .def("nine",          &pyd::ef::nine         , boost::python::return_value_policy<boost::python::copy_const_reference>()) .staticmethod("nine")
    .def("ten",           &pyd::ef::ten          , boost::python::return_value_policy<boost::python::copy_const_reference>()) .staticmethod("ten")
    .def("twenty",        &pyd::ef::twenty       , boost::python::return_value_policy<boost::python::copy_const_reference>()) .staticmethod("twenty")
    .def("thirty",        &pyd::ef::thirty       , boost::python::return_value_policy<boost::python::copy_const_reference>()) .staticmethod("thirty")
    .def("forty",         &pyd::ef::forty        , boost::python::return_value_policy<boost::python::copy_const_reference>()) .staticmethod("forty")
    .def("fifty",         &pyd::ef::fifty        , boost::python::return_value_policy<boost::python::copy_const_reference>()) .staticmethod("fifty")
    .def("hundred",       &pyd::ef::hundred      , boost::python::return_value_policy<boost::python::copy_const_reference>()) .staticmethod("hundred")
    .def("two_hundred",   &pyd::ef::two_hundred  , boost::python::return_value_policy<boost::python::copy_const_reference>()) .staticmethod("two_hundred")
    .def("three_hundred", &pyd::ef::three_hundred, boost::python::return_value_policy<boost::python::copy_const_reference>()) .staticmethod("three_hundred")
    .def("four_hundred",  &pyd::ef::four_hundred , boost::python::return_value_policy<boost::python::copy_const_reference>()) .staticmethod("four_hundred")
    .def("five_hundred",  &pyd::ef::five_hundred , boost::python::return_value_policy<boost::python::copy_const_reference>()) .staticmethod("five_hundred")
    .def("thousand",      &pyd::ef::thousand     , boost::python::return_value_policy<boost::python::copy_const_reference>()) .staticmethod("thousand")
    .def("two_k",         &pyd::ef::two_k        , boost::python::return_value_policy<boost::python::copy_const_reference>()) .staticmethod("two_k")
    .def("three_k",       &pyd::ef::three_k      , boost::python::return_value_policy<boost::python::copy_const_reference>()) .staticmethod("three_k")
    .def("four_k",        &pyd::ef::four_k       , boost::python::return_value_policy<boost::python::copy_const_reference>()) .staticmethod("four_k")
    .def("five_k",        &pyd::ef::five_k       , boost::python::return_value_policy<boost::python::copy_const_reference>()) .staticmethod("five_k")
    .def("ten_k",         &pyd::ef::ten_k        , boost::python::return_value_policy<boost::python::copy_const_reference>()) .staticmethod("ten_k")
    .def("twenty_k",      &pyd::ef::twenty_k     , boost::python::return_value_policy<boost::python::copy_const_reference>()) .staticmethod("twenty_k")
    .def("thirty_k",      &pyd::ef::thirty_k     , boost::python::return_value_policy<boost::python::copy_const_reference>()) .staticmethod("thirty_k")
    .def("forty_k",       &pyd::ef::forty_k      , boost::python::return_value_policy<boost::python::copy_const_reference>()) .staticmethod("forty_k")
    .def("fifty_k",       &pyd::ef::fifty_k      , boost::python::return_value_policy<boost::python::copy_const_reference>()) .staticmethod("fifty_k")
    .def("hundred_k",     &pyd::ef::hundred_k    , boost::python::return_value_policy<boost::python::copy_const_reference>()) .staticmethod("hundred_k")
    .def("million",       &pyd::ef::million      , boost::python::return_value_policy<boost::python::copy_const_reference>()) .staticmethod("million")
    .def("ten_M",         &pyd::ef::ten_M        , boost::python::return_value_policy<boost::python::copy_const_reference>()) .staticmethod("ten_M")
    .def("hundred_M",     &pyd::ef::hundred_M    , boost::python::return_value_policy<boost::python::copy_const_reference>()) .staticmethod("hundred_M")
    .def("billion",       &pyd::ef::billion      , boost::python::return_value_policy<boost::python::copy_const_reference>()) .staticmethod("billion")
    .def("trillion",      &pyd::ef::trillion     , boost::python::return_value_policy<boost::python::copy_const_reference>()) .staticmethod("trillion")
    .def("googol",        &pyd::ef::googol       , boost::python::return_value_policy<boost::python::copy_const_reference>()) .staticmethod("googol")
    .def("int32max",      &pyd::ef::int32max     , boost::python::return_value_policy<boost::python::copy_const_reference>()) .staticmethod("int32max")
    .def("int32min",      &pyd::ef::int32min     , boost::python::return_value_policy<boost::python::copy_const_reference>()) .staticmethod("int32min")
    .def("int64max",      &pyd::ef::int64max     , boost::python::return_value_policy<boost::python::copy_const_reference>()) .staticmethod("int64max")
    .def("int64min",      &pyd::ef::int64min     , boost::python::return_value_policy<boost::python::copy_const_reference>()) .staticmethod("int64min")
    .def("one_minus",     &pyd::ef::one_minus    , boost::python::return_value_policy<boost::python::copy_const_reference>()) .staticmethod("one_minus")
    .def("tenth",         &pyd::ef::tenth        , boost::python::return_value_policy<boost::python::copy_const_reference>()) .staticmethod("tenth")
    .def("eighth",        &pyd::ef::eighth       , boost::python::return_value_policy<boost::python::copy_const_reference>()) .staticmethod("eighth")
    .def("fifth",         &pyd::ef::fifth        , boost::python::return_value_policy<boost::python::copy_const_reference>()) .staticmethod("fifth")
    .def("quarter",       &pyd::ef::quarter      , boost::python::return_value_policy<boost::python::copy_const_reference>()) .staticmethod("quarter")
    .def("third",         &pyd::ef::third        , boost::python::return_value_policy<boost::python::copy_const_reference>()) .staticmethod("third")
    .def("two_third",     &pyd::ef::two_third    , boost::python::return_value_policy<boost::python::copy_const_reference>()) .staticmethod("two_third")
    .def("four_third",    &pyd::ef::four_third   , boost::python::return_value_policy<boost::python::copy_const_reference>()) .staticmethod("four_third")
    .def("three_half",    &pyd::ef::three_half   , boost::python::return_value_policy<boost::python::copy_const_reference>()) .staticmethod("three_half")
    .def("sqrt2",         &pyd::ef::sqrt2        , boost::python::return_value_policy<boost::python::copy_const_reference>()) .staticmethod("sqrt2")
    .def("sqrt3",         &pyd::ef::sqrt3        , boost::python::return_value_policy<boost::python::copy_const_reference>()) .staticmethod("sqrt3")
    .def("pi",            &pyd::ef::pi           , boost::python::return_value_policy<boost::python::copy_const_reference>()) .staticmethod("pi")
    .def("pi_half",       &pyd::ef::pi_half      , boost::python::return_value_policy<boost::python::copy_const_reference>()) .staticmethod("pi_half")
    .def("pi_quarter",    &pyd::ef::pi_quarter   , boost::python::return_value_policy<boost::python::copy_const_reference>()) .staticmethod("pi_quarter")
    .def("pi_squared",    &pyd::ef::pi_squared   , boost::python::return_value_policy<boost::python::copy_const_reference>()) .staticmethod("pi_squared")
    .def("two_pi",        &pyd::ef::two_pi       , boost::python::return_value_policy<boost::python::copy_const_reference>()) .staticmethod("two_pi")
    .def("sqrt_pi",       &pyd::ef::sqrt_pi      , boost::python::return_value_policy<boost::python::copy_const_reference>()) .staticmethod("sqrt_pi")
    .def("degree",        &pyd::ef::degree       , boost::python::return_value_policy<boost::python::copy_const_reference>()) .staticmethod("degree")
    .def("exp1",          &pyd::ef::exp1         , boost::python::return_value_policy<boost::python::copy_const_reference>()) .staticmethod("exp1")
    .def("ln2",           &pyd::ef::ln2          , boost::python::return_value_policy<boost::python::copy_const_reference>()) .staticmethod("ln2")
    .def("ln3",           &pyd::ef::ln3          , boost::python::return_value_policy<boost::python::copy_const_reference>()) .staticmethod("ln3")
    .def("ln10",          &pyd::ef::ln10         , boost::python::return_value_policy<boost::python::copy_const_reference>()) .staticmethod("ln10")
    .def("log10_2",       &pyd::ef::log10_2      , boost::python::return_value_policy<boost::python::copy_const_reference>()) .staticmethod("log10_2")
    .def("golden_ratio",  &pyd::ef::golden_ratio , boost::python::return_value_policy<boost::python::copy_const_reference>()) .staticmethod("golden_ratio")
    .def("euler_gamma",   &pyd::ef::euler_gamma  , boost::python::return_value_policy<boost::python::copy_const_reference>()) .staticmethod("euler_gamma")
    .def("catalan",       &pyd::ef::catalan      , boost::python::return_value_policy<boost::python::copy_const_reference>()) .staticmethod("catalan")
    .def("khinchin",      &pyd::ef::khinchin     , boost::python::return_value_policy<boost::python::copy_const_reference>()) .staticmethod("khinchin")
    .def("glaisher",      &pyd::ef::glaisher     , boost::python::return_value_policy<boost::python::copy_const_reference>()) .staticmethod("glaisher")

    // <functions/elementary/elementary_math.h>
    .def("max_iteration", &pyd::ef::max_iteration)  .staticmethod("max_iteration")
    .def("floor",         &pyd::ef::floor)          .staticmethod("floor")
    .def("ceil",          &pyd::ef::ceil)           .staticmethod("ceil")
    .def("sgn",           &pyd::ef::sgn)            .staticmethod("sgn")

    .def("isnan",         static_cast<bool(*)(const pyd::e_float&)>   (&pyd::ef::isnan))
    .def("isnan",         static_cast<bool(*)(const pyd::ef_complex&)>(&pyd::ef::isnan))
    .staticmethod("isnan")
    .def("isfinite",      static_cast<bool(*)(const pyd::e_float&)>   (&pyd::ef::isfinite))
    .def("isfinite",      static_cast<bool(*)(const pyd::ef_complex&)>(&pyd::ef::isfinite))
    .staticmethod("isfinite")
    .def("isinf",         static_cast<bool(*)(const pyd::e_float&)>   (&pyd::ef::isinf))
    .def("isinf",         static_cast<bool(*)(const pyd::ef_complex&)>(&pyd::ef::isinf))
    .staticmethod("isinf")
    .def("isneg",         static_cast<bool(*)(const pyd::e_float&)>   (&pyd::ef::isneg))
    .def("isneg",         static_cast<bool(*)(const pyd::ef_complex&)>(&pyd::ef::isneg))
    .staticmethod("isneg")
    .def("ispos",         static_cast<bool(*)(const pyd::e_float&)>   (&pyd::ef::ispos))
    .def("ispos",         static_cast<bool(*)(const pyd::ef_complex&)>(&pyd::ef::ispos))
    .staticmethod("ispos")
    .def("isint",         static_cast<bool(*)(const pyd::e_float&)>   (&pyd::ef::isint))
    .def("isint",         static_cast<bool(*)(const pyd::ef_complex&)>(&pyd::ef::isint))
    .staticmethod("isint")
    .def("isone",         static_cast<bool(*)(const pyd::e_float&)>   (&pyd::ef::isone))
    .def("isone",         static_cast<bool(*)(const pyd::ef_complex&)>(&pyd::ef::isone))
    .staticmethod("isone")
    .def("iszero",        static_cast<bool(*)(const pyd::e_float&)>   (&pyd::ef::iszero))
    .def("iszero",        static_cast<bool(*)(const pyd::ef_complex&)>(&pyd::ef::iszero))
    .staticmethod("iszero")

    .def("tol",           &pyd::ef::tol)            .staticmethod("tol")
    .def("fabs",          &pyd::ef::fabs)           .staticmethod("fabs")
    .def("abs",           &pyd::ef::abs)            .staticmethod("abs")
    .def("real",          &pyd::ef::real)           .staticmethod("real")
    .def("imag",          &pyd::ef::imag)           .staticmethod("imag")
    .def("integer_part",  &pyd::ef::integer_part)   .staticmethod("integer_part")
    .def("decimal_part",  &pyd::ef::decimal_part)   .staticmethod("decimal_part")
    .def("to_parts",      &pyd::ef::to_parts)       .staticmethod("to_parts")

    .def("to_double",        static_cast<double(*)(const pyd::e_float&)>   (&pyd::ef::to_double))
    .def("to_double",        static_cast<double(*)(const pyd::ef_complex&)>(&pyd::ef::to_double))
    .staticmethod("to_double")
    .def("to_int64",         static_cast<INT64 (*)(const pyd::e_float&)>   (&pyd::ef::to_int64))
    .def("to_int64",         static_cast<INT64 (*)(const pyd::ef_complex&)>(&pyd::ef::to_int64))
    .staticmethod("to_int64")
    .def("to_int32",         static_cast<INT32 (*)(const pyd::e_float&)>   (&pyd::ef::to_int32))
    .def("to_int32",         static_cast<INT32 (*)(const pyd::ef_complex&)>(&pyd::ef::to_int32))
    .staticmethod("to_int32")
    .def("small_arg",        static_cast<bool  (*)(const pyd::e_float&)>   (&pyd::ef::small_arg))
    .def("small_arg",        static_cast<bool  (*)(const pyd::ef_complex&)>(&pyd::ef::small_arg))
    .staticmethod("small_arg")
    .def("large_arg",        static_cast<bool  (*)(const pyd::e_float&)>   (&pyd::ef::large_arg))
    .def("large_arg",        static_cast<bool  (*)(const pyd::ef_complex&)>(&pyd::ef::large_arg))
    .staticmethod("large_arg")
    .def("near_one",         static_cast<bool  (*)(const pyd::e_float&)>   (&pyd::ef::near_one))
    .def("near_one",         static_cast<bool  (*)(const pyd::ef_complex&)>(&pyd::ef::near_one))
    .staticmethod("near_one")
    .def("near_int",         static_cast<bool  (*)(const pyd::e_float&)>   (&pyd::ef::near_int))
    .def("near_int",         static_cast<bool  (*)(const pyd::ef_complex&)>(&pyd::ef::near_int))
    .staticmethod("near_int")

    // <functions/elementary/elementary_trig.h>
    .def("sincos", &pyd::ef::sincos) .staticmethod("sincos")
    .def("sin",    &pyd::ef::sin)    .staticmethod("sin")
    .def("cos",    &pyd::ef::cos)    .staticmethod("cos")
    .def("tan",    &pyd::ef::tan)    .staticmethod("tan")
    .def("csc",    &pyd::ef::csc)    .staticmethod("csc")
    .def("sec",    &pyd::ef::sec)    .staticmethod("sec")
    .def("cot",    &pyd::ef::cot)    .staticmethod("cot")
    .def("asin",   &pyd::ef::asin)   .staticmethod("asin")
    .def("acos",   &pyd::ef::acos)   .staticmethod("acos")
    .def("atan",   &pyd::ef::atan)   .staticmethod("atan")
    .def("atan2",  &pyd::ef::atan2)  .staticmethod("atan2")

    // <functions/elementary/elementary_trans.h>
    .def("pow2",     &pyd::ef::pow2)     .staticmethod("pow2")
    .def("pown",     &pyd::ef::pown)     .staticmethod("pown")
    .def("inv",      &pyd::ef::inv)      .staticmethod("inv")
    .def("sqrt",     &pyd::ef::sqrt)     .staticmethod("sqrt")
    .def("cbrt",     &pyd::ef::cbrt)     .staticmethod("cbrt")
    .def("rootn",    &pyd::ef::rootn)    .staticmethod("rootn")
    .def("exp",      &pyd::ef::exp)      .staticmethod("exp")
    .def("log",      &pyd::ef::log)      .staticmethod("log")
    .def("log10",    &pyd::ef::log10)    .staticmethod("log10")
    .def("loga",     &pyd::ef::loga)     .staticmethod("loga")
    .def("log1p",    &pyd::ef::log1p)    .staticmethod("log1p")
    .def("log1p1m2", &pyd::ef::log1p1m2) .staticmethod("log1p1m2")
    .def("pow",      &pyd::ef::pow)      .staticmethod("pow")
    .def("sinhcosh", &pyd::ef::sinhcosh) .staticmethod("sinhcosh")
    .def("sinh",     &pyd::ef::sinh)     .staticmethod("sinh")
    .def("cosh",     &pyd::ef::cosh)     .staticmethod("cosh")
    .def("tanh",     &pyd::ef::tanh)     .staticmethod("tanh")
    .def("asinh",    &pyd::ef::asinh)    .staticmethod("asinh")
    .def("acosh",    &pyd::ef::acosh)    .staticmethod("acosh")
    .def("atanh",    &pyd::ef::atanh)    .staticmethod("atanh")

    // <functions/bessel/airy.h>
    .def("airy_a",       &pyd::ef::airy_a)        .staticmethod("airy_a")
    .def("airy_a_prime", &pyd::ef::airy_a_prime)  .staticmethod("airy_a_prime")
    .def("airy_b",       &pyd::ef::airy_b)        .staticmethod("airy_b")
    .def("airy_b_prime", &pyd::ef::airy_b_prime)  .staticmethod("airy_b_prime")
    .def("airy_a_zero",  &pyd::ef::airy_a_zero)   .staticmethod("airy_a_zero")
    .def("airy_b_zero",  &pyd::ef::airy_b_zero)   .staticmethod("airy_b_zero")

    // <functions/bessel/bessel.h>
    .def("cyl_bessel_i", static_cast<pyd::e_float(*)(const pyd::e_float&, const pyd::e_float&)>(&pyd::ef::cyl_bessel_i))
    .def("cyl_bessel_i", static_cast<pyd::e_float(*)(const INT32, const pyd::e_float&)>        (&pyd::ef::cyl_bessel_i))
    .staticmethod("cyl_bessel_i")
    .def("cyl_bessel_j", static_cast<pyd::e_float(*)(const pyd::e_float&, const pyd::e_float&)>(&pyd::ef::cyl_bessel_j))
    .def("cyl_bessel_j", static_cast<pyd::e_float(*)(const INT32, const pyd::e_float&)>        (&pyd::ef::cyl_bessel_j))
    .staticmethod("cyl_bessel_j")
    .def("cyl_bessel_k", static_cast<pyd::e_float(*)(const pyd::e_float&, const pyd::e_float&)>(&pyd::ef::cyl_bessel_k))
    .def("cyl_bessel_k", static_cast<pyd::e_float(*)(const INT32, const pyd::e_float&)>        (&pyd::ef::cyl_bessel_k))
    .staticmethod("cyl_bessel_k")
    .def("cyl_bessel_y", static_cast<pyd::e_float(*)(const pyd::e_float&, const pyd::e_float&)>(&pyd::ef::cyl_bessel_y))
    .def("cyl_bessel_y", static_cast<pyd::e_float(*)(const INT32, const pyd::e_float&)>        (&pyd::ef::cyl_bessel_y))
    .staticmethod("cyl_bessel_y")
    .def("cyl_bessel_i_prime", static_cast<pyd::e_float(*)(const pyd::e_float&, const pyd::e_float&)>(&pyd::ef::cyl_bessel_i_prime))
    .def("cyl_bessel_i_prime", static_cast<pyd::e_float(*)(const INT32, const pyd::e_float&)>        (&pyd::ef::cyl_bessel_i_prime))
    .staticmethod("cyl_bessel_i_prime")
    .def("cyl_bessel_j_prime", static_cast<pyd::e_float(*)(const pyd::e_float&, const pyd::e_float&)>(&pyd::ef::cyl_bessel_j_prime))
    .def("cyl_bessel_j_prime", static_cast<pyd::e_float(*)(const INT32, const pyd::e_float&)>        (&pyd::ef::cyl_bessel_j_prime))
    .staticmethod("cyl_bessel_j_prime")
    .def("cyl_bessel_k_prime", static_cast<pyd::e_float(*)(const pyd::e_float&, const pyd::e_float&)>(&pyd::ef::cyl_bessel_k_prime))
    .def("cyl_bessel_k_prime", static_cast<pyd::e_float(*)(const INT32, const pyd::e_float&)>        (&pyd::ef::cyl_bessel_k_prime))
    .staticmethod("cyl_bessel_k_prime")
    .def("cyl_bessel_y_prime", static_cast<pyd::e_float(*)(const pyd::e_float&, const pyd::e_float&)>(&pyd::ef::cyl_bessel_y_prime))
    .def("cyl_bessel_y_prime", static_cast<pyd::e_float(*)(const INT32, const pyd::e_float&)>        (&pyd::ef::cyl_bessel_y_prime))
    .staticmethod("cyl_bessel_y_prime")
    .def("cyl_bessel_j_zero", static_cast<boost::python::list(*)(const pyd::e_float&, const UINT32)>(&pyd::ef::cyl_bessel_j_zero))
    .def("cyl_bessel_j_zero", static_cast<boost::python::list(*)(const INT32, const UINT32)>        (&pyd::ef::cyl_bessel_j_zero))
    .staticmethod("cyl_bessel_j_zero")

    // <functions/elliptic/elliptic.h>
    .def("comp_ellint_1", &pyd::ef::comp_ellint_1) .staticmethod("comp_ellint_1")
    .def("ellint_1",      &pyd::ef::ellint_1)      .staticmethod("ellint_1")
    .def("comp_ellint_2", &pyd::ef::comp_ellint_2) .staticmethod("comp_ellint_2")
    .def("ellint_2",      &pyd::ef::ellint_2)      .staticmethod("ellint_2")

    // <functions/gamma/gamma.h>
    .def("gamma",                &pyd::ef::gamma)                     .staticmethod("gamma")
    .def("gamma_near_n",         &pyd::ef::gamma_near_n)              .staticmethod("gamma_near_n")
    .def("incomplete_gamma",     &pyd::ef::incomplete_gamma)          .staticmethod("incomplete_gamma")
    .def("gen_incomplete_gamma", &pyd::ef::gen_incomplete_gamma)      .staticmethod("gen_incomplete_gamma")
    .def("beta",                 &pyd::ef::beta)                      .staticmethod("beta")
    .def("incomplete_beta",      &pyd::ef::incomplete_beta)           .staticmethod("incomplete_beta")
    .def("factorial",            &pyd::ef::factorial)                 .staticmethod("factorial")
    .def("factorial2",           &pyd::ef::factorial2)                .staticmethod("factorial2")
    .def("binomial", static_cast<pyd::e_float(*)(const UINT32,  const UINT32)>             (&pyd::ef::binomial))
    .def("binomial", static_cast<pyd::e_float(*)(const UINT32,  const pyd::e_float&)>      (&pyd::ef::binomial))
    .def("binomial", static_cast<pyd::e_float(*)(const pyd::e_float&, const UINT32)>       (&pyd::ef::binomial))
    .def("binomial", static_cast<pyd::e_float(*)(const pyd::e_float&, const pyd::e_float&)>(&pyd::ef::binomial))
    .staticmethod("binomial")
    .def("pochhammer", static_cast<pyd::e_float(*)(const pyd::e_float&, const UINT32)>       (&pyd::ef::pochhammer))
    .def("pochhammer", static_cast<pyd::e_float(*)(const pyd::e_float&, const pyd::e_float&)>(&pyd::ef::pochhammer))
    .staticmethod("pochhammer")

    // <functions/gamma/poly_gamma.h>
    .def("poly_gamma", static_cast<pyd::e_float(*)(const pyd::e_float&)>             (&pyd::ef::poly_gamma))
    .def("poly_gamma", static_cast<pyd::e_float(*)(const INT32, const pyd::e_float&)>(&pyd::ef::poly_gamma))
    .staticmethod("poly_gamma")

    // <functions/hypergeometric/hermite.h>
    // and
    // <functions/polynomials/polynomials.h>
    .def("hermite", static_cast<pyd::e_float(*)(const pyd::e_float&, const pyd::e_float&)>                      (&pyd::ef::hermite))
    .def("hermite", static_cast<pyd::e_float(*)(const INT32,         const pyd::e_float&)>                      (&pyd::ef::hermite))
    .def("hermite", static_cast<pyd::e_float(*)(const UINT32,        const pyd::e_float&, boost::python::list&)>(&pyd::ef::hermite))
    .staticmethod("hermite")

    // <functions/hypergeometric/hypergeometric.h>
    .def("hyperg_0f0",     &pyd::ef::hyperg_0f0)     .staticmethod("hyperg_0f0")
    .def("hyperg_0f1",     &pyd::ef::hyperg_0f1)     .staticmethod("hyperg_0f1")
    .def("hyperg_0f1_reg", &pyd::ef::hyperg_0f1_reg) .staticmethod("hyperg_0f1_reg")
    .def("hyperg_1f0",     &pyd::ef::hyperg_1f0)     .staticmethod("hyperg_1f0")
    .def("hyperg_1f1",     &pyd::ef::hyperg_1f1)     .staticmethod("hyperg_1f1")
    .def("hyperg_1f1_reg", &pyd::ef::hyperg_1f1_reg) .staticmethod("hyperg_1f1_reg")
    .def("hyperg_2f0",     &pyd::ef::hyperg_2f0)     .staticmethod("hyperg_2f0")
    .def("hyperg_2f1",     &pyd::ef::hyperg_2f1)     .staticmethod("hyperg_2f1")
    .def("hyperg_2f1_reg", &pyd::ef::hyperg_2f1_reg) .staticmethod("hyperg_2f1_reg")
    .def("hyperg_pfq",     &pyd::ef::hyperg_pfq)     .staticmethod("hyperg_pfq")
    .def("conf_hyperg",    &pyd::ef::conf_hyperg)    .staticmethod("conf_hyperg")
    .def("hyperg",         &pyd::ef::hyperg)         .staticmethod("hyperg")

    // <functions/hypergeometric/laguerre.h>
    // and
    // <functions/polynomials/polynomials.h>
    .def("laguerre", static_cast<pyd::e_float(*)(const pyd::e_float&, const pyd::e_float&)>                      (&pyd::ef::laguerre))
    .def("laguerre", static_cast<pyd::e_float(*)(const pyd::e_float&, const pyd::e_float&, const pyd::e_float&)> (&pyd::ef::laguerre))
    .def("laguerre", static_cast<pyd::e_float(*)(const INT32,         const pyd::e_float&, const pyd::e_float&)> (&pyd::ef::laguerre))
    .def("laguerre", static_cast<pyd::e_float(*)(const INT32,         const INT32,         const pyd::e_float&)> (&pyd::ef::laguerre))
    .def("laguerre", static_cast<pyd::e_float(*)(const INT32,         const pyd::e_float&)>                      (&pyd::ef::laguerre))
    .def("laguerre", static_cast<pyd::e_float(*)(const UINT32,        const pyd::e_float&, boost::python::list&)>(&pyd::ef::laguerre))
    .staticmethod("laguerre")

    // <functions/hypergeometric/legendre.h>
    // and
    // <functions/polynomials/polynomials.h>
    .def("legendre_p", static_cast<pyd::e_float(*)(const pyd::e_float&, const pyd::e_float&)>                     (&pyd::ef::legendre_p))
    .def("legendre_p", static_cast<pyd::e_float(*)(const pyd::e_float&, const pyd::e_float&, const pyd::e_float&)>(&pyd::ef::legendre_p))
    .def("legendre_p", static_cast<pyd::e_float(*)(const pyd::e_float&, const INT32,         const pyd::e_float&)>(&pyd::ef::legendre_p))
    .def("legendre_p", static_cast<pyd::e_float(*)(const INT32,         const pyd::e_float&, const pyd::e_float&)>(&pyd::ef::legendre_p))
    .def("legendre_p", static_cast<pyd::e_float(*)(const INT32,         const INT32,         const pyd::e_float&)>(&pyd::ef::legendre_p))
    .def("legendre_p", static_cast<pyd::e_float(*)(const INT32,         const pyd::e_float&)>                     (&pyd::ef::legendre_p))
    .staticmethod("legendre_p")
    .def("legendre_q", static_cast<pyd::e_float(*)(const pyd::e_float&, const pyd::e_float&)>                     (&pyd::ef::legendre_q))
    .def("legendre_q", static_cast<pyd::e_float(*)(const pyd::e_float&, const pyd::e_float&, const pyd::e_float&)>(&pyd::ef::legendre_q))
    .def("legendre_q", static_cast<pyd::e_float(*)(const pyd::e_float&, const INT32,         const pyd::e_float&)>(&pyd::ef::legendre_q))
    .def("legendre_q", static_cast<pyd::e_float(*)(const INT32,         const pyd::e_float&, const pyd::e_float&)>(&pyd::ef::legendre_q))
    .def("legendre_q", static_cast<pyd::e_float(*)(const INT32,         const INT32,         const pyd::e_float&)>(&pyd::ef::legendre_q))
    .def("legendre_q", static_cast<pyd::e_float(*)(const INT32,         const pyd::e_float&)>                     (&pyd::ef::legendre_q))
    .staticmethod("legendre_q")

    // <functions/hypergeometric/parabolic_cylinder_d.h>
    .def("weber_d", &pyd::ef::weber_d) .staticmethod("weber_d")

    // <functions/integer/integer.h>
    .def("euler",     &pyd::ef::euler)     .staticmethod("euler")
    .def("bernoulli", &pyd::ef::bernoulli) .staticmethod("bernoulli")
    .def("stirling2", &pyd::ef::stirling2) .staticmethod("stirling2")

    // <functions/integer/prime.h>
    .def("prime", &pyd::ef::prime) .staticmethod("prime")

    // <functions/polynomials/polynomials.h>
    .def("chebyshev_t", static_cast<pyd::e_float(*)(const INT32,  const pyd::e_float&)>                      (&pyd::ef::chebyshev_t))
    .def("chebyshev_t", static_cast<pyd::e_float(*)(const UINT32, const pyd::e_float&, boost::python::list&)>(&pyd::ef::chebyshev_t))
    .staticmethod("chebyshev_t")
    .def("chebyshev_u", static_cast<pyd::e_float(*)(const INT32,  const pyd::e_float&)>                      (&pyd::ef::chebyshev_u))
    .def("chebyshev_u", static_cast<pyd::e_float(*)(const UINT32, const pyd::e_float&, boost::python::list&)>(&pyd::ef::chebyshev_u))
    .staticmethod("chebyshev_u")

    // <functions/zeta/polylog.h>
    .def("poly_logarithm", &pyd::ef::poly_logarithm)  .staticmethod("poly_logarithm")

    // <functions/zeta/zeta.h>
    .def("riemann_zeta", static_cast<pyd::e_float(*)(const INT32)>        (&pyd::ef::riemann_zeta))
    .def("riemann_zeta", static_cast<pyd::e_float(*)(const pyd::e_float&)>(&pyd::ef::riemann_zeta))
    .staticmethod("riemann_zeta")
    .def("hurwitz_zeta", static_cast<pyd::e_float(*)(const pyd::e_float&, const INT32)>        (&pyd::ef::hurwitz_zeta))
    .def("hurwitz_zeta", static_cast<pyd::e_float(*)(const pyd::e_float&, const pyd::e_float&)>(&pyd::ef::hurwitz_zeta))
    .staticmethod("hurwitz_zeta")
    ;

  boost::python::class_<e_float_pyd::pyd::efz>("efz", "efz interface structure")
    // <functions/elementary/elementary_complex.h>
    .def("norm",     &pyd::efz::norm)      .staticmethod("norm")
    .def("abs",      &pyd::efz::abs)       .staticmethod("abs")
    .def("arg",      &pyd::efz::arg)       .staticmethod("arg")
    .def("real",     &pyd::efz::real)      .staticmethod("real")
    .def("imag",     &pyd::efz::imag)      .staticmethod("imag")
    .def("conj",     &pyd::efz::conj)      .staticmethod("conj")
    .def("iz",       &pyd::efz::iz)        .staticmethod("iz")
    .def("polar",    &pyd::efz::polar)     .staticmethod("polar")
    .def("sin",      &pyd::efz::sin)       .staticmethod("sin")
    .def("cos",      &pyd::efz::cos)       .staticmethod("cos")
    .def("tan",      &pyd::efz::tan)       .staticmethod("tan")
    .def("sincos",   &pyd::efz::sincos)    .staticmethod("sincos")
    .def("csc",      &pyd::efz::csc)       .staticmethod("csc")
    .def("sec",      &pyd::efz::sec)       .staticmethod("sec")
    .def("cot",      &pyd::efz::cot)       .staticmethod("cot")
    .def("asin",     &pyd::efz::asin)      .staticmethod("asin")
    .def("acos",     &pyd::efz::acos)      .staticmethod("acos")
    .def("atan",     &pyd::efz::atan)      .staticmethod("atan")
    .def("inv",      &pyd::efz::inv)       .staticmethod("inv")
    .def("sqrt",     &pyd::efz::sqrt)      .staticmethod("sqrt")
    .def("exp",      &pyd::efz::exp)       .staticmethod("exp")
    .def("log",      &pyd::efz::log)       .staticmethod("log")
    .def("log10",    &pyd::efz::log10)     .staticmethod("log10")
    .def("loga",     &pyd::efz::loga)      .staticmethod("loga")
    .def("pown",     &pyd::efz::pown)      .staticmethod("pown")
    .def("pow",      &pyd::efz::pow)       .staticmethod("pow")
    .def("rootn",    &pyd::efz::rootn)     .staticmethod("rootn")
    .def("sinh",     &pyd::efz::sinh)      .staticmethod("sinh")
    .def("cosh",     &pyd::efz::cosh)      .staticmethod("cosh")
    .def("tanh",     &pyd::efz::tanh)      .staticmethod("tanh")
    .def("sinhcosh", &pyd::efz::sinhcosh)  .staticmethod("sinhcosh")
    .def("asinh",    &pyd::efz::asinh)     .staticmethod("asinh")
    .def("acosh",    &pyd::efz::acosh)     .staticmethod("acosh")
    .def("atanh",    &pyd::efz::atanh)     .staticmethod("atanh")

    // <functions/gamma/gamma.h>
    .def("gamma",      &pyd::efz::gamma)     .staticmethod("gamma")
    .def("beta",       &pyd::efz::beta)      .staticmethod("beta")
    .def("pochhammer", static_cast<pyd::ef_complex(*)(const pyd::ef_complex&, const UINT32)>          (&pyd::efz::pochhammer))
    .def("pochhammer", static_cast<pyd::ef_complex(*)(const pyd::ef_complex&, const pyd::ef_complex&)>(&pyd::efz::pochhammer))
    .staticmethod("pochhammer")

      // <functions/polynomials/polynomials.h>
    .def("chebyshev_t", static_cast<pyd::ef_complex(*)(const INT32,  const pyd::ef_complex&)>                      (&pyd::efz::chebyshev_t))
    .def("chebyshev_t", static_cast<pyd::ef_complex(*)(const UINT32, const pyd::ef_complex&, boost::python::list&)>(&pyd::efz::chebyshev_t))
    .staticmethod("chebyshev_t")
    .def("chebyshev_u", static_cast<pyd::ef_complex(*)(const INT32,  const pyd::ef_complex&)>                      (&pyd::efz::chebyshev_u))
    .def("chebyshev_u", static_cast<pyd::ef_complex(*)(const UINT32, const pyd::ef_complex&, boost::python::list&)>(&pyd::efz::chebyshev_u))
    .staticmethod("chebyshev_u")
    .def("hermite",     static_cast<pyd::ef_complex(*)(const INT32,  const pyd::ef_complex&)>                      (&pyd::efz::hermite))
    .def("hermite",     static_cast<pyd::ef_complex(*)(const UINT32, const pyd::ef_complex&, boost::python::list&)>(&pyd::efz::hermite))
    .staticmethod("hermite")
    .def("laguerre",    static_cast<pyd::ef_complex(*)(const INT32,  const pyd::ef_complex&)>                      (&pyd::efz::laguerre))
    .def("laguerre",    static_cast<pyd::ef_complex(*)(const UINT32, const pyd::ef_complex&, boost::python::list&)>(&pyd::efz::laguerre))
    .staticmethod("laguerre")

      // <functions/zeta/zeta.h>
    .def("riemann_zeta", &pyd::efz::riemann_zeta)  .staticmethod("riemann_zeta")
    .def("hurwitz_zeta", static_cast<pyd::ef_complex(*)(const pyd::ef_complex&, const INT32)>           (&pyd::efz::hurwitz_zeta))
    .def("hurwitz_zeta", static_cast<pyd::ef_complex(*)(const pyd::ef_complex&, const pyd::ef_complex&)>(&pyd::efz::hurwitz_zeta))
    .staticmethod("hurwitz_zeta")
    ;
}
