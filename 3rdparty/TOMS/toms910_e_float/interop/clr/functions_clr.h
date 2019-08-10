
#ifndef _FUNCTIONS_CLR_2010_01_06_H_
  #define _FUNCTIONS_CLR_2010_01_06_H_

  #include <functions/functions.h>
  #include <interop/clr/e_float_clr.h>
  #include <interop/clr/e_float_complex_clr.h>
  #include <interop/clr/functions_base_clr.h>

  namespace e_float_clr
  {
    namespace cli
    {
      public value struct ef sealed : public functions_base
      {
        // <functions/constants/constants.h>
        static cli::e_float^ zero         (void) { return gcnew cli::e_float(::ef::zero         ()); }
        static cli::e_float^ one          (void) { return gcnew cli::e_float(::ef::one          ()); }
        static cli::e_float^ half         (void) { return gcnew cli::e_float(::ef::half         ()); }
        static cli::e_float^ value_min    (void) { return gcnew cli::e_float(::ef::value_min    ()); }
        static cli::e_float^ value_max    (void) { return gcnew cli::e_float(::ef::value_max    ()); }
        static cli::e_float^ value_eps    (void) { return gcnew cli::e_float(::ef::value_eps    ()); }
        static cli::e_float^ value_inf    (void) { return gcnew cli::e_float(::ef::value_inf    ()); }
        static cli::e_float^ value_nan    (void) { return gcnew cli::e_float(::ef::value_nan    ()); }
        static cli::e_float^ two          (void) { return gcnew cli::e_float(::ef::two          ()); }
        static cli::e_float^ three        (void) { return gcnew cli::e_float(::ef::three        ()); }
        static cli::e_float^ four         (void) { return gcnew cli::e_float(::ef::four         ()); }
        static cli::e_float^ five         (void) { return gcnew cli::e_float(::ef::five         ()); }
        static cli::e_float^ six          (void) { return gcnew cli::e_float(::ef::six          ()); }
        static cli::e_float^ seven        (void) { return gcnew cli::e_float(::ef::seven        ()); }
        static cli::e_float^ eight        (void) { return gcnew cli::e_float(::ef::eight        ()); }
        static cli::e_float^ nine         (void) { return gcnew cli::e_float(::ef::nine         ()); }
        static cli::e_float^ ten          (void) { return gcnew cli::e_float(::ef::ten          ()); }
        static cli::e_float^ twenty       (void) { return gcnew cli::e_float(::ef::twenty       ()); }
        static cli::e_float^ thirty       (void) { return gcnew cli::e_float(::ef::thirty       ()); }
        static cli::e_float^ forty        (void) { return gcnew cli::e_float(::ef::forty        ()); }
        static cli::e_float^ fifty        (void) { return gcnew cli::e_float(::ef::fifty        ()); }
        static cli::e_float^ hundred      (void) { return gcnew cli::e_float(::ef::hundred      ()); }
        static cli::e_float^ two_hundred  (void) { return gcnew cli::e_float(::ef::two_hundred  ()); }
        static cli::e_float^ three_hundred(void) { return gcnew cli::e_float(::ef::three_hundred()); }
        static cli::e_float^ four_hundred (void) { return gcnew cli::e_float(::ef::four_hundred ()); }
        static cli::e_float^ five_hundred (void) { return gcnew cli::e_float(::ef::five_hundred ()); }
        static cli::e_float^ thousand     (void) { return gcnew cli::e_float(::ef::thousand     ()); }
        static cli::e_float^ two_k        (void) { return gcnew cli::e_float(::ef::two_k        ()); }
        static cli::e_float^ three_k      (void) { return gcnew cli::e_float(::ef::three_k      ()); }
        static cli::e_float^ four_k       (void) { return gcnew cli::e_float(::ef::four_k       ()); }
        static cli::e_float^ five_k       (void) { return gcnew cli::e_float(::ef::five_k       ()); }
        static cli::e_float^ ten_k        (void) { return gcnew cli::e_float(::ef::ten_k        ()); }
        static cli::e_float^ twenty_k     (void) { return gcnew cli::e_float(::ef::twenty_k     ()); }
        static cli::e_float^ thirty_k     (void) { return gcnew cli::e_float(::ef::thirty_k     ()); }
        static cli::e_float^ forty_k      (void) { return gcnew cli::e_float(::ef::forty_k      ()); }
        static cli::e_float^ fifty_k      (void) { return gcnew cli::e_float(::ef::fifty_k      ()); }
        static cli::e_float^ hundred_k    (void) { return gcnew cli::e_float(::ef::hundred_k    ()); }
        static cli::e_float^ million      (void) { return gcnew cli::e_float(::ef::million      ()); }
        static cli::e_float^ ten_M        (void) { return gcnew cli::e_float(::ef::ten_M        ()); }
        static cli::e_float^ hundred_M    (void) { return gcnew cli::e_float(::ef::hundred_M    ()); }
        static cli::e_float^ billion      (void) { return gcnew cli::e_float(::ef::billion      ()); }
        static cli::e_float^ trillion     (void) { return gcnew cli::e_float(::ef::trillion     ()); }
        static cli::e_float^ googol       (void) { return gcnew cli::e_float(::ef::googol       ()); }
        static cli::e_float^ int32max     (void) { return gcnew cli::e_float(::ef::int32max     ()); }
        static cli::e_float^ int32min     (void) { return gcnew cli::e_float(::ef::int32min     ()); }
        static cli::e_float^ int64max     (void) { return gcnew cli::e_float(::ef::int64max     ()); }
        static cli::e_float^ int64min     (void) { return gcnew cli::e_float(::ef::int64min     ()); }
        static cli::e_float^ one_minus    (void) { return gcnew cli::e_float(::ef::one_minus    ()); }
        static cli::e_float^ tenth        (void) { return gcnew cli::e_float(::ef::tenth        ()); }
        static cli::e_float^ eighth       (void) { return gcnew cli::e_float(::ef::eighth       ()); }
        static cli::e_float^ fifth        (void) { return gcnew cli::e_float(::ef::fifth        ()); }
        static cli::e_float^ quarter      (void) { return gcnew cli::e_float(::ef::quarter      ()); }
        static cli::e_float^ third        (void) { return gcnew cli::e_float(::ef::third        ()); }
        static cli::e_float^ two_third    (void) { return gcnew cli::e_float(::ef::two_third    ()); }
        static cli::e_float^ four_third   (void) { return gcnew cli::e_float(::ef::four_third   ()); }
        static cli::e_float^ three_half   (void) { return gcnew cli::e_float(::ef::three_half   ()); }
        static cli::e_float^ sqrt2        (void) { return gcnew cli::e_float(::ef::sqrt2        ()); }
        static cli::e_float^ sqrt3        (void) { return gcnew cli::e_float(::ef::sqrt3        ()); }
        static cli::e_float^ pi           (void) { return gcnew cli::e_float(::ef::pi           ()); }
        static cli::e_float^ pi_half      (void) { return gcnew cli::e_float(::ef::pi_half      ()); }
        static cli::e_float^ pi_quarter   (void) { return gcnew cli::e_float(::ef::pi_quarter   ()); }
        static cli::e_float^ pi_squared   (void) { return gcnew cli::e_float(::ef::pi_squared   ()); }
        static cli::e_float^ two_pi       (void) { return gcnew cli::e_float(::ef::two_pi       ()); }
        static cli::e_float^ sqrt_pi      (void) { return gcnew cli::e_float(::ef::sqrt_pi      ()); }
        static cli::e_float^ degree       (void) { return gcnew cli::e_float(::ef::degree       ()); }
        static cli::e_float^ exp1         (void) { return gcnew cli::e_float(::ef::exp1         ()); }
        static cli::e_float^ ln2          (void) { return gcnew cli::e_float(::ef::ln2          ()); }
        static cli::e_float^ ln3          (void) { return gcnew cli::e_float(::ef::ln3          ()); }
        static cli::e_float^ ln10         (void) { return gcnew cli::e_float(::ef::ln10         ()); }
        static cli::e_float^ log10_2      (void) { return gcnew cli::e_float(::ef::log10_2      ()); }
        static cli::e_float^ golden_ratio (void) { return gcnew cli::e_float(::ef::golden_ratio ()); }
        static cli::e_float^ euler_gamma  (void) { return gcnew cli::e_float(::ef::euler_gamma  ()); }
        static cli::e_float^ catalan      (void) { return gcnew cli::e_float(::ef::catalan      ()); }
        static cli::e_float^ khinchin     (void) { return gcnew cli::e_float(::ef::khinchin     ()); }
        static cli::e_float^ glaisher     (void) { return gcnew cli::e_float(::ef::glaisher     ()); }

        // <functions/elementary/elementary_math.h>
        static INT32 max_iteration(void) { return ::ef::max_iteration(); }

        static cli::e_float^ floor(cli::e_float^ x) { return gcnew cli::e_float(::ef::floor(x)); }
        static cli::e_float^ ceil (cli::e_float^ x) { return gcnew cli::e_float(::ef::ceil (x)); }
        static INT32         sgn  (cli::e_float^ x) { return ::ef::sgn(x); }

        static bool isnan   (cli::e_float^ x)          { return ::ef::isnan(x); }
        static bool isnan   (cli::ef_complex^ z)       { return ::ef::isnan(z); }
        static bool isfinite(cli::e_float^ x)        { return ::ef::isfinite(x); }
        static bool isfinite(cli::ef_complex^ z)     { return ::ef::isfinite(z); }
        static bool isinf   (cli::e_float^ x)          { return ::ef::isinf(x); }
        static bool isinf   (cli::ef_complex^ z)       { return ::ef::isinf(z); }
        static bool isneg   (cli::e_float^ x)          { return ::ef::isneg(x); }
        static bool isneg   (cli::ef_complex^ z)       { return ::ef::isneg(z); }
        static bool ispos   (cli::e_float^ x)          { return ::ef::ispos(x); }
        static bool ispos   (cli::ef_complex^ z)       { return ::ef::ispos(z); }
        static bool isint   (cli::e_float^ x)          { return ::ef::isint(x); }
        static bool isint   (cli::ef_complex^ z)       { return ::ef::isint(z); }
        static bool isone   (cli::e_float^ x)          { return ::ef::isone(x); }
        static bool isone   (cli::ef_complex^ z)       { return ::ef::isone(z); }
        static bool iszero  (cli::e_float^ x)          { return ::ef::iszero(x); }
        static bool iszero  (cli::ef_complex^ z)       { return ::ef::iszero(z); }

        static INT64 tol(void) { return ::ef::tol(); }

        static cli::e_float^ fabs(cli::e_float^ x)   { return gcnew cli::e_float(::ef::fabs(x)); }
        static cli::e_float^ abs (cli::e_float^ x)   { return gcnew cli::e_float(::ef::fabs(x)); }
        static cli::e_float^ real(cli::e_float^ x)   { return x; }
        static cli::e_float^ imag(cli::e_float^ x)   { static_cast<void>(x); return cli::ef::zero(); }

        static cli::e_float^ integer_part(cli::e_float^ x) { return gcnew cli::e_float(::ef::integer_part(x)); }
        static cli::e_float^ decimal_part(cli::e_float^ x) { return gcnew cli::e_float(::ef::decimal_part(x)); }

        static void to_parts(cli::e_float^ x, double& mantissa, INT64& exponent) { ::ef::to_parts(x, mantissa, exponent); }

        static double to_double(cli::e_float^ x)    { return ::ef::to_double(x); }
        static double to_double(cli::ef_complex^ z) { return ::ef::to_double(z); }

        static INT64 to_int64(cli::e_float^ x)      { return ::ef::to_int64(x); }
        static INT64 to_int64(cli::ef_complex^ z)   { return ::ef::to_int64(z); }
        static INT32 to_int32(cli::e_float^ x)      { return ::ef::to_int32(x); }
        static INT32 to_int32(cli::ef_complex^ z)   { return ::ef::to_int32(z); }

        static bool small_arg(cli::e_float^ x)      { return ::ef::small_arg(x); }
        static bool small_arg(cli::ef_complex^ z)   { return ::ef::small_arg(z); }
        static bool large_arg(cli::e_float^ x)      { return ::ef::large_arg(x); }
        static bool large_arg(cli::ef_complex^ z)   { return ::ef::large_arg(z); }
        static bool near_one (cli::e_float^ x)      { return ::ef::near_one(x); }
        static bool near_one (cli::ef_complex^ z)   { return ::ef::near_one(z); }
        static bool near_int (cli::e_float^ x)      { return ::ef::near_int(x); }
        static bool near_int (cli::ef_complex^ z)   { return ::ef::near_int(z); }

        // <functions/elementary/elementary_trig.h>
        static void          sincos(cli::e_float^ x, cli::e_float^ s, cli::e_float^ c) { ::ef::sincos(x, const_cast< ::e_float*>(s->my_ptr()), const_cast< ::e_float*>(c->my_ptr())); }
        static cli::e_float^ sin   (cli::e_float^ x) { return gcnew cli::e_float(::ef::sin (x)); }
        static cli::e_float^ cos   (cli::e_float^ x) { return gcnew cli::e_float(::ef::cos (x)); }
        static cli::e_float^ tan   (cli::e_float^ x) { return gcnew cli::e_float(::ef::tan (x)); }
        static cli::e_float^ csc   (cli::e_float^ x) { return gcnew cli::e_float(::ef::csc (x)); }
        static cli::e_float^ sec   (cli::e_float^ x) { return gcnew cli::e_float(::ef::sec (x)); }
        static cli::e_float^ cot   (cli::e_float^ x) { return gcnew cli::e_float(::ef::cot (x)); }
        static cli::e_float^ asin  (cli::e_float^ x) { return gcnew cli::e_float(::ef::asin(x)); }
        static cli::e_float^ acos  (cli::e_float^ x) { return gcnew cli::e_float(::ef::acos(x)); }
        static cli::e_float^ atan  (cli::e_float^ x) { return gcnew cli::e_float(::ef::atan(x)); }
        static cli::e_float^ atan2 (cli::e_float^ y, cli::e_float^ x) { return gcnew cli::e_float(::ef::atan2(y, x)); }

        // <functions/elementary/elementary_trans.h>
        static cli::e_float^ pow2    (const INT64 p)                    { return gcnew cli::e_float(::ef::pow2(p)); }
        static cli::e_float^ pown    (cli::e_float^ x, const INT64 p)   { return gcnew cli::e_float(::ef::pown(x, p)); }
        static cli::e_float^ inv     (cli::e_float^ x)                  { return gcnew cli::e_float(::ef::inv(x)); }
        static cli::e_float^ sqrt    (cli::e_float^ x)                  { return gcnew cli::e_float(::ef::sqrt(x)); }
        static cli::e_float^ cbrt    (cli::e_float^ x)                  { return gcnew cli::e_float(::ef::cbrt(x)); }
        static cli::e_float^ rootn   (cli::e_float^ x, const INT32 p)   { return gcnew cli::e_float(::ef::rootn(x, p)); }
        static cli::e_float^ exp     (cli::e_float^ x)                  { return gcnew cli::e_float(::ef::exp(x)); }
        static cli::e_float^ log     (cli::e_float^ x)                  { return gcnew cli::e_float(::ef::log(x)); }
        static cli::e_float^ log10   (cli::e_float^ x)                  { return gcnew cli::e_float(::ef::log10(x)); }
        static cli::e_float^ loga    (cli::e_float^ a, cli::e_float^ x) { return gcnew cli::e_float(::ef::loga(a, x)); }
        static cli::e_float^ log1p   (cli::e_float^ x)                  { return gcnew cli::e_float(::ef::log1p(x)); }
        static cli::e_float^ log1p1m2(cli::e_float^ x)                  { return gcnew cli::e_float(::ef::log1p1m2(x)); }
        static cli::e_float^ pow     (cli::e_float^ x, cli::e_float^ a) { return gcnew cli::e_float(::ef::pow(x, a)); }
        static void          sinhcosh(cli::e_float^ x, cli::e_float^ s, cli::e_float^ c) { ::ef::sinhcosh(x, const_cast< ::e_float*>(s->my_ptr()), const_cast< ::e_float*>(c->my_ptr())); }
        static cli::e_float^ sinh    (cli::e_float^ x)                  { return gcnew cli::e_float(::ef::sinh(x)); }
        static cli::e_float^ cosh    (cli::e_float^ x)                  { return gcnew cli::e_float(::ef::cosh(x)); }
        static cli::e_float^ tanh    (cli::e_float^ x)                  { return gcnew cli::e_float(::ef::tanh(x)); }
        static cli::e_float^ asinh   (cli::e_float^ x)                  { return gcnew cli::e_float(::ef::asinh(x)); }
        static cli::e_float^ acosh   (cli::e_float^ x)                  { return gcnew cli::e_float(::ef::acosh(x)); }
        static cli::e_float^ atanh   (cli::e_float^ x)                  { return gcnew cli::e_float(::ef::atanh(x)); }

        // <functions/bessel/airy.h>
        static cli::e_float^ airy_a      (cli::e_float^ x) { return gcnew cli::e_float(::ef::airy_a(x)); }
        static cli::e_float^ airy_a_prime(cli::e_float^ x) { return gcnew cli::e_float(::ef::airy_a_prime(x)); }
        static cli::e_float^ airy_b      (cli::e_float^ x) { return gcnew cli::e_float(::ef::airy_b(x)); }
        static cli::e_float^ airy_b_prime(cli::e_float^ x) { return gcnew cli::e_float(::ef::airy_b_prime(x)); }

        static System::Collections::Generic::List<cli::e_float^>^ airy_a_zero(const UINT32 k);
        static System::Collections::Generic::List<cli::e_float^>^ airy_b_zero(const UINT32 k);

        // <functions/bessel/bessel.h>
        static cli::e_float^ cyl_bessel_i(cli::e_float^ v, cli::e_float^ x) { return gcnew cli::e_float(::ef::cyl_bessel_i(v, x)); }
        static cli::e_float^ cyl_bessel_i(const INT32 n, cli::e_float^ x)   { return gcnew cli::e_float(::ef::cyl_bessel_i(n, x)); }
        static cli::e_float^ cyl_bessel_j(cli::e_float^ v, cli::e_float^ x) { return gcnew cli::e_float(::ef::cyl_bessel_j(v, x)); }
        static cli::e_float^ cyl_bessel_j(const INT32 n, cli::e_float^ x)   { return gcnew cli::e_float(::ef::cyl_bessel_j(n, x)); }
        static cli::e_float^ cyl_bessel_k(cli::e_float^ v, cli::e_float^ x) { return gcnew cli::e_float(::ef::cyl_bessel_k(v, x)); }
        static cli::e_float^ cyl_bessel_k(const INT32 n, cli::e_float^ x)   { return gcnew cli::e_float(::ef::cyl_bessel_k(n, x)); }
        static cli::e_float^ cyl_bessel_y(cli::e_float^ v, cli::e_float^ x) { return gcnew cli::e_float(::ef::cyl_bessel_y(v, x)); }
        static cli::e_float^ cyl_bessel_y(const INT32 n, cli::e_float^ x)   { return gcnew cli::e_float(::ef::cyl_bessel_y(n, x)); }

        static cli::e_float^ cyl_bessel_i_prime(cli::e_float^ v, cli::e_float^ x) { return gcnew cli::e_float(::ef::cyl_bessel_i_prime(v, x)); }
        static cli::e_float^ cyl_bessel_i_prime(const INT32 n, cli::e_float^ x)   { return gcnew cli::e_float(::ef::cyl_bessel_i_prime(n, x)); }
        static cli::e_float^ cyl_bessel_j_prime(cli::e_float^ v, cli::e_float^ x) { return gcnew cli::e_float(::ef::cyl_bessel_j_prime(v, x)); }
        static cli::e_float^ cyl_bessel_j_prime(const INT32 n, cli::e_float^ x)   { return gcnew cli::e_float(::ef::cyl_bessel_j_prime(n, x)); }
        static cli::e_float^ cyl_bessel_k_prime(cli::e_float^ v, cli::e_float^ x) { return gcnew cli::e_float(::ef::cyl_bessel_k_prime(v, x)); }
        static cli::e_float^ cyl_bessel_k_prime(const INT32 n, cli::e_float^ x)   { return gcnew cli::e_float(::ef::cyl_bessel_k_prime(n, x)); }
        static cli::e_float^ cyl_bessel_y_prime(cli::e_float^ v, cli::e_float^ x) { return gcnew cli::e_float(::ef::cyl_bessel_y_prime(v, x)); }
        static cli::e_float^ cyl_bessel_y_prime(const INT32 n, cli::e_float^ x)   { return gcnew cli::e_float(::ef::cyl_bessel_y_prime(n, x)); }

        static System::Collections::Generic::List<cli::e_float^>^ cyl_bessel_j_zero(const INT32 n,   const UINT32 k);
        static System::Collections::Generic::List<cli::e_float^>^ cyl_bessel_j_zero(cli::e_float^ v, const UINT32 k);

        // <functions/elliptic/elliptic.h>
        static cli::e_float^ comp_ellint_1(cli::e_float^ m)                    { return gcnew cli::e_float(::ef::comp_ellint_1(m)); }
        static cli::e_float^      ellint_1(cli::e_float^ m, cli::e_float^ phi) { return gcnew cli::e_float(::ef::ellint_1(m, phi)); }
        static cli::e_float^ comp_ellint_2(cli::e_float^ m)                    { return gcnew cli::e_float(::ef::comp_ellint_2(m)); }
        static cli::e_float^      ellint_2(cli::e_float^ m, cli::e_float^ phi) { return gcnew cli::e_float(::ef::ellint_2(m, phi)); }

        // <functions/gamma/gamma.h>
        static cli::e_float^ gamma               (cli::e_float^ x)                                     { return gcnew cli::e_float(::ef::gamma(x)); } 
        static cli::e_float^ gamma_near_n        (const INT32 n, cli::e_float^ x)                      { return gcnew cli::e_float(::ef::gamma_near_n(n, x));}
        static cli::e_float^ incomplete_gamma    (cli::e_float^ a, cli::e_float^ x)                    { return gcnew cli::e_float(::ef::incomplete_gamma(a, x)); }
        static cli::e_float^ gen_incomplete_gamma(cli::e_float^ a, cli::e_float^ x0, cli::e_float^ x1) { return gcnew cli::e_float(::ef::gen_incomplete_gamma(a, x0, x1)); }
        static cli::e_float^ beta                (cli::e_float^ a, cli::e_float^ b)                    { return gcnew cli::e_float(::ef::beta(a, b)); }
        static cli::e_float^ incomplete_beta     (cli::e_float^ x, cli::e_float^ a, cli::e_float^ b)   { return gcnew cli::e_float(::ef::incomplete_beta(x, a, b)); }
        static cli::e_float^ factorial           (const UINT32 n)                                      { return gcnew cli::e_float(::ef::factorial(n)); }
        static cli::e_float^ factorial2          (const  INT32 n)                                      { return gcnew cli::e_float(::ef::factorial2(n)); }
        static cli::e_float^ binomial            (const UINT32 n, const UINT32 k)                      { return gcnew cli::e_float(::ef::binomial(n, k)); }
        static cli::e_float^ binomial            (const UINT32 n, cli::e_float^ y)                     { return gcnew cli::e_float(::ef::binomial(n, y)); }
        static cli::e_float^ binomial            (cli::e_float^ x, const UINT32 k)                     { return gcnew cli::e_float(::ef::binomial(x, k)); }
        static cli::e_float^ binomial            (cli::e_float^ x, cli::e_float^ y)                    { return gcnew cli::e_float(::ef::binomial(x, y)); }
        static cli::e_float^ pochhammer          (cli::e_float^ x, const UINT32 n)                     { return gcnew cli::e_float(::ef::pochhammer(x, n)); }
        static cli::e_float^ pochhammer          (cli::e_float^ x, cli::e_float^ a)                    { return gcnew cli::e_float(::ef::pochhammer(x, a)); }

        // <functions/gamma/poly_gamma.h>
        static cli::e_float^ poly_gamma(cli::e_float^ x)                { return gcnew cli::e_float(::ef::poly_gamma(x)); }
        static cli::e_float^ poly_gamma(const INT32 n, cli::e_float^ x) { return gcnew cli::e_float(::ef::poly_gamma(n, x)); }

        // <functions/hypergeometric/hermite.h>
        static cli::e_float^ hermite(cli::e_float^ v, cli::e_float^ x) { return gcnew cli::e_float(::ef::hermite(v, x)); }

        // <functions/hypergeometric/hypergeometric.h>
        static cli::e_float^ hyperg_0f0    (cli::e_float^ x)                                                    { return gcnew cli::e_float(::ef::hyperg_0f0(x)); }
        static cli::e_float^ hyperg_0f1    (cli::e_float^ b, cli::e_float^ x)                                   { return gcnew cli::e_float(::ef::hyperg_0f1(b, x)); }
        static cli::e_float^ hyperg_0f1_reg(cli::e_float^ a, cli::e_float^ x)                                   { return gcnew cli::e_float(::ef::hyperg_0f1_reg(a, x)); }
        static cli::e_float^ hyperg_1f0    (cli::e_float^ a, cli::e_float^ x)                                   { return gcnew cli::e_float(::ef::hyperg_1f0(a, x)); }
        static cli::e_float^ hyperg_1f1    (cli::e_float^ a, cli::e_float^ b, cli::e_float^ x)                  { return gcnew cli::e_float(::ef::hyperg_1f1(a, b, x)); }
        static cli::e_float^ hyperg_1f1_reg(cli::e_float^ a, cli::e_float^ b, cli::e_float^ x)                  { return gcnew cli::e_float(::ef::hyperg_1f1_reg(a, b, x)); }
        static cli::e_float^ hyperg_2f0    (cli::e_float^ a, cli::e_float^ b, cli::e_float^ x)                  { return gcnew cli::e_float(::ef::hyperg_2f0(a, b, x)); }
        static cli::e_float^ hyperg_2f1    (cli::e_float^ a, cli::e_float^ b, cli::e_float^ c, cli::e_float^ x) { return gcnew cli::e_float(::ef::hyperg_2f1(a, b, c, x)); }
        static cli::e_float^ hyperg_2f1_reg(cli::e_float^ a, cli::e_float^ b, cli::e_float^ c, cli::e_float^ x) { return gcnew cli::e_float(::ef::hyperg_2f1_reg(a, b, c, x)); }
        static cli::e_float^ hyperg_pfq    (System::Collections::Generic::List<cli::e_float^>^ a, System::Collections::Generic::List<cli::e_float^>^ b, cli::e_float^ x);
        static cli::e_float^ conf_hyperg   (cli::e_float^ a, cli::e_float^ c, cli::e_float^ x)                  { return cli::ef::hyperg_1f1(a, c, x); }
        static cli::e_float^      hyperg   (cli::e_float^ a, cli::e_float^ b, cli::e_float^ c, cli::e_float^ x) { return cli::ef::hyperg_2f1(a, b, c, x); }

        // <functions/hypergeometric/laguerre.h>
        static cli::e_float^ laguerre(cli::e_float^ v, cli::e_float^ x)                  { return gcnew cli::e_float(::ef::laguerre(v, x)); }
        static cli::e_float^ laguerre(cli::e_float^ v, cli::e_float^ L, cli::e_float^ x) { return gcnew cli::e_float(::ef::laguerre(v, L, x)); }
        static cli::e_float^ laguerre(const INT32 n,   cli::e_float^ L, cli::e_float^ x) { return gcnew cli::e_float(::ef::laguerre(n, L, x)); }
        static cli::e_float^ laguerre(const INT32 n,   const INT32 m,   cli::e_float^ x) { return gcnew cli::e_float(::ef::laguerre(n, m, x)); }

        // <functions/hypergeometric/legendre.h>
        static cli::e_float^ legendre_p(cli::e_float^ v, cli::e_float^ x)                  { return gcnew cli::e_float(::ef::legendre_p(v, x)); }
        static cli::e_float^ legendre_p(cli::e_float^ v, cli::e_float^ u, cli::e_float^ x) { return gcnew cli::e_float(::ef::legendre_p(v, u, x)); }
        static cli::e_float^ legendre_p(cli::e_float^ v, const INT32   m, cli::e_float^ x) { return gcnew cli::e_float(::ef::legendre_p(v, m, x)); }
        static cli::e_float^ legendre_p(const INT32   n, cli::e_float^ u, cli::e_float^ x) { return gcnew cli::e_float(::ef::legendre_p(n, u, x)); }
        static cli::e_float^ legendre_p(const INT32   n, const INT32   m, cli::e_float^ x) { return gcnew cli::e_float(::ef::legendre_p(n, m, x)); }
        static cli::e_float^ legendre_q(cli::e_float^ v, cli::e_float^ x)                  { return gcnew cli::e_float(::ef::legendre_q(v, x)); }
        static cli::e_float^ legendre_q(cli::e_float^ v, cli::e_float^ u, cli::e_float^ x) { return gcnew cli::e_float(::ef::legendre_q(v, u, x)); }
        static cli::e_float^ legendre_q(cli::e_float^ v, const INT32   m, cli::e_float^ x) { return gcnew cli::e_float(::ef::legendre_q(v, m, x)); }
        static cli::e_float^ legendre_q(const INT32   n, cli::e_float^ u, cli::e_float^ x) { return gcnew cli::e_float(::ef::legendre_q(n, u, x)); }
        static cli::e_float^ legendre_q(const INT32   n, const INT32   m, cli::e_float^ x) { return gcnew cli::e_float(::ef::legendre_q(n, m, x)); }

        // <functions/hypergeometric/parabolic_cylinder_d.h>
        static cli::e_float^ weber_d(cli::e_float^ v, cli::e_float^ x) { return gcnew cli::e_float(::ef::weber_d(v, x)); }

        // <functions/integer/integer.h>
        static cli::e_float^ euler    (const UINT32 n)                 { return gcnew cli::e_float(::ef::euler(n)); }
        static cli::e_float^ bernoulli(const UINT32 n)                 { return gcnew cli::e_float(::ef::bernoulli(n)); }
        static cli::e_float^ stirling2(const UINT32 n, const UINT32 k) { return gcnew cli::e_float(::ef::stirling2(n, k)); }

        // <functions/integer/prime.h>
        System::Collections::Generic::List<System::UInt32>^ prime(const UINT32 n);

        // <functions/polynomials/polynomials.h>
        static cli::e_float^ chebyshev_t(const INT32 n, cli::e_float^ x) { return gcnew cli::e_float(::ef::chebyshev_t(n, x)); }
        static cli::e_float^ chebyshev_u(const INT32 n, cli::e_float^ x) { return gcnew cli::e_float(::ef::chebyshev_u(n, x)); }
        static cli::e_float^ hermite    (const INT32 n, cli::e_float^ x) { return gcnew cli::e_float(::ef::hermite    (n, x)); }
        static cli::e_float^ laguerre   (const INT32 n, cli::e_float^ x) { return gcnew cli::e_float(::ef::laguerre   (n, x)); }
        static cli::e_float^ legendre_p (const INT32 n, cli::e_float^ x) { return gcnew cli::e_float(::ef::legendre_p (n, x)); }
        static cli::e_float^ legendre_q (const INT32 n, cli::e_float^ x) { return gcnew cli::e_float(::ef::legendre_q (n, x)); }

        static cli::e_float^ chebyshev_t(const UINT32 n, cli::e_float^ x, System::Collections::Generic::List<cli::e_float^>^ vp);
        static cli::e_float^ chebyshev_u(const UINT32 n, cli::e_float^ x, System::Collections::Generic::List<cli::e_float^>^ vp);
        static cli::e_float^ hermite    (const UINT32 n, cli::e_float^ x, System::Collections::Generic::List<cli::e_float^>^ vp);
        static cli::e_float^ laguerre   (const UINT32 n, cli::e_float^ x, System::Collections::Generic::List<cli::e_float^>^ vp);

        // <functions/zeta/polylog.h>
        static cli::e_float^ poly_logarithm(const INT32 n, cli::e_float^ x) { return gcnew cli::e_float(::ef::poly_logarithm(n, x)); }

        // <functions/zeta/zeta.h>
        static cli::e_float^ riemann_zeta(const INT32 n)                    { return gcnew cli::e_float(::ef::riemann_zeta(n)); }
        static cli::e_float^ riemann_zeta(cli::e_float^ s)                  { return gcnew cli::e_float(::ef::riemann_zeta(s)); }
        static cli::e_float^ hurwitz_zeta(cli::e_float^ s, const INT32 n)   { return gcnew cli::e_float(::ef::hurwitz_zeta(s, n)); }
        static cli::e_float^ hurwitz_zeta(cli::e_float^ s, cli::e_float^ a) { return gcnew cli::e_float(::ef::hurwitz_zeta(s, a)); }
      };
    }
  }

#endif // _FUNCTIONS_CLR_2010_01_06_H_
