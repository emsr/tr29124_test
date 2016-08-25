
namespace std
{
  // Use narrow structs for  return types.
  // Prefer returns to pointer or ref arguments.
  // This will work swimmingly with structured bindings.
  //
  // I think all the basic math functions should be constexpr.
  // Error haqndling - the old ones have global error errno (from C).
  // The specfuns throw.
  //   - should these and the old basic ones have throwing versions
  //     - since we don't like global error reporting
  //     - exception is part of the signature(?)
  //     - only if peoplr can flip back to another or no error reporting...
  // We could do like filesystem and some others and double the api with
  // error_code return args.  In math, errors really are exceptional where
  // libs with this return failure is quite often an option.  In math failure
  // needs to be figured out, fixed, cleaned up after.

  // Log to arbitrary base - the inverse of pow(base, x).
  float logf(float base, float x);
  double log(double base, double x);
  long double logl(long double base, long double x);


  // Trigonometric functions

  // Combined sine and cosine.
  template<typename Tp>
    struct sincos_t
    {
      Tp sin_value;
      Tp cos_value;
    };

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
    }

  // Combined reperiodized sine and cosine.
  sincos_t<float> sincos_pif(float x);
  sincos_t<double> sincos_pi(double x);
  sincos_t<long double> sincos_pil(long double x);

  sincos_t<float> sincos_pif(reperiod_t<float x);
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

  // This is essentially a poor man's complex.
  // Conversion?
  template<typename Tp>
    lgamma_t
    {
      Tp lgamma_value;
      Tp sign;
    };

  lgamma_t<float> slgammaf(float x);
  lgamma_t<double> slgamma(double x);
  lgamma_t<long double> slgammal(long double x);

} // namespace std
