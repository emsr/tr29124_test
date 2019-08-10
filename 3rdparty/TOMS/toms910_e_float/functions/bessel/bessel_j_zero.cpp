
#include <numeric>

#include <functions/bessel/bessel.h>
#include <functions/bessel/airy.h>
#include <utility/util_coefficient_expansion.h>
#include <utility/util_find_root_bisect.h>
#include <utility/util_find_root_newton_raphson.h>
#include <utility/util_interpolate.h>

namespace BesselJZero
{
  static const double& two_thirds(void)
  {
    static const double two_over_three = static_cast<double>(2.0) / static_cast<double>(3.0);
    return two_over_three;
  }

  class RootZetaZ : public Util::FindRootBisect<double>
  {
  private:

    static const double& my_tol(void) { static const double val = static_cast<double>(0.001);   return val; }

    const double zeta;

    const RootZetaZ& operator=(const RootZetaZ&);

  public:

    RootZetaZ(const double& var_min,
              const double& var_max,
              const double& var_tol,
              const double& var_zeta) : Util::FindRootBisect<double>(var_min, var_max, var_tol),
                                        zeta(var_zeta) { }
  private:

    virtual double my_function(const double& x) const
    {
      return   ::sqrt((x * x) - static_cast<double>(1.0))
             - ::acos(static_cast<double>(1.0) / x)
             - (two_thirds() * ::pow(zeta, static_cast<double>(1.5)));
    }
  };

  class NewtonRaphsonBesselJv : public Util::FindRootNewtonRaphson<e_float>
  {
  private:

    const e_float my_v;

  public:
  
    NewtonRaphsonBesselJv(const e_float& lo,
                          const e_float& hi,
                          const e_float& tol,
                          const e_float& v) : Util::FindRootNewtonRaphson<e_float>(lo, hi, tol),
                                              my_v(v) { }

    virtual ~NewtonRaphsonBesselJv() { }

  private:

    virtual e_float my_function  (const e_float& x) const { return ef::cyl_bessel_j(my_v, x); }
    virtual e_float my_derivative(const e_float& x) const { return ef::cyl_bessel_j_prime(my_v, x); }
  };

  class NewtonRaphsonBesselJn : public Util::FindRootNewtonRaphson<e_float>
  {
  private:

    const INT32 my_n;

  public:
  
    NewtonRaphsonBesselJn(const e_float& lo,
                          const e_float& hi,
                          const e_float& tol,
                          const INT32    n) : Util::FindRootNewtonRaphson<e_float>(lo, hi, tol),
                                              my_n(n) { }

    virtual ~NewtonRaphsonBesselJn() { }

  private:

    virtual e_float my_function  (const e_float& x) const { return ef::cyl_bessel_j(my_n, x); }
    virtual e_float my_derivative(const e_float& x) const { return ef::cyl_bessel_j_prime(my_n, x); }
  };

  class Quadratic : public Util::FunctionOperation<double>
  {
  private:
  
    Quadratic(const Quadratic&);
    const Quadratic& operator=(const Quadratic&);

    const   double a;
    const   double b;
    const   double c;
    mutable double rm;
    mutable double rp;

  public:

    Quadratic(const double& A, const double& B, const double& C) : a (A),
                                                                   b (B),
                                                                   c (C),
                                                                   rm(static_cast<double>(0.0)),
                                                                   rp(static_cast<double>(0.0)) { }

    virtual ~Quadratic() { }

  private:

    virtual double my_operation(void) const
    {
      const double sq = sqrt((b * b) - (static_cast<double>(4.0) * (a * c)));

      rm = (-b - sq) / (static_cast<double>(2.0) * a);
      rp = (-b + sq) / (static_cast<double>(2.0) * a);

      op_ok = ef::isfinite(rm) || ef::isfinite(rp);

      return static_cast<double>(0.0);
    }

    virtual double my_function(const double& x) const { return ((a * a) * (x * x)) + ((b * x) + c); } // NOCOVER_LINE

  public:

    const double& root_pos(void) const { return rp; }
    const double& root_neg(void) const { return rm; }
  };

  static double first_zero_AS_Eq_9_5_14(const double nu)
  {
    // Obtain the estimate of the first zero of Jnu. The estimate is calculated from
    // Abramowitz and Stegun 9.5.14, page 371.
    
    const double one_over_nu_pow_two_third = ::pow(nu, -two_thirds());

    static const std::tr1::array<double, 6u> AS_Eq_9_5_14_coef =
    {{
      static_cast<double>( 1.0),
      static_cast<double>( 1.85575),
      static_cast<double>( 1.03315),
      static_cast<double>(-0.00397),
      static_cast<double>(-0.0908),
      static_cast<double>( 0.043)
    }};

    return std::accumulate(AS_Eq_9_5_14_coef.begin(),
                           AS_Eq_9_5_14_coef.end(),
                           static_cast<double>(0.0),
                           Util::coefficient_expansion<double, double>(one_over_nu_pow_two_third, nu));
  }

  static double sth_zero_large_nu_AS_Eq_9_5_26(const double nu, const UINT32 s)
  {
    // Obtain the estimate of the s'th zero of Jnu for nu larger than about one.
    // The estimate is computed from Abramowitz and Stegun 9.5.22 and 9.5.26, page 371.
    
    // The inversion of z as a function of zeta, as described in the text following
    // A&S Eq. 9.5.26, is accomplished by performing a Taylor expansion of Eq. 9.3.39
    // for large z and solving the resulting quadratic equation taking the positive root
    // of the quadratic.
    // In other words: (2/3)(-zeta)^(3/2) approx = z + 1/2z - pi/2
    // which leads to: z^2 - [(2/3)(-zeta)^(3/2) + pi/2]z + 1/2 = 0,
    // The positive root of the quadratic is used. Subsequently the bisection method is
    // used to find a refined, more exact value for the root z as a function of zeta.

    static const double pi_half = static_cast<double>(1.57079632679489661923);
    static const double one     = static_cast<double>(1.0);
    static const double half    = static_cast<double>(0.5);

    const double zeta = ::pow(nu, -two_thirds()) * (-AiryZero::ai_estimate_sth_zero(s));

    const double b = -((two_thirds() * ::pow(zeta, static_cast<double>(1.5))) + pi_half);

    const BesselJZero::Quadratic qz(one, b, half);
    qz.operation();

    const double z0   = qz.root_pos();
    const double zmin = z0 - one;

    // Dynamically adjust the tolerance of the root-finding based on the order.
    // Large orders require a tighter tolerance.

    static const double tol_max = static_cast<double>(1.0e-3);
    static const double tol_min = static_cast<double>(1.0e-9);
    
    double tol = one / (nu * static_cast<double>(10.0));
    
    if(tol < tol_min) { tol = tol_min; }
    if(tol > tol_max) { tol = tol_max; }

    RootZetaZ rzz(std::max(one, zmin), z0 + one, tol, zeta);

    // Calculate the root z as a function of zeta.
    const double z = rzz.operation();

    const double zsq_minus_one      = (z * z) - one;
    const double zsq_minus_one_sqrt = ::sqrt(zsq_minus_one);

    const double h = ::pow((static_cast<double>(4.0) * zeta) / zsq_minus_one, static_cast<double>(0.25));

    const double b0_term_5_24 = static_cast<double>(5.0) / (static_cast<double>(24.0) * (zsq_minus_one * zsq_minus_one_sqrt));
    const double b0_term_1_8  =                      one / (static_cast<double>( 8.0) *  zsq_minus_one_sqrt);
    const double b0_term_5_48 = static_cast<double>(5.0) / (static_cast<double>(48.0) * (zeta * zeta));

    const double b0 = -b0_term_5_48 + ((b0_term_5_24 + b0_term_1_8) / ::sqrt(zeta));

    const double f1 = ((z * (h * h)) * b0) / static_cast<double>(2.0);

    return (nu * z) + (f1 / nu);
  }
  
  static double sth_zero_small_nu_AS_Eq_9_5_12(const double nu, const UINT32 s)
  {
    // Obtain an estimate of the s'th zero of Jnu for nu smaller than about one.
    // The estimate is computed from Abramowitz and Stegun 9.5.12, page 371.

    const double mu           = (nu * nu) * static_cast<double>(4.0);
    const double mu_squared   = mu * mu;
    const double mu_cubed     = mu * mu_squared;
    const double mu_minus_one = mu - static_cast<double>(1.0);

    const std::tr1::array<double, 4u> AS_Eq_9_5_12_coef =
    {{
      -mu_minus_one,
      -((static_cast<double>(4.0) * mu_minus_one) * (( static_cast<double>(7.0) * mu) - static_cast<double>(31.0))) / static_cast<double>(3.0),
      -((static_cast<double>(32.0) * mu_minus_one) * (((static_cast<double>(83.0) * mu_squared) - (static_cast<double>(982.0) * mu)) + static_cast<double>(3779.0))) / static_cast<double>(15.0),
      -((static_cast<double>(64.0) * mu_minus_one) * (((static_cast<double>(6949.0) * mu_cubed) - (static_cast<double>(153855.0) * mu_squared)) + ((static_cast<double>(1585743.0) * mu) - static_cast<double>(6277237.0)))) / static_cast<double>(105.0)
    }};

    const double beta                        = static_cast<double>(3.14159265358979323) * ((static_cast<double>(s + 1u) + (nu / static_cast<double>(2.0))) - static_cast<double>(0.25));
    const double one_over_eight_beta         = static_cast<double>(1.0) / (static_cast<double>(8.0) * beta);
    const double one_over_eight_beta_squared = one_over_eight_beta * one_over_eight_beta;

    return std::accumulate(AS_Eq_9_5_12_coef.begin(),
                           AS_Eq_9_5_12_coef.end(),
                           beta,
                           Util::coefficient_expansion<double, double>(one_over_eight_beta_squared, one_over_eight_beta));
  }
}

std::deque<e_float> ef::cyl_bessel_j_zero(const e_float& v, const UINT32 k)
{
  std::deque<e_float> zeros;

  if(ef::isneg(v))
  {
    return zeros;
  }

  const double nu = ef::to_double(v);

  // Interpolation data for the value of the first zero of Jnu for 0 < nu < 1.
  static const std::tr1::array<Util::point<double>, 6u> initial_guess_first_zero_data =
  {{
    Util::point<double>(static_cast<double>(0.0), static_cast<double>(2.4048)),
    Util::point<double>(static_cast<double>(0.2), static_cast<double>(2.7071)),
    Util::point<double>(static_cast<double>(0.4), static_cast<double>(2.9988)),
    Util::point<double>(static_cast<double>(0.6), static_cast<double>(3.2825)),
    Util::point<double>(static_cast<double>(0.8), static_cast<double>(3.5598)),
    Util::point<double>(static_cast<double>(1.0), static_cast<double>(3.8317)),
  }};

  static const std::vector<Util::point<double> > initial_guess_first_zero_points(initial_guess_first_zero_data.begin(),
                                                                                 initial_guess_first_zero_data.end());

  // Obtain an estimate for the first zero of Jnu. Interpolation is used for values
  // of nu ranging from 0 < nu < 1. Abramowitz and Stegun Eq. 9.5.14 is used for
  // values of nu for which nu > 1.

  double initial_guess = nu < static_cast<double>(1.2) ? Util::linear_interpolate<double>::interpolate(nu, initial_guess_first_zero_points)
                                                       : BesselJZero::first_zero_AS_Eq_9_5_14(nu);

  for(UINT32 i = static_cast<UINT32>(0u); i < k; i++)
  {
    // Compute the Newton-Raphson iteration.
    const BesselJZero::NewtonRaphsonBesselJv nr_jv(e_float(initial_guess - static_cast<double>(0.5)),
                                                   e_float(initial_guess + static_cast<double>(0.5)),
                                                   std::numeric_limits<e_float>::epsilon(),
                                                   v);

    const e_float root = nr_jv.operation();

    if(nr_jv.success())
    {
      zeros.push_back(root);
    }
    else
    {
      break;
    }

    initial_guess = nu < static_cast<double>(1.2) ? BesselJZero::sth_zero_small_nu_AS_Eq_9_5_12(nu, static_cast<UINT32>(i + 1u))
                                                  : BesselJZero::sth_zero_large_nu_AS_Eq_9_5_26(nu, static_cast<UINT32>(i + 1u));
  }

  return zeros;
}

std::deque<e_float> ef::cyl_bessel_j_zero(const INT32 n, const UINT32 k)
{
  std::deque<e_float> zeros;

  if(n < static_cast<INT32>(0))
  {
    return zeros;
  }

  static const std::tr1::array<double, 6u> first_few_zeros_data =
  {{
    static_cast<double>(2.4048),
    static_cast<double>(3.8317),
    static_cast<double>(5.1356),
    static_cast<double>(6.3802),
    static_cast<double>(7.5883),
    static_cast<double>(8.7715)
  }};

  const double nu = static_cast<double>(n);

  // Obtain an estimate for the first zero of Jn. Tabulated values are used for
  // small n. Abramowitz and Stegun Equation 9.5.14 is used for larger values of n.

  double initial_guess = static_cast<std::size_t>(n) < first_few_zeros_data.size() ? first_few_zeros_data[static_cast<std::size_t>(n)]
                                                                                   : BesselJZero::first_zero_AS_Eq_9_5_14(nu);

  for(UINT32 i = static_cast<UINT32>(0u); i < k; i++)
  {
    // Compute the Newton-Raphson iteration.

    const BesselJZero::NewtonRaphsonBesselJn nr_jn(e_float(initial_guess - static_cast<double>(0.5)),
                                                   e_float(initial_guess + static_cast<double>(0.5)),
                                                   std::numeric_limits<e_float>::epsilon(),
                                                   n);

    const e_float root = nr_jn.operation();

    if(nr_jn.success())
    {
      zeros.push_back(root);
    }
    else
    {
      break;
    }

    initial_guess = n <= static_cast<INT32>(1) ? BesselJZero::sth_zero_small_nu_AS_Eq_9_5_12(nu, static_cast<UINT32>(i + 1u))
                                               : BesselJZero::sth_zero_large_nu_AS_Eq_9_5_26(nu, static_cast<UINT32>(i + 1u));
  }

  return zeros;
}
