
#include <vector>
#include <numeric>

#include <e_float/e_float.h>
#include <functions/bessel/airy.h>
#include <functions/constants/constants.h>
#include <functions/elementary/elementary.h>
#include <utility/util_coefficient_expansion.h>
#include <utility/util_find_root_newton_raphson.h>

namespace AiryZero
{
  struct compute
  {
  private:
  
    compute(const compute&);
    const compute& operator=(const compute&);

  public:

    compute() { }
    virtual ~compute() { }

  public:

    std::deque<e_float> compute_zeros(const UINT32 k) const;
    double fz(const std::size_t i) const;

  private:

    virtual double zsd(void) const = 0;
    virtual Util::FindRootNewtonRaphson<e_float>* newton_raphson_generator(const e_float& lo,
                                                                           const e_float& hi,
                                                                           const e_float& tol) const = 0;
    virtual const std::vector<double>& AS_Table_10_13_data(void) const = 0;
    virtual const std::vector<double>& AS_Eq_10_4_105_coef(void) const = 0;
  };
  
  struct Ai : public compute
  {
  private:

    class NewtonRaphsonAiryAi : public Util::FindRootNewtonRaphson<e_float>
    {
    public:
    
      NewtonRaphsonAiryAi(const e_float& lo,
                          const e_float& hi,
                          const e_float& tol) : Util::FindRootNewtonRaphson<e_float>(lo, hi, tol) { }

      virtual ~NewtonRaphsonAiryAi() { }

    private:

      virtual e_float my_function  (const e_float& x) const { return ef::airy_a(x); }
      virtual e_float my_derivative(const e_float& x) const { return ef::airy_a_prime(x); }
    };

  public:

    Ai() { }
    virtual ~Ai() { }

  private:

    virtual double zsd(void) const { return static_cast<double>(1.0); }

    virtual Util::FindRootNewtonRaphson<e_float>* newton_raphson_generator(const e_float& lo,
                                                                           const e_float& hi,
                                                                           const e_float& tol) const
    {
      return static_cast<Util::FindRootNewtonRaphson<e_float>*>(new NewtonRaphsonAiryAi(lo, hi, tol));
    }

    virtual const std::vector<double>& AS_Table_10_13_data(void) const
    {
      static const std::tr1::array<double, 10u> a =
      {{
        static_cast<double>( -2.33810741),
        static_cast<double>( -4.08794944),
        static_cast<double>( -5.52055983),
        static_cast<double>( -6.78670809),
        static_cast<double>( -7.94413359),
        static_cast<double>( -9.02265085),
        static_cast<double>(-10.04017434),
        static_cast<double>(-11.00852430),
        static_cast<double>(-11.93601556),
        static_cast<double>(-12.82877675)
      }};

      static const std::vector<double> v(a.begin(), a.end());

      return v;
    }
    
    virtual const std::vector<double>& AS_Eq_10_4_105_coef(void) const
    {
      static const std::tr1::array<double, 6u> a =
      {{
        static_cast<double>(           1.0),
        static_cast<double>(           5.0) / static_cast<double>(       48.0),
        static_cast<double>(          -5.0) / static_cast<double>(       36.0),
        static_cast<double>(       77125.0) / static_cast<double>(    82944.0),
        static_cast<double>(  -108056875.0) / static_cast<double>(  6967296.0),
        static_cast<double>(162375596875.0) / static_cast<double>(334430208.0)
      }};

      static const std::vector<double> v(a.begin(), a.end());

      return v;
    }
  };

  struct Bi : public compute
  {
  private:

    class NewtonRaphsonAiryBi : public Util::FindRootNewtonRaphson<e_float>
    {
    public:
    
      NewtonRaphsonAiryBi(const e_float& lo,
                            const e_float& hi,
                            const e_float& tol) : Util::FindRootNewtonRaphson<e_float>(lo, hi, tol) { }

      virtual ~NewtonRaphsonAiryBi() { }

    private:

      virtual e_float my_function  (const e_float& x) const { return ef::airy_b(x); }
      virtual e_float my_derivative(const e_float& x) const { return ef::airy_b_prime(x); }
    };

  public:

    Bi() { }
    virtual ~Bi() { }

  private:

    virtual double zsd(void) const { return static_cast<double>(3.0); }

    virtual Util::FindRootNewtonRaphson<e_float>* newton_raphson_generator(const e_float& lo,
                                                                           const e_float& hi,
                                                                           const e_float& tol) const
    {
      return static_cast<Util::FindRootNewtonRaphson<e_float>*>(new NewtonRaphsonAiryBi(lo, hi, tol));
    }

    virtual const std::vector<double>& AS_Table_10_13_data(void) const
    {
      static const std::tr1::array<double, 10u> a =
      {{
        static_cast<double>( -1.17371322),
        static_cast<double>( -3.27109330),
        static_cast<double>( -4.83073784),
        static_cast<double>( -6.16985213),
        static_cast<double>( -7.37676208),
        static_cast<double>( -8.49194885),
        static_cast<double>( -9.53819438),
        static_cast<double>(-10.52991351),
        static_cast<double>(-11.47695335),
        static_cast<double>(-12.38641714)
      }};

      static const std::vector<double> v(a.begin(), a.end());

      return v;
    }
    
    virtual const std::vector<double>& AS_Eq_10_4_105_coef(void) const
    {
      static const std::tr1::array<double, 6u> a =
      {{
        static_cast<double>(           1.0),
        static_cast<double>(          -7.0) / static_cast<double>(       48.0),
        static_cast<double>(          35.0) / static_cast<double>(      288.0),
        static_cast<double>(     -181223.0) / static_cast<double>(   207360.0),
        static_cast<double>(    18683371.0) / static_cast<double>(  1244160.0),
        static_cast<double>(-91145884361.0) / static_cast<double>(191102976.0)
      }};

      static const std::vector<double> v(a.begin(), a.end());

      return v;
    }
  };
}

double AiryZero::ai_estimate_sth_zero(const UINT32 s)
{
  // Returns the double estimate of the s'th negative zero of of Ai.
  return AiryZero::Ai().fz(static_cast<std::size_t>(s));
}

double AiryZero::compute::fz(const std::size_t i) const
{
  // Obtain the estimate of the first few zeros from the tabulated data in
  // Abramowitz and Stegun Table 10.13, page 478.

  if(i < AS_Table_10_13_data().size())
  {
    return AS_Table_10_13_data()[i];
  }
  else
  {
    // Implements the function f(z) as shown in Abramowitz and Stegun Eq. 10.4.105, page 450.

    static const double three_pi = static_cast<double>(9.424777960769379715);

    // Compute f(z) from Abramowitz and Stegun Eq. 10.4.94 - 10.4.101, page 450.
    // Note that the call to the function zsd(...) uses the virtual mechanism.

    const double s = static_cast<double>(static_cast<UINT32>(i + 1u));
    const double z = (((static_cast<double>(4.0) * s) - zsd()) * three_pi) / static_cast<double>(8.0);
    const double one_over_z_squared = static_cast<double>(1.0) / (z * z);

    return  -::pow(z, static_cast<double>(2.0) / static_cast<double>(3.0))
           * std::accumulate(AS_Eq_10_4_105_coef().begin(),
                             AS_Eq_10_4_105_coef().end(),
                             static_cast<double>(0.0),
                             Util::coefficient_expansion<double, double>(one_over_z_squared));
  }
}

std::deque<e_float> AiryZero::compute::compute_zeros(const UINT32 k) const
{
  // Compute a series of Airy function zeros using Newton-Raphson iteration.

  std::deque<e_float> zeros;

  double range_for_root_search = static_cast<double>(0.1);

  for(std::size_t i = static_cast<std::size_t>(0u); i < static_cast<std::size_t>(k); i++)
  {
    const double initial_estimate = fz(i);

    // Compute the Newton-Raphson iteration. Note that the calls to the functions
    // f(...) and df(...), which respectively return a function and its derivative,
    // do make use of the virtual mechanism.

    const Util::FindRootNewtonRaphson<e_float>* p_nr = newton_raphson_generator(e_float(initial_estimate - static_cast<double>(0.5)),
                                                                                e_float(initial_estimate + static_cast<double>(0.5)),
                                                                                std::numeric_limits<e_float>::epsilon());

    const e_float root      = p_nr->operation();
    const bool    b_success = p_nr->success();

    delete p_nr;

    if(b_success)
    {
      zeros.push_back(root);
    }
    else
    {
      break;
    }

    if(zeros.size() > static_cast<std::deque<e_float>::size_type>(1u))
    {
      // Dynamically adjust the search range for the next root based on the distance
      // between the previous two roots.
      const e_float distance_between_roots = ef::fabs(zeros.back() - *(zeros.end() - 2u));

      range_for_root_search = static_cast<double>(0.5) * ef::to_double(distance_between_roots);
    }
  }

  return zeros;
}

std::deque<e_float> ef::airy_a_zero(const UINT32 k) { return AiryZero::Ai().compute_zeros(k); }
std::deque<e_float> ef::airy_b_zero(const UINT32 k) { return AiryZero::Bi().compute_zeros(k); }
