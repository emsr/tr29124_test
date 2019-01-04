
#include <functional>
#include <numeric>

#include <functions/bessel/bessel_asymp.h>
#include <functions/bessel/bessel_recursion.h>
#include <functions/tables/tables.h>
#include <utility/util_interpolate.h>

namespace BesselAsymp
{
  // Define function objects types for the zeta-Debye product sums in
  // the asymptotic Airy-type Bessel expansions.

  struct multiplies_zeta_debye : public std::multiplies<e_float>
  {
  private:

    const multiplies_zeta_debye& operator=(const multiplies_zeta_debye&);

  public:

    virtual ~multiplies_zeta_debye() { }

    const bool  has_phase;
          INT32 s_phase;

    multiplies_zeta_debye(const multiplies_zeta_debye& mzd) : has_phase(mzd.has_phase),
                                                              s_phase  (mzd.s_phase) { }

    multiplies_zeta_debye(const bool hp) : has_phase(hp),
                                           s_phase  (static_cast<INT32>(0)) { }

    e_float operator()(const e_float& u, const e_float& v)
    {
      return !neg_phase() ? (u * v) : -(u * v);
    }

  private:

    virtual INT32 my_sp(void) const { return s_phase; }

    bool neg_phase(void)
    {
      const bool phase_is_neg = static_cast<INT32>(static_cast<INT32>(my_sp() / static_cast<INT32>(2)) % static_cast<INT32>(2)) != static_cast<INT32>(0);

      ++s_phase;

      return has_phase && phase_is_neg;
    }
  };

  struct multiplies_zeta_debye_mu : public multiplies_zeta_debye
  {
    multiplies_zeta_debye_mu(const bool hp) : multiplies_zeta_debye(hp) { }
  };

  struct multiplies_zeta_debye_lambda : public multiplies_zeta_debye
  {
  public:

    multiplies_zeta_debye_lambda(const bool hp) : multiplies_zeta_debye(hp) { }

  private:

    virtual INT32 my_sp(void) const { return static_cast<INT32>(s_phase + static_cast<INT32>(3)); }
  };

  e_float UniformAsyBase::UniformAsymptotic(const e_float& v, const e_float& x) const
  {
    // Abramowitz and Stegun 9.3.35 - 9.3.42, page 368.

    const e_float one_over_nu   = ef::one() / v;
    const e_float z             = x * one_over_nu;
    const e_float sqrt_zsq_term = ef::sqrt(ef::fabs((z * z) - ef::one()));

    const e_float two_third_zeta_pow_3_2 = z_is_gt_one() ? (sqrt_zsq_term - ef::acos(ef::one() / z))
                                                         : (ef::log((ef::one() + sqrt_zsq_term) / z) - sqrt_zsq_term);

    const e_float zeta_pow_three_half          = (two_third_zeta_pow_3_2 * static_cast<INT32>(3)) / static_cast<INT32>(2);
    const e_float sqrt_zeta                    = ef::cbrt(zeta_pow_three_half);
    const e_float one_over_zeta_pow_three_half = ef::one() / zeta_pow_three_half;

    std::deque<e_float> v_t;
    std::deque<e_float> v_debye      (static_cast<std::size_t>(1u), ef::one());
    std::deque<e_float> v_zeta       (static_cast<std::size_t>(1u), ef::one());
    std::deque<e_float> v_zeta_mu    (static_cast<std::size_t>(1u), ef::one());
    std::deque<e_float> v_zeta_lambda(static_cast<std::size_t>(1u), ef::one());

    const e_float t                   = ef::one() / sqrt_zsq_term;
    const e_float one_over_nu_squared = one_over_nu * one_over_nu;

    static const INT32 max_half_k = static_cast<INT32>(Tables::A144618().size() / 2u);

    e_float sum_ak                = ef::zero();
    e_float sum_bk                = ef::zero();
    e_float one_over_nu_pow_two_k = ef::one();

    for(INT32 k = static_cast<INT32>(0); k < max_half_k; k++)
    {
      const INT32 two_k_plus_one = static_cast<INT32>((static_cast<INT32>(2) * k) + static_cast<INT32>(1));

      while(v_zeta.size() < static_cast<std::size_t>(two_k_plus_one + static_cast<INT32>(1)))
      {
        v_zeta.push_back       (v_zeta.back() * one_over_zeta_pow_three_half);
        v_zeta_mu.push_back    (v_zeta.back() * BesselAsymp::MuS()    [v_zeta_mu.size()]);
        v_zeta_lambda.push_back(v_zeta.back() * BesselAsymp::LambdaS()[v_zeta_lambda.size()]);
        v_debye.push_back      (BesselAsymp::DebyeU(static_cast<INT32>(v_debye.size()), t, v_t, z_is_gt_one()));
      }

      const e_float term_ak = std::inner_product(v_zeta_mu.begin(),
                                                 v_zeta_mu.end()  - static_cast<std::size_t>(1u),
                                                 v_debye.rbegin() + static_cast<std::size_t>(1u),
                                                 ef::zero(),
                                                 std::plus<e_float>(),
                                                 multiplies_zeta_debye_mu(z_is_gt_one())) * one_over_nu_pow_two_k;

      const e_float term_bk = (std::inner_product(v_zeta_lambda.begin(),
                                                  v_zeta_lambda.end(),
                                                  v_debye.rbegin(),
                                                  ef::zero(),
                                                  std::plus<e_float>(),
                                                  multiplies_zeta_debye_lambda(z_is_gt_one())) * one_over_nu_pow_two_k) / -sqrt_zeta;

      if((term_ak.order() < -ef::tol()) && (term_bk.order() < -ef::tol()))
      {
        break;
      }

      sum_ak += term_ak;
      sum_bk += term_bk;

      one_over_nu_pow_two_k *= one_over_nu_squared;
    }

    const e_float nu_pow_one_third = ef::cbrt(v);
    const e_float nu_pow_two_third = nu_pow_one_third * nu_pow_one_third;
    const e_float zeta             = sqrt_zeta * sqrt_zeta;

    const e_float nu_pow_two_third_zeta = nu_pow_two_third * (z_is_gt_one() ? -zeta : zeta);

    const e_float factor_ak = my_airy      (nu_pow_two_third_zeta) /  nu_pow_one_third;
          e_float factor_bk = my_airy_prime(nu_pow_two_third_zeta) / (nu_pow_two_third * v);

    const e_float result = (ef::sqrt2() * ef::sqrt(sqrt_zeta / sqrt_zsq_term)) * ((factor_ak * sum_ak) + (factor_bk * sum_bk));

    return !is_bess_y() ? result : -result;
  }

  #if(0) // Not used.
  e_float DebyeAsyBase::DebyeAsymptotic(const e_float& v, const e_float& x) const
  {
    // Abramowitz and Stegun 9.3.7 - 9.3.13, page 366.

    const e_float alpha = ef::acosh(v / x);
    const e_float tanha = ef::tanh(alpha);
    const e_float cotha = ef::one() / tanha;

    const e_float one_over_nu = ef::one() / v;

    e_float sum               = ef::one();
    e_float one_over_nu_pow_k = one_over_nu;

    std::deque<e_float> vt;

    for(INT32 k = static_cast<INT32>(1); k < static_cast<INT32>(Tables::A144618().size()); k++)
    {
      const e_float term = BesselRecursion::DebyeU(k, cotha, vt) * one_over_nu_pow_k;

      if(term.order() < -ef::tol())
      {
        break;
      }

      const bool negate_term = ((k % static_cast<INT32>(2)) != static_cast<INT32>(0)) && is_bess_y();

      !negate_term ? sum += term : sum -= term;

      one_over_nu_pow_k *= one_over_nu;
    }

    const e_float ta = tanha - alpha;

    const e_float tanha_alpha = !is_bess_y() ? ta : -ta;

    const e_float result = (ef::pow(ef::exp(tanha_alpha), v) * sum) / ef::sqrt((pi_term() * v) * tanha);

    return !is_bess_y() ? result : -result;
  }
  #endif

  INT32 UniformAsy_z_gt_one::lower(const double& xd) const
  {
    static_cast<void>(xd);

    static const UniformAsy_z_gt_one ua;

    static const double zd  = static_cast<double>(0.0);

    static const INT32 n_lo = ua.upper(zd);

    return n_lo;
  }

  INT32 UniformAsy_z_gt_one::upper(const double& xd) const
  {
    static const std::tr1::array<Util::point<double>, 24u > pt_data =
    {{
      Util::point<double>(static_cast<double>(       500.0), static_cast<double>(       500.0 -    80.0)),
      Util::point<double>(static_cast<double>(      1000.0), static_cast<double>(      1000.0 -   150.0)),
      Util::point<double>(static_cast<double>(      2000.0), static_cast<double>(      2000.0 -   300.0)),
      Util::point<double>(static_cast<double>(      4000.0), static_cast<double>(      4000.0 -   700.0)),
      Util::point<double>(static_cast<double>(      6000.0), static_cast<double>(      6000.0 -  1000.0)),
      Util::point<double>(static_cast<double>(      8000.0), static_cast<double>(      8000.0 -  1500.0)),
      Util::point<double>(static_cast<double>(     10000.0), static_cast<double>(     10000.0 -  2000.0)),
      Util::point<double>(static_cast<double>(     20000.0), static_cast<double>(     20000.0 -  2000.0)),
      Util::point<double>(static_cast<double>(     30000.0), static_cast<double>(     30000.0 -  2000.0)),
      Util::point<double>(static_cast<double>(     40000.0), static_cast<double>(     40000.0 -  2000.0)),
      Util::point<double>(static_cast<double>(     50000.0), static_cast<double>(     50000.0 -  2000.0)),
      Util::point<double>(static_cast<double>(    100000.0), static_cast<double>(    100000.0 -  2000.0)),
      Util::point<double>(static_cast<double>(    200000.0), static_cast<double>(    200000.0 -  2500.0)),
      Util::point<double>(static_cast<double>(    300000.0), static_cast<double>(    300000.0 -  3000.0)),
      Util::point<double>(static_cast<double>(   1500000.0), static_cast<double>(   1500000.0 -  3000.0)),
      Util::point<double>(static_cast<double>(   2000000.0), static_cast<double>(   2000000.0 -  3500.0)),
      Util::point<double>(static_cast<double>(   3000000.0), static_cast<double>(   3000000.0 -  4000.0)),
      Util::point<double>(static_cast<double>(   5000000.0), static_cast<double>(   5000000.0 -  4000.0)),
      Util::point<double>(static_cast<double>(  10000000.0), static_cast<double>(  10000000.0 -  6000.0)),
      Util::point<double>(static_cast<double>(  50000000.0), static_cast<double>(  50000000.0 -  8000.0)),
      Util::point<double>(static_cast<double>( 100000000.0), static_cast<double>( 100000000.0 -  8000.0)),
      Util::point<double>(static_cast<double>(1000000000.0), static_cast<double>(1000000000.0 - 22000.0)),
      Util::point<double>(static_cast<double>(1500000000.0), static_cast<double>(1500000000.0 - 24000.0)),
      Util::point<double>(static_cast<double>(2000000000.0), static_cast<double>(2000000000.0 - 26000.0))
    }};

    static const std::vector<Util::point<double> > points(pt_data.begin(), pt_data.end());

    return static_cast<INT32>(static_cast<INT64>(Util::linear_interpolate<double>::interpolate(xd, points)));
  }

  INT32 UniformAsy_z_lt_one::lower(const double& xd) const
  {
    static const std::tr1::array<Util::point<double>, 24u > pt_data =
    {{
      Util::point<double>(static_cast<double>(       500.0), static_cast<double>(       500.0 +    80.0)),
      Util::point<double>(static_cast<double>(      1000.0), static_cast<double>(      1000.0 +   150.0)),
      Util::point<double>(static_cast<double>(      2000.0), static_cast<double>(      2000.0 +   300.0)),
      Util::point<double>(static_cast<double>(      4000.0), static_cast<double>(      4000.0 +   700.0)),
      Util::point<double>(static_cast<double>(      6000.0), static_cast<double>(      6000.0 +  1000.0)),
      Util::point<double>(static_cast<double>(      8000.0), static_cast<double>(      8000.0 +  1500.0)),
      Util::point<double>(static_cast<double>(     10000.0), static_cast<double>(     10000.0 +  2000.0)),
      Util::point<double>(static_cast<double>(     20000.0), static_cast<double>(     20000.0 +  2000.0)),
      Util::point<double>(static_cast<double>(     30000.0), static_cast<double>(     30000.0 +  2000.0)),
      Util::point<double>(static_cast<double>(     40000.0), static_cast<double>(     40000.0 +  2000.0)),
      Util::point<double>(static_cast<double>(     50000.0), static_cast<double>(     50000.0 +  2000.0)),
      Util::point<double>(static_cast<double>(    100000.0), static_cast<double>(    100000.0 +  2000.0)),
      Util::point<double>(static_cast<double>(    200000.0), static_cast<double>(    200000.0 +  2500.0)),
      Util::point<double>(static_cast<double>(    300000.0), static_cast<double>(    300000.0 +  3000.0)),
      Util::point<double>(static_cast<double>(   1500000.0), static_cast<double>(   1500000.0 +  3000.0)),
      Util::point<double>(static_cast<double>(   2000000.0), static_cast<double>(   2000000.0 +  3500.0)),
      Util::point<double>(static_cast<double>(   3000000.0), static_cast<double>(   3000000.0 +  4000.0)),
      Util::point<double>(static_cast<double>(   5000000.0), static_cast<double>(   5000000.0 +  4000.0)),
      Util::point<double>(static_cast<double>(  10000000.0), static_cast<double>(  10000000.0 +  6000.0)),
      Util::point<double>(static_cast<double>(  50000000.0), static_cast<double>(  50000000.0 +  8000.0)),
      Util::point<double>(static_cast<double>( 100000000.0), static_cast<double>( 100000000.0 +  8000.0)),
      Util::point<double>(static_cast<double>(1000000000.0), static_cast<double>(1000000000.0 + 22000.0)),
      Util::point<double>(static_cast<double>(1500000000.0), static_cast<double>(1500000000.0 + 24000.0)),
      Util::point<double>(static_cast<double>(2000000000.0), static_cast<double>(2000000000.0 + 26000.0))
    }};

    static const std::vector<Util::point<double> > points(pt_data.begin(), pt_data.end());

    return static_cast<INT32>(static_cast<INT64>(Util::linear_interpolate<double>::interpolate(xd, points)));
  }

  INT32 UniformAsy_z_lt_one::upper(const double& xd) const
  {
    static_cast<void>(xd);

    return static_cast<INT32>(2000000000L);
  }
}

const std::vector<e_float>& BesselAsymp::LambdaS(void)
{
  static std::vector<e_float> lambda_s;
  
  if(lambda_s.empty())
  {
    lambda_s.resize(Tables::A144618().size());

    e_float one_over_k_fact_144_pow_k = ef::one();

    lambda_s[0] = ef::one();

    for(INT32 k = static_cast<INT32>(1u); k < static_cast<INT32>(lambda_s.size()); k++)
    {
      one_over_k_fact_144_pow_k /= static_cast<INT32>(k * static_cast<INT32>(144));

      lambda_s[k] = one_over_k_fact_144_pow_k;

      for(INT32 m  = static_cast<INT32>((static_cast<INT32>(2) * k) + static_cast<INT32>(1));
                m <= static_cast<INT32>((static_cast<INT32>(6) * k) - static_cast<INT32>(1));
                m += static_cast<INT32>(2)
         )
      {
        lambda_s[k] *= m;
      }
    }
  }

  static const std::vector<e_float> lambda_s_values(lambda_s.begin(), lambda_s.end());

  return lambda_s_values;
}

const std::vector<e_float>& BesselAsymp::MuS(void)
{
  static std::vector<e_float> mu_s;

  if(mu_s.empty())
  {
    mu_s.resize(LambdaS().size());

    mu_s[0] = ef::one();

    for(INT32 s = static_cast<INT32>(0); s < static_cast<INT32>(LambdaS().size()); s++)
    {
      const INT32 six_s           = static_cast<INT32>(static_cast<INT32>(6) * s);
      const INT32 six_s_plus_one  = static_cast<INT32>(six_s + static_cast<INT32>(1));
      const INT32 six_s_minus_one = static_cast<INT32>(six_s - static_cast<INT32>(1));

      mu_s[s] = -(LambdaS()[s] * six_s_plus_one) / six_s_minus_one;
    }
  }

  static const std::vector<e_float> mu_s_values(mu_s.begin(), mu_s.end());

  return mu_s_values;
}

e_float BesselAsymp::DebyeU(const INT32 n, const e_float& t, std::deque<e_float>& vt, const bool has_phase)
{
  if(vt.empty())
  {
    vt.push_back(ef::one());
  }

  const INT32 n_first_t = static_cast<INT32>(n / static_cast<INT32>(2));
  const INT32 n_poly_t  = static_cast<INT32>(n + static_cast<INT32>(1));

  const e_float t_sq = t * t;

  while(static_cast<INT32>(n_first_t + n_poly_t) >= static_cast<INT32>(vt.size()))
  {
    const bool neg_phase = has_phase && ((static_cast<INT32>(vt.size()) % static_cast<std::size_t>(2u)) != static_cast<INT32>(0));

    const e_float t_sq_2n = ef::fabs(vt.back()) * t_sq;

    vt.push_back(!neg_phase ? t_sq_2n : -t_sq_2n);
  }

  const e_float sum_top = std::inner_product(Tables::A144617()[n]().begin(),
                                             Tables::A144617()[n]().end(),
                                             vt.begin() + static_cast<std::size_t>(n_first_t),
                                             ef::zero());

  const e_float sum = sum_top / Tables::A144618()[n]();

  const bool b_n_is_even = ((n % static_cast<INT32>(2)) == static_cast<INT32>(0));

  return (b_n_is_even ? sum : (sum * t));
}
