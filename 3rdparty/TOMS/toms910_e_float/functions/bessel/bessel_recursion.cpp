
#include <deque>
#include <numeric>
#include <algorithm>

#include <e_float/e_float.h>
#include <functions/bessel/bessel_recursion.h>
#include <functions/bessel/bessel_recursion_order.h>
#include <functions/constants/constants.h>
#include <functions/elementary/elementary.h>
#include <functions/gamma/gamma.h>
#include <functions/tables/tables.h>
#include <utility/util_interpolate.h>
#include <utility/util_point.h>

e_float BesselRecursion::RecurJn(const INT32 n, const e_float& x, std::deque<e_float>* const pJn)
{
  const bool b_neg = ef::isneg(x);
  
  const e_float xx = !b_neg ? x : -x;

  const double xd = ef::to_double(xx);

  // Find the recursion order which is necessary for downward recursion and
  // ensure that the value of the order is even numbered.
  const INT32 N0 = BesselRecursionOrder::RecursionStartOrderJ0(   xd);
  const INT32 Nn = BesselRecursionOrder::RecursionStartOrderJn(n, xd);
  const INT32 N2 = static_cast<INT32>(n + static_cast<INT32>(2));
  const INT32 Nm = static_cast<INT32>(static_cast<INT32>(2) * (std::max(N2, std::max(N0, Nn)) / static_cast<INT32>(2)));

  // Use recursion and normalization with a Neumann sum.

  e_float Jm_p2 = ef::zero();
  e_float Jm_p1 = ef::one();
  e_float Jm;
  e_float Jmn;

  e_float sum_m = ef::two();

  const e_float two_over_x = ef::two() / xx;

  const bool has_jn_recur = pJn != static_cast<std::deque<e_float>* const>(0u);

  if(has_jn_recur) { pJn->resize(Nm); }

  // Do the downward recursion of Jn:
  //
  //                  Jn+1
  //   Jn = [ 2 (n+1) ---- ] - Jn+2
  //                   x 

  for(INT32 m = static_cast<INT32>(Nm - static_cast<INT32>(1)); m >= static_cast<INT32>(0); m--)
  {
    Jm    = ((Jm_p1 * two_over_x) * static_cast<INT32>(m + static_cast<INT32>(1))) - Jm_p2;
    Jm_p2 = Jm_p1;
    Jm_p1 = Jm;

    if(has_jn_recur) { pJn->at(m) = Jm; }

    // Store the value of the Bessel function which has the sought order.
    if(n == m)
    {
      Jmn = Jm;
    }

    // Do the normalization using a Neumann sum which is:
    // 1 = J_0(x) + 2 * J_2(x) + 2 * J_4(x) + 2 * J_6(x) + ...
    if((m % static_cast<INT32>(2)) == static_cast<INT32>(0))
    {
      sum_m += (m == static_cast<INT32>(0) ? Jm : Jm * static_cast<INT32>(2));
    }
  }

  const e_float norm = ef::one() / sum_m;

  if(has_jn_recur)
  {
    std::transform(pJn->begin(),
                   pJn->end(),
                   pJn->begin(),
                   std::bind2nd(std::multiplies<e_float>(), norm));

    if(b_neg)
    {
      bool b_negate_term = false;

      for(std::size_t i = static_cast<std::size_t>(0u); i < pJn->size(); i++)
      {
        if(b_negate_term) { pJn->at(i) = -pJn->at(i); }

        b_negate_term = !b_negate_term;
      }
    }
  }

  if(b_neg)
  {
    if(static_cast<INT32>(n % static_cast<INT32>(2)) != static_cast<INT32>(0))
    {
      Jmn = -Jmn;
    }
  }

  return Jmn * norm;
}

namespace BesselRecursion
{
  template<typename T> static inline T RecurJvTemplateNormOne(const e_float& v, const T& x, std::deque<T>* const pJv)
  {
    using ef::abs;
    using efz::abs;
    using ef::pow;
    using efz::pow;

    const double xd = ef::to_double(abs(x));
    const INT32  n  = ef::to_int32(v);

    // Find the recursion order which is necessary for downward recursion and
    // ensure that the value of the order is even numbered.
    const INT32 N0 = BesselRecursionOrder::RecursionStartOrderJ0(xd);
    const INT32 Nn = BesselRecursionOrder::RecursionStartOrderJn(n, xd);
    const INT32 N2 = static_cast<INT32>(n + static_cast<INT32>(2));
    const INT32 Nm = static_cast<INT32>(static_cast<INT32>(2) * (std::max(N2, std::max(N0, Nn)) / static_cast<INT32>(2)));

    // Use recursion and normalization with a Neumann sum.

    const e_float v_frac = ef::decimal_part(v);

    INT32   k                   = static_cast<INT32>(Nm / static_cast<INT32>(2));
    e_float one_over_k_fact     = ef::one() / ef::factorial(static_cast<UINT32>(k));
    e_float k_plus_v_frac       = v_frac + k;
    e_float gamma_k_plus_v_frac = ef::gamma(k_plus_v_frac);

    T Jv_p2 = ef::zero();
    T Jv_p1 = ef::one();
    T Jv;
    T Jvn;

    const e_float v_frac_plus_two_k_over_k_fact = (v_frac + static_cast<INT32>(k * static_cast<INT32>(2))) * one_over_k_fact;

    T sum_v_frac = v_frac_plus_two_k_over_k_fact * gamma_k_plus_v_frac;

    const bool has_jv_recur = pJv != static_cast<std::deque<T>* const>(0u);

    if(has_jv_recur) { pJv->resize(Nm); }

    // Do the downward recursion of Jv:
    //
    //                  Jv+1
    //   Jv = [ 2 (v+1) ---- ] - Jv+2
    //                   x 

    const T two_over_x = ef::two() / x;

    for(INT32 m = static_cast<INT32>(Nm - static_cast<INT32>(1)); m >= static_cast<INT32>(0); m--)
    {
      const e_float n_plus_v_frac_plus_one = v_frac + static_cast<INT32>(m + static_cast<INT32>(1));

      // Downward recursion for Jv_frac.
      Jv    = ((n_plus_v_frac_plus_one * Jv_p1) * two_over_x) - Jv_p2;
      Jv_p2 = Jv_p1;
      Jv_p1 = Jv;

      if(has_jv_recur) { pJv->at(m) = Jv; }

      // Store the value of the Bessel function which has the sought order.
      if(n == m)
      {
        Jvn = Jv;
      }

      // Do the normalization using a Neumann sum which is:
      // (x/2)^v = Sum_k { [((v + 2k) gamma(v + k)) / k!] * J_v+2k }
      if((m % static_cast<INT32>(2)) == static_cast<INT32>(0))
      {
        one_over_k_fact     *= k--;
        gamma_k_plus_v_frac /= --k_plus_v_frac;

        // Increment the normalization sum for Jv_frac.
        sum_v_frac += (((v_frac + m) * gamma_k_plus_v_frac) * Jv) * one_over_k_fact;
      }
    }

    const T x_half_pow_v_frac = pow(x / static_cast<INT32>(2), T(v_frac));

    // Compute the normalizations using the Neumann sums in combination with (x/2)^v.
    const T norm =  x_half_pow_v_frac / sum_v_frac;

    if(has_jv_recur)
    {
      std::transform(pJv->begin(),
                     pJv->end(),
                     pJv->begin(),
                     std::bind2nd(std::multiplies<T>(), norm));
    }

    // Normalize In_v_frac.
    return Jvn * norm;
  }

  template<typename T> static inline T RecurJvTemplateNormCos(const e_float& v, const T& x, std::deque<T>* const pJv)
  {
    using ef::abs;
    using efz::abs;
    using ef::cos;
    using efz::cos;
    using ef::pow;
    using efz::pow;

    const double xd = ef::to_double(abs(x));
    const INT32  n  = ef::to_int32(v);

    // Find the recursion order which is necessary for downward recursion and
    // ensure that the value of the order is even numbered.
    const INT32 N0 = BesselRecursionOrder::RecursionStartOrderJ0(   xd);
    const INT32 Nn = BesselRecursionOrder::RecursionStartOrderJn(n, xd);
    const INT32 N2 = static_cast<INT32>(n + static_cast<INT32>(2));
    const INT32 Nm = static_cast<INT32>(static_cast<INT32>(2) * (std::max(N2, std::max(N0, Nn)) / static_cast<INT32>(2)));

    // Use recursion and normalization with a Neumann sum.

    const e_float v_frac = ef::decimal_part(v);

    e_float m_fact                  = ef::factorial(static_cast<UINT32>(Nm));
    e_float v_frac_plus_m           = v_frac + Nm;
    e_float two_v_frac_plus_m       = v_frac + v_frac_plus_m;
    e_float gamma_two_v_frac_plus_m = ef::gamma(two_v_frac_plus_m);

    T Jv_p2 = ef::zero();
    T Jv_p1 = ef::one();
    T Jv;
    T Jvn;

    bool b_neg_sum_term = static_cast<INT32>(static_cast<INT32>(Nm / static_cast<INT32>(2)) % static_cast<INT32>(2)) != static_cast<INT32>(0);

    T term = (v_frac_plus_m * gamma_two_v_frac_plus_m) / m_fact;

    T sum_v_frac = !b_neg_sum_term ? term : -term;

    const bool has_jv_recur = pJv != static_cast<std::deque<T>* const>(0u);

    if(has_jv_recur) { pJv->resize(Nm); }

    // Do the downward recursion of Jv:
    //
    //                  Jv+1
    //   Jv = [ 2 (v+1) ---- ] - Jv+2
    //                   x 

    const T two_over_x = ef::two() / x;

    for(INT32 m = static_cast<INT32>(Nm - static_cast<INT32>(1)); m >= static_cast<INT32>(0); m--)
    {
      // Downward recursion for Jv_frac.
      Jv    = ((v_frac_plus_m * Jv_p1) * two_over_x) - Jv_p2;
      Jv_p2 = Jv_p1;
      Jv_p1 = Jv;

      if(has_jv_recur) { pJv->at(m) = Jv; }

      // Store the value of the Bessel function which has the sought order.
      if(n == m)
      {
        Jvn = Jv;
      }

      m_fact /= static_cast<INT32>(m + static_cast<INT32>(1));
      gamma_two_v_frac_plus_m /= --two_v_frac_plus_m;
      --v_frac_plus_m;

      // Do the normalization using a Neumann sum which is:
      // cos(z) { [(z/2)^v] [gamma(2v)/gamma(v)] } = Sum_m { (-1)^m [(v + 2m) J_v+2m gamma(2v + 2m)] / (2m)! }
      if((m % static_cast<INT32>(2)) == static_cast<INT32>(0))
      {
        b_neg_sum_term = static_cast<INT32>(static_cast<INT32>(m / static_cast<INT32>(2)) % static_cast<INT32>(2)) != static_cast<INT32>(0);

        term = ((v_frac_plus_m * Jv) * gamma_two_v_frac_plus_m) / m_fact;

        // Increment the normalization sum for Jv_frac.
        !b_neg_sum_term ? sum_v_frac += term : sum_v_frac -= term;
      }
    }

    // Compute the normalization.

    const T       x_half_pow_v_frac = pow(x / static_cast<INT32>(2), T(v_frac));
    const e_float gamma_term        = ef::gamma(v_frac * static_cast<INT32>(2)) / ef::gamma(v_frac);

    const T norm = ((cos(x) * x_half_pow_v_frac) * gamma_term) / sum_v_frac;

    if(has_jv_recur)
    {
      std::transform(pJv->begin(),
                     pJv->end(),
                     pJv->begin(),
                     std::bind2nd(std::multiplies<T>(), norm));
    }

    // Normalize In_v_frac.
    return Jvn * norm;
  }
}

e_float BesselRecursion::RecurJv(const e_float& v, const e_float& x, std::deque<e_float>* const pJv)
{
  if(ef::isint(v))
  {
    return RecurJn(static_cast<INT32>(ef::to_int64(v)), x, pJv);
  }

  return RecurJvTemplateNormOne<e_float>(v, x, pJv);
}

ef_complex BesselRecursion::RecurJv(const e_float& v, const ef_complex& z, std::deque<ef_complex>* const pJv)
{
/*
  if(ef::fabs(z.imag()) < (ef::fabs(z.real()) / static_cast<INT32>(4)))
  {
    return RecurJvTemplateNormOne<ef_complex>(v, z, pJv);
  }
*/
  return RecurJvTemplateNormCos<ef_complex>(v, z, pJv);
}

e_float BesselRecursion::RecurIn(const INT32 n, const e_float& x, std::deque<e_float>* const pIn)
{
  const double xd = ef::to_double(x);

  // Find the recursion order which is necessary for downward recursion and
  // ensure that the value of the order is even numbered.
  const INT32 N0 = BesselRecursionOrder::RecursionStartOrderI0(   xd);
  const INT32 Nn = BesselRecursionOrder::RecursionStartOrderIn(n, xd);
  const INT32 N2 = static_cast<INT32>(n + static_cast<INT32>(2));
  const INT32 Nm = static_cast<INT32>(static_cast<INT32>(2) * (std::max(N2, std::max(N0, Nn)) / static_cast<INT32>(2)));

  // Use recursion and normalization with a Neumann sum.

  e_float Im_p2 = ef::zero();
  e_float Im_p1 = ef::one();
  e_float Im;
  e_float Imn;

  e_float sum_m = ef::two();

  const e_float two_over_x = ef::two() / x;

  const bool has_in_recur = pIn != static_cast<std::deque<e_float>* const>(0u);

  if(has_in_recur) { pIn->resize(Nm); }

  // Do the downward recursion of Jn:
  //
  //                  In+1
  //   In = [ 2 (n+1) ---- ] + In+2
  //                   x 

  for(INT32 m = static_cast<INT32>(Nm - static_cast<INT32>(1)); m >= static_cast<INT32>(0); m--)
  {
    Im    = ((Im_p1 * two_over_x) * static_cast<INT32>(m + static_cast<INT32>(1))) + Im_p2;
    Im_p2 = Im_p1;
    Im_p1 = Im;

    if(has_in_recur) { pIn->at(m) = Im; }

    // Store the value of the Bessel function which has the sought order.
    if(n == m)
    {
      Imn = Im;
    }

    // Do the normalization using a Neumann sum which is:
    // exp(x) = I_0(x) + 2 * I_1(x) + 2 * I_2(x) + 2 * I_3(x) + ...
    sum_m += (m == static_cast<INT32>(0) ? Im : Im * static_cast<INT32>(2));
  }

  const e_float norm = ef::exp(x) / sum_m;

  if(has_in_recur)
  {
    std::transform(pIn->begin(),
                   pIn->end(),
                   pIn->begin(),
                   std::bind2nd(std::multiplies<e_float>(), norm));
  }

  return Imn * norm;
}

e_float BesselRecursion::RecurIv(const e_float& v, const e_float& x, std::deque<e_float>* const pIv)
{
  if(ef::isint(v))
  {
    return RecurIn(static_cast<INT32>(ef::to_int64(v)), x, pIv);
  }

  const double xd = ef::to_double(x);
  const INT32  n  = ef::to_int32(v);

  // Find the recursion order which is necessary for downward recursion and
  // ensure that the value of the order is even numbered.
  const INT32 N0 = BesselRecursionOrder::RecursionStartOrderI0(   xd);
  const INT32 Nn = BesselRecursionOrder::RecursionStartOrderIn(n, xd);
  const INT32 N2 = static_cast<INT32>(n + static_cast<INT32>(2));
  const INT32 Nm = static_cast<INT32>(static_cast<INT32>(2) * (std::max(N2, std::max(N0, Nn)) / static_cast<INT32>(2)));

  // Use recursion and normalization with a Neumann sum.

  const e_float v_frac = ef::decimal_part(v);

  e_float m_plus_v_frac           = v_frac + Nm;
  e_float m_plus_two_v_frac       = m_plus_v_frac + v_frac;
  e_float gamma_m_plus_two_v_frac = ef::gamma(m_plus_two_v_frac);
  e_float one_over_m_fact         = ef::one() / ef::factorial(static_cast<UINT32>(Nm));

  e_float Iv_p2 = ef::zero();
  e_float Iv_p1 = ef::one();
  e_float Iv;
  e_float Ivn;

  e_float n_plus_one_plus_v_frac = m_plus_v_frac;

  e_float sum_v_frac = (gamma_m_plus_two_v_frac * m_plus_v_frac) * one_over_m_fact;


  const bool has_iv_recur = pIv != static_cast<std::deque<e_float>* const>(0u);

  if(has_iv_recur) { pIv->resize(Nm); }

  // Do the downward recursion of In:
  //
  //                  Iv+1
  //   In = [ 2 (v+1) ---- ] + Iv+2
  //                   x 

  const e_float two_over_x = ef::two() / x;

  for(INT32 m = static_cast<INT32>(Nm - static_cast<INT32>(1)); m >= static_cast<INT32>(0); m--)
  {
    Iv    = ((Iv_p1 * two_over_x) * n_plus_one_plus_v_frac) + Iv_p2;
    Iv_p2 = Iv_p1;
    Iv_p1 = Iv;

    if(has_iv_recur) { pIv->at(m) = Iv; }

    // Store the value of the Bessel function which has the sought order.
    if(n == m)
    {
      Ivn = Iv;
    }

    // Do the normalization using a Neumann sum which is:
    // [(e^x) ((x/2)^v)] * [gamma(2v)/gamma(v)] = Sum_m { [((v + m) gamma(2v + m)) / m!] * I_v+m }

    gamma_m_plus_two_v_frac /= --m_plus_two_v_frac;
    --m_plus_v_frac;
    one_over_m_fact *= static_cast<INT32>(m + static_cast<INT32>(1));

    --n_plus_one_plus_v_frac;

    sum_v_frac += ((gamma_m_plus_two_v_frac * m_plus_v_frac) * one_over_m_fact) * Iv;
  }

  const e_float x_half_pow_v_frac = ef::pow(x / static_cast<INT32>(2), v_frac);

  // Compute the normalization.
  const e_float factor = ((ef::exp(x) * x_half_pow_v_frac) * gamma_m_plus_two_v_frac) / ef::gamma(v_frac);
  const e_float norm   = factor / sum_v_frac;

  if(has_iv_recur)
  {
    std::transform(pIv->begin(),
                   pIv->end(),
                   pIv->begin(),
                   std::bind2nd(std::multiplies<e_float>(), norm));
  }

  // Normalize In_v_frac.
  return Ivn * norm;
}
