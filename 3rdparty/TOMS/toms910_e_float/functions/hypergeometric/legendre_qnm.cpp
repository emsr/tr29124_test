
#include <e_float/e_float.h>
#include <functions/constants/constants.h>
#include <functions/elementary/elementary.h>
#include <functions/gamma/gamma.h>
#include <functions/hypergeometric/legendre.h>
#include <functions/hypergeometric/legendre_xnm.h>
#include <functions/polynomials/polynomials.h>

namespace Legendre_Series
{
  class LegendreQnm : public LegendreXnm
  {
  public:

    LegendreQnm() { }

    virtual ~LegendreQnm() { }

  private:

    virtual e_float AtIdenticallyZero        (const INT32 n, const INT32 m) const;
    virtual e_float AtIdenticallyOne         (const INT32 n, const INT32 m) const { static_cast<void>(n); static_cast<void>(m); return std::numeric_limits<e_float>::quiet_NaN(); } // NOCOVER_LINE
    virtual e_float AtReflectNegativeDegree  (const INT32 n, const INT32 m, const e_float& x) const;
    virtual e_float AtReflectNegativeArgument(const INT32 n, const INT32 m, const e_float& x) const;

    virtual bool AtResultForDegreeAndOrderIsZero(const INT32 n, const INT32 m) const { static_cast<void>(n); static_cast<void>(m); return false; }

    virtual e_float L_00(const e_float& x) const;
    virtual e_float L_10(const e_float& x) const;
    virtual e_float L_01(const e_float& x) const;
    virtual e_float L_11(const e_float& x) const;

    virtual e_float MyRecursion(const INT32 n, const INT32 m, const e_float& x) const;

    virtual e_float Legendre_n(const INT32 n, const e_float& x) const { return ef::legendre_q(n, x); } // NOCOVER_LINE
  };
}

e_float Legendre_Series::LegendreQnm::AtIdenticallyZero(const INT32 n, const INT32 m) const
{
  const INT32 n_plus_m        = static_cast<INT32>(n + m);
  const INT32 n_plus_m_is_odd = static_cast<INT32>(n_plus_m % static_cast<INT32>(2)) != static_cast<INT32>(0);

  if(n_plus_m_is_odd)
  {
    const INT32 n_plus_m_plus_one          = static_cast<INT32>(n_plus_m + static_cast<INT32>(1));
    const INT32 n_plus_m_plus_one_over_two = static_cast<INT32>(n_plus_m_plus_one / static_cast<INT32>(2));
    const bool  b_negate                   = static_cast<INT32>(n_plus_m_plus_one_over_two % static_cast<INT32>(2)) != static_cast<INT32>(0);

    const UINT32 n_plus_m_minus_one = static_cast<UINT32>(static_cast<UINT32>(n_plus_m) - static_cast<UINT32>(1u));
    const UINT32 n_minus_m          = static_cast<UINT32>(n - m);

    const e_float Qnm = ef::factorial2(n_plus_m_minus_one) / ef::factorial2(n_minus_m);

    return !b_negate ? Qnm : -Qnm; 
  }
  else
  {
    return ef::zero();
  }
}

e_float Legendre_Series::LegendreQnm::AtReflectNegativeDegree(const INT32 n, const INT32 m, const e_float& x) const
{
  static_cast<void>(n);
  static_cast<void>(m);
  static_cast<void>(x);
  return std::numeric_limits<e_float>::quiet_NaN();
}

e_float Legendre_Series::LegendreQnm::AtReflectNegativeArgument(const INT32 n, const INT32 m, const e_float& x) const
{
  const INT32 n_plus_m_plus_one = static_cast<INT32>(static_cast<INT32>(n + m) + static_cast<INT32>(1));

  const bool b_negate = static_cast<INT32>(n_plus_m_plus_one % static_cast<INT32>(2)) != static_cast<INT32>(0);

  const e_float Qnm = MyLegendre(n, m, -x);

  return (!b_negate ? Qnm : -Qnm);
}

e_float Legendre_Series::LegendreQnm::L_00(const e_float& x) const
{
  return ef::legendre_q(static_cast<INT32>(0), x);
}

e_float Legendre_Series::LegendreQnm::L_10(const e_float& x) const
{
  return ef::legendre_q(static_cast<INT32>(1), x);
}

e_float Legendre_Series::LegendreQnm::L_01(const e_float& x) const
{
  return -ef::one() / ef::sqrt(ef::one() - (x * x));
}

e_float Legendre_Series::LegendreQnm::L_11(const e_float& x) const
{
  const e_float log_term = ef::small_arg(x) ?  ef::log1p1m2(x)
                                            : (ef::log((ef::one() + x) / (ef::one() - x)) / static_cast<INT32>(2));

  const e_float one_minus_x_squared = ef::one() - (x * x);

  return -ef::sqrt(one_minus_x_squared) * (log_term + (x / one_minus_x_squared));
}

e_float Legendre_Series::LegendreQnm::MyRecursion(const INT32 n, const INT32 m, const e_float& x) const
{
  // First recur the degree n upward for two columns of the order m.
  const e_float Qn0 = Legendre_n(n, x);

  e_float Qn1_m2 = L_01(x);
  e_float Qn1_m1 = L_11(x);
  e_float Qn1;

  const INT32 m1 = static_cast<INT32>(1);

  for(INT32 k = static_cast<INT32>(2); k <= n; k++)
  {
    // Use upward recursion of the degree.
    const INT32 two_n_minus_one    = static_cast<INT32>((k * static_cast<INT32>(2)) - static_cast<INT32>(1));
    const INT32 m_plus_n_minus_one = static_cast<INT32>(static_cast<INT32>(k + m1) - static_cast<INT32>(1));
    const INT32 n_minus_m          = static_cast<INT32>(k - m1);

    Qn1 = (((two_n_minus_one * x) * Qn1_m1) - (m_plus_n_minus_one * Qn1_m2)) / n_minus_m;

    Qn1_m2 = Qn1_m1;
    Qn1_m1 = Qn1;
  }

  // Use rightward (upward) recursion of the order.
  const e_float x_squared = x * x;

  const e_float sqrt_x_term = ef::sqrt(ef::one() - x_squared);

  e_float Qnm_mm2 = Qn0;
  e_float Qnm_mm1 = Qn1;
  e_float Qnm_m;

  for(INT32 j = static_cast<INT32>(2); j <= m; j++)
  {
    const INT32 n_plus_m_minus_one = static_cast<INT32>(static_cast<INT32>(n + j) - static_cast<INT32>(1));
    const INT32 n_minus_m_plus_two = static_cast<INT32>(static_cast<INT32>(n - j) + static_cast<INT32>(2));
    const INT32 two_m_minus_two    = static_cast<INT32>(static_cast<INT32>(j - static_cast<INT32>(1)) * static_cast<INT32>(2));

    const e_float term1 = -((two_m_minus_two * x) * Qnm_mm1) / sqrt_x_term;

    const e_float term2 = -((n_plus_m_minus_one * Qnm_mm2) * n_minus_m_plus_two);

    Qnm_m = term1 + term2;

    Qnm_mm2 = Qnm_mm1;
    Qnm_mm1 = Qnm_m;
  }

  return Qnm_m;
}

e_float ef::legendre_q(const INT32 n, const INT32 m, const e_float& x)
{
  static const Legendre_Series::LegendreQnm L_Qnm;
  return L_Qnm.MyLegendre(n, m, x);
}
