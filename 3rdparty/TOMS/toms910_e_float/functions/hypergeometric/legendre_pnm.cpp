
#include <functions/constants/constants.h>
#include <functions/elementary/elementary.h>
#include <functions/gamma/gamma.h>
#include <functions/hypergeometric/legendre.h>
#include <functions/hypergeometric/legendre_xnm.h>
#include <functions/polynomials/polynomials.h>

namespace Legendre_Series
{
  class LegendrePnm : public LegendreXnm
  {
  public:

    LegendrePnm() { }

    virtual ~LegendrePnm() { }

  private:

    virtual e_float AtIdenticallyZero        (const INT32 n, const INT32 m) const;
    virtual e_float AtIdenticallyOne         (const INT32 n, const INT32 m) const { static_cast<void>(n); static_cast<void>(m); return ef::zero(); } // NOCOVER_LINE
    virtual e_float AtReflectNegativeDegree  (const INT32 n, const INT32 m, const e_float& x) const;
    virtual e_float AtReflectNegativeArgument(const INT32 n, const INT32 m, const e_float& x) const;

    virtual bool AtResultForDegreeAndOrderIsZero(const INT32 n, const INT32 m) const { return static_cast<INT32>(m < static_cast<INT32>(0) ? -m : m) > n; }

    virtual e_float L_00(const e_float& x) const;
    virtual e_float L_10(const e_float& x) const;
    virtual e_float L_01(const e_float& x) const;
    virtual e_float L_11(const e_float& x) const;

    virtual e_float MyRecursion(const INT32 n, const INT32 m, const e_float& x) const;

    virtual e_float Legendre_n(const INT32 n, const e_float& x) const { return ef::legendre_p(n, x); } // NOCOVER_LINE
  };
}

e_float Legendre_Series::LegendrePnm::AtIdenticallyZero(const INT32 n, const INT32 m) const
{
  const INT32 n_plus_m         = static_cast<INT32>(n + m);
  const bool  n_plus_m_is_even = static_cast<INT32>(n_plus_m % static_cast<INT32>(2)) == static_cast<INT32>(0);

  if(n_plus_m_is_even)
  {
    const INT32 n_plus_m_over_two = static_cast<INT32>(n_plus_m / static_cast<INT32>(2));
    const bool  b_negate          = static_cast<INT32>(n_plus_m_over_two % static_cast<INT32>(2)) != static_cast<INT32>(0);

    const UINT32 n_plus_m_minus_one = static_cast<UINT32>(static_cast<UINT32>(n_plus_m) - static_cast<UINT32>(1u));
    const UINT32 n_minus_m          = static_cast<UINT32>(n - m);

    const e_float Pnm = ef::factorial2(n_plus_m_minus_one) / ef::factorial2(n_minus_m);

    return !b_negate ? Pnm : -Pnm; 
  }
  else
  {
    return ef::zero();
  }
}

e_float Legendre_Series::LegendrePnm::AtReflectNegativeDegree(const INT32 n, const INT32 m, const e_float& x) const
{
  // Handle negative degree with n < -1 using
  // Abramowitz & Stegun 8.2.1, P_[-n-1](z) = P_[n](z).
  return MyLegendre(static_cast<INT32>(static_cast<INT32>(-n) - static_cast<INT32>(1)), m, x);
}

e_float Legendre_Series::LegendrePnm::AtReflectNegativeArgument(const INT32 n, const INT32 m, const e_float& x) const
{
  const INT32 n_minus_m = static_cast<INT32>(n - m);

  const bool b_negate = static_cast<INT32>(n_minus_m % static_cast<INT32>(2)) != static_cast<INT32>(0);

  const e_float Pnm = MyLegendre(n, m, -x);

  return !b_negate ? Pnm : -Pnm;
}

e_float Legendre_Series::LegendrePnm::L_00(const e_float& x) const
{
  return ef::legendre_p(static_cast<INT32>(0), x);
}

e_float Legendre_Series::LegendrePnm::L_10(const e_float& x) const
{
  return ef::legendre_p(static_cast<INT32>(1), x);
}

e_float Legendre_Series::LegendrePnm::L_01(const e_float& x) const
{
  static_cast<void>(x);
  return ef::zero();
}

e_float Legendre_Series::LegendrePnm::L_11(const e_float& x) const
{
  return -ef::sqrt(ef::fabs(ef::one() - (x * x)));
}

e_float Legendre_Series::LegendrePnm::MyRecursion(const INT32 n, const INT32 m, const e_float& x) const
{
  // Calculate Pn,n and Pn+1,n and use upward recursion of the degree.

  INT32 nr    = m;
  INT32 two_n = static_cast<INT32>(2) * nr;

  const bool m_is_odd = ((m % static_cast<INT32>(2)) != static_cast<INT32>(0));

  const e_float one_minus_x_squared            = ef::one() - (x * x);
  const e_float one_minus_x_squared_pow_m_half = ef::pown(one_minus_x_squared, static_cast<INT64>(m / static_cast<INT32>(2)));
  const e_float two_n_minus_one_fact2          = ef::factorial2(two_n - static_cast<INT32>(1));

  e_float Pnm_m2 = two_n_minus_one_fact2 * one_minus_x_squared_pow_m_half;
  
  if(m_is_odd)
  {
    Pnm_m2 = -Pnm_m2 * ef::sqrt(one_minus_x_squared);
  }

  if(n == m)
  {
    return Pnm_m2;
  }

  e_float Pnm_m1 = x * (two_n + static_cast<INT32>(1)) * Pnm_m2;

  if(n == static_cast<INT32>(m + static_cast<INT32>(1)))
  {
    return Pnm_m1;
  }

  e_float Pnm;

  for(nr = static_cast<INT32>(m + 2); nr <= n; nr++)
  {
    // Use upward recursion of the degree.
    const INT32 two_n_minus_one    = (nr * static_cast<INT32>(2)) - static_cast<INT32>(1);
    const INT32 n_plus_m_minus_one = static_cast<INT32>(nr + m) - static_cast<INT32>(1);
    const INT32 n_minus_m          = static_cast<INT32>(nr - m);

    Pnm = (((two_n_minus_one * x) * Pnm_m1) - (n_plus_m_minus_one * Pnm_m2)) / n_minus_m;

    Pnm_m2 = Pnm_m1;
    Pnm_m1 = Pnm;
  }

  return Pnm;
}

e_float ef::legendre_p(const INT32 n, const INT32 m, const e_float& x)
{
  static const Legendre_Series::LegendrePnm L_Pnm;
  return L_Pnm.MyLegendre(n, m, x);
}
