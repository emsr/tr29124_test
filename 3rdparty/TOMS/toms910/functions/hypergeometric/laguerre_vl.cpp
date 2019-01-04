
#include <functions/constants/constants.h>
#include <functions/elementary/elementary.h>
#include <functions/gamma/gamma.h>
#include <functions/hypergeometric/hypergeometric.h>
#include <functions/hypergeometric/laguerre.h>
#include <functions/polynomials/polynomials.h>

e_float ef::laguerre(const e_float& v, const e_float& x)
{
  if(ef::isint(v))
  {
    return ef::laguerre(ef::to_int32(v), x);
  }
  else
  {
    return ef::hyperg_1f1(-v, ef::one(), x);
  }
}

e_float ef::laguerre(const e_float& v, const e_float& L, const e_float& x)
{
  if(ef::isint(v))
  {
    return ef::laguerre(ef::to_int32(v), L, x);
  }
  else
  {
    const e_float nu_plus_one = v + ef::one();

    const e_float H1F1 = ef::hyperg_1f1_reg(-v, L + ef::one(), x);

    return (ef::gamma(L + nu_plus_one) * H1F1) / ef::gamma(nu_plus_one);
  }
}

e_float ef::laguerre(const INT32 n, const e_float& L, const e_float& x)
{
  if(ef::isint(L))
  {
    return ef::laguerre(n, static_cast<INT32>(ef::to_int64(L)), x);
  }

  if(n < static_cast<INT32>(0))
  {
    const e_float nu_plus_one = n + ef::one();

    const e_float H1F1 = ef::hyperg_1f1_reg(e_float(static_cast<UINT32>(-n)), L + ef::one(), x);

    return (ef::gamma(L + nu_plus_one) * H1F1) / ef::gamma(nu_plus_one);
  }

  const e_float L_minus_one = L - ef::one();

  e_float Ln_m2 = ef::one();
  e_float Ln_m1 = (L + ef::one()) - x;
  e_float Ln;

  if(n == static_cast<INT32>(0))
  {
    return Ln_m2;
  }

  if(n == static_cast<INT32>(1))
  {
    return Ln_m1;
  }

  for(INT32 i = static_cast<INT32>(2); i <= n; i++)
  {
    const INT32 two_i = static_cast<INT32>(static_cast<INT32>(2) * i);

    Ln = ((((two_i + L_minus_one) - x) * Ln_m1) - ((i + L_minus_one) * Ln_m2)) / i;

    Ln_m2 = Ln_m1;
    Ln_m1 = Ln;
  }

  return Ln;
}

e_float ef::laguerre(const INT32 n, const INT32 m, const e_float& x)
{
  if(n < static_cast<INT32>(0))
  {
    const e_float nu_plus_one = n + ef::one();

    const e_float H1F1 = ef::hyperg_1f1_reg(e_float(static_cast<UINT32>(-n)), e_float(m) + ef::one(), x);

    return (ef::gamma(m + nu_plus_one) * H1F1) / ef::gamma(nu_plus_one);
  }

  const e_float m_minus_one = e_float(m - static_cast<INT32>(1));

  e_float Ln_m2 = ef::one();
  e_float Ln_m1 = (m + ef::one()) - x;
  e_float Ln;

  if(n == static_cast<INT32>(0))
  {
    return Ln_m2;
  }
  else if(n == static_cast<INT32>(1))
  {
    return Ln_m1;
  }
  else
  {
    for(INT32 i = static_cast<INT32>(2); i <= n; i++)
    {
      const INT32 two_i = static_cast<INT32>(static_cast<INT32>(2) * i);

      Ln = ((((two_i + m_minus_one) - x) * Ln_m1) - ((i + m_minus_one) * Ln_m2)) / i;

      Ln_m2 = Ln_m1;
      Ln_m1 = Ln;
    }

    return Ln;
  }
}
