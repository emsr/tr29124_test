
#include <functions/complex/e_float_complex.h>
#include <functions/constants/constants.h>
#include <functions/elementary/elementary.h>
#include <functions/polynomials/polynomials.h>
#include <functions/hypergeometric/hermite.h>
#include <functions/hypergeometric/hypergeometric.h>

namespace OrthogonalPoly_Series
{
  typedef enum enumPolyType
  {
    ChebyshevT_Type = 1,
    ChebyshevU_Type = 2,
    LaguerreL_Type  = 3,
    HermiteH_Type   = 4
  }
  tPolyType;

  template<typename T> static inline T OrthogonalPoly_Template(const T& x, const UINT32 n, const tPolyType type, std::vector<T>* const vp = static_cast<std::vector<T>*>(0u))
  {
    // Compute the value of an orthogonal polynomial of one of the following types:
    // Chebyshev 1st, Chebyshev 2nd, Laguerre, or Hermite

    if(vp != static_cast<std::vector<T>*>(0u))
    {
      vp->clear();
      vp->reserve(static_cast<std::size_t>(n + 1u));
    }

    T y0 = ef::one();

    if(vp != static_cast<std::vector<T>*>(0u))
    {
      vp->push_back(y0);
    }

    if(n == static_cast<UINT32>(0u))
    {
      return y0;
    }
    
    T y1;

    if(type == ChebyshevT_Type)
    {
      y1 = x;
    }
    else if(type == LaguerreL_Type)
    {
      y1 = ef::one() - x;
    }
    else
    {
      y1 = x * static_cast<UINT32>(2u);
    }

    if(vp != static_cast<std::vector<T>*>(0u))
    {
      vp->push_back(y1);
    }

    if(n == static_cast<UINT32>(1u))
    {
      return y1;
    }

    T a = ef::two();
    T b = ef::zero();
    T c = ef::one();

    T yk;
  
    // Calculate higher orders using the recurrence relation.
    // The direction of stability is upward recurrence.
    for(INT32 k = static_cast<INT32>(2); k <= static_cast<INT32>(n); k++)
    {
      if(type == LaguerreL_Type)
      {
        a = -ef::one() / k;
        b =  ef::two() + a;
        c =  ef::one() + a;
      }
      else if(type == HermiteH_Type)
      {
        c = ef::two() * (k - ef::one());
      }
      
      yk = (((a * x) + b) * y1) - (c * y0);
      
      y0 = y1;
      y1 = yk;

      if(vp != static_cast<std::vector<T>*>(0u))
      {
        vp->push_back(yk);
      }
    }

    return yk;
  }
}

e_float ef::chebyshev_t(const INT32 n, const e_float& x)
{
  if(ef::isneg(x))
  {
    const bool b_negate = ((n % static_cast<INT32>(2)) != static_cast<INT32>(0));

    const e_float y = chebyshev_t(n, -x);

    return (!b_negate ? y : -y);
  }

  if(n < static_cast<INT32>(0))
  {
    const INT32 nn = static_cast<INT32>(-n);

    return chebyshev_t(nn, x);
  }
  else
  {
    return OrthogonalPoly_Series::OrthogonalPoly_Template(x, static_cast<UINT32>(n), OrthogonalPoly_Series::ChebyshevT_Type);
  }
}

e_float ef::chebyshev_u(const INT32 n, const e_float& x)
{
  if(ef::isneg(x))
  {
    const bool b_negate = ((n % static_cast<INT32>(2)) != static_cast<INT32>(0));

    const e_float y = chebyshev_u(n, -x);

    return (!b_negate ? y : -y);
  }

  if(n < static_cast<INT32>(0))
  {
    if(n == static_cast<INT32>(-2))
    {
      return ef::one();
    }
    else if(n == static_cast<INT32>(-1))
    {
      return ef::zero();
    }
    else
    {
      const INT32 n_minus_two = static_cast<INT32>(static_cast<INT32>(-n) - static_cast<INT32>(2));

      return -chebyshev_u(n_minus_two, x);
    }
  }
  else
  {
    return OrthogonalPoly_Series::OrthogonalPoly_Template(x, static_cast<UINT32>(n), OrthogonalPoly_Series::ChebyshevU_Type);
  }
}

e_float ef::hermite(const INT32 n, const e_float& x)
{
  if(n < static_cast<INT32>(0))
  {
    return ef::hermite(e_float(n), x);
  }
  else
  {
    return OrthogonalPoly_Series::OrthogonalPoly_Template(x, static_cast<UINT32>(n), OrthogonalPoly_Series::HermiteH_Type);
  }
}

e_float ef::laguerre(const INT32 n, const e_float& x)
{
  if(n < static_cast<INT32>(0))
  {
    // Use the representation of the Laguerre function in terms of the
    // hypergeometric 1F1 function: L_v(x) = 1F1(-v, 1, x).
    return ef::hyperg_1f1(e_float(static_cast<INT32>(-n)), ef::one(), x);
  }
  else if(n == static_cast<INT32>(0))
  {
    return ef::one();
  }
  else if(n == static_cast<INT32>(1))
  {
    return ef::one() - x;
  }

  if(ef::isneg(x))
  {
    // Use the representation of the Laguerre function in terms of the
    // hypergeometric 1F1 function: L_v(x) = 1F1(-v, 1, x).
    return ef::hyperg_1f1(e_float(static_cast<INT32>(-n)), ef::one(), x);
  }
  else
  {
    return OrthogonalPoly_Series::OrthogonalPoly_Template(x, static_cast<UINT32>(n), OrthogonalPoly_Series::LaguerreL_Type);
  }
}

ef_complex efz::chebyshev_t(const INT32 n, const ef_complex& z)
{
  if(n < static_cast<INT32>(0))
  {
    const INT32 nn = static_cast<INT32>(-n);

    return chebyshev_t(nn, z);
  }
  else
  {
    return OrthogonalPoly_Series::OrthogonalPoly_Template(z, static_cast<UINT32>(n), OrthogonalPoly_Series::ChebyshevT_Type);
  }
}

ef_complex efz::chebyshev_u(const INT32 n, const ef_complex& z)
{
  if(n < static_cast<INT32>(0))
  {
    if(n == static_cast<INT32>(-2))
    {
      return ef::one();
    }
    else if(n == static_cast<INT32>(-1))
    {
      return ef::zero();
    }
    else
    {
      const INT32 n_minus_two = static_cast<INT32>(static_cast<INT32>(-n) - static_cast<INT32>(2));

      return -chebyshev_u(n_minus_two, z);
    }
  }
  else
  {
    return OrthogonalPoly_Series::OrthogonalPoly_Template(z, static_cast<UINT32>(n), OrthogonalPoly_Series::ChebyshevU_Type);
  }
}

ef_complex efz::hermite(const INT32 n, const ef_complex& z)
{
  if(n < static_cast<INT32>(0))
  {
    // TBD: Negative order not yet implemented for complex Hermite polynomials.
    return ef::zero();
  }
  else
  {
    return OrthogonalPoly_Series::OrthogonalPoly_Template(z, static_cast<UINT32>(n), OrthogonalPoly_Series::HermiteH_Type);
  }
}

ef_complex efz::laguerre(const INT32 n, const ef_complex& z)
{
  if(ef::isneg(z))
  {
    // TBD: Negative argument not yet implemented for complex Laguerre polynomials.
    return ef::zero();
  }


  if(n < static_cast<INT32>(0))
  {
    // TBD: Negative order not yet implemented for complex Laguerre polynomials.
    return ef::zero();
  }
  else if(n == static_cast<INT32>(0))
  {
    return ef::one();
  }
  else if(n == static_cast<INT32>(1))
  {
    return ef::one() - z;
  }

  return OrthogonalPoly_Series::OrthogonalPoly_Template(z, static_cast<UINT32>(n), OrthogonalPoly_Series::LaguerreL_Type);
}

// NOCOVER_BLK_BEG
e_float ef::chebyshev_t(const UINT32 n, const e_float& x, std::vector<e_float>* const vp) { return OrthogonalPoly_Series::OrthogonalPoly_Template(x, static_cast<INT32>(n), OrthogonalPoly_Series::ChebyshevT_Type, vp); }
e_float ef::chebyshev_u(const UINT32 n, const e_float& x, std::vector<e_float>* const vp) { return OrthogonalPoly_Series::OrthogonalPoly_Template(x, static_cast<INT32>(n), OrthogonalPoly_Series::ChebyshevU_Type, vp); }
e_float ef::hermite    (const UINT32 n, const e_float& x, std::vector<e_float>* const vp) { return OrthogonalPoly_Series::OrthogonalPoly_Template(x, static_cast<INT32>(n), OrthogonalPoly_Series::HermiteH_Type,   vp); }
e_float ef::laguerre   (const UINT32 n, const e_float& x, std::vector<e_float>* const vp) { return OrthogonalPoly_Series::OrthogonalPoly_Template(x, static_cast<INT32>(n), OrthogonalPoly_Series::LaguerreL_Type,  vp); }

ef_complex efz::chebyshev_t(const UINT32 n, const ef_complex& z, std::vector<ef_complex>* const vp) { return OrthogonalPoly_Series::OrthogonalPoly_Template(z, static_cast<INT32>(n), OrthogonalPoly_Series::ChebyshevT_Type, vp); }
ef_complex efz::chebyshev_u(const UINT32 n, const ef_complex& z, std::vector<ef_complex>* const vp) { return OrthogonalPoly_Series::OrthogonalPoly_Template(z, static_cast<INT32>(n), OrthogonalPoly_Series::ChebyshevU_Type, vp); }
ef_complex efz::hermite    (const UINT32 n, const ef_complex& z, std::vector<ef_complex>* const vp) { return OrthogonalPoly_Series::OrthogonalPoly_Template(z, static_cast<INT32>(n), OrthogonalPoly_Series::HermiteH_Type,   vp); }
ef_complex efz::laguerre   (const UINT32 n, const ef_complex& z, std::vector<ef_complex>* const vp) { return OrthogonalPoly_Series::OrthogonalPoly_Template(z, static_cast<INT32>(n), OrthogonalPoly_Series::LaguerreL_Type,  vp); }
// NOCOVER_BLK_END
