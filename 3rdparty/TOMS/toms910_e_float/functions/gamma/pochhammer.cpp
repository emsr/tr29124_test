#include <functions/complex/e_float_complex.h>
#include <functions/elementary/elementary.h>
#include <functions/gamma/gamma.h>

namespace Pochhammer_Series
{
  template<typename T> static inline T Pochhammer_Template(const T& x, const UINT32 n)
  {
    using ef::real;
    using efz::real;
    using ef::gamma;
    using efz::gamma;

    if(n == static_cast<UINT32>(0u))
    {
      return ef::one();
    }
    else if(n == static_cast<UINT32>(1u))
    {
      return x;
    }
    else
    {
      if(n < static_cast<UINT32>(50u))
      {
        T val (x);
        T term(x);

        for(UINT32 i = static_cast<UINT32>(1u); i < n; i++)
        {
          val *= ++term;
        }

        return val;
      }
      else
      {
        return gamma(x + T(n)) / gamma(x);
      }
    }
  }
}

e_float ef::pochhammer(const e_float& x, const UINT32 n)
{
  return Pochhammer_Series::Pochhammer_Template<e_float>(x, n);
}

ef_complex efz::pochhammer(const ef_complex& x, const UINT32 n)
{
  return Pochhammer_Series::Pochhammer_Template<ef_complex>(x, n);
}

e_float ef::pochhammer(const e_float& x, const e_float& a)
{
  return ef::gamma(x + a) / ef::gamma(x);
}

ef_complex efz::pochhammer(const ef_complex& z, const ef_complex& a)
{
  return efz::gamma(z + a) / efz::gamma(z);
}

