
#ifndef _LEGENDRE_XV_2009_03_29_H_
  #define _LEGENDRE_XV_2009_03_29_H_

  #include <functions/hypergeometric/legendre_x.h>

  namespace Legendre_Series
  {
    class LegendreXv : public LegendreX
    {
    protected:

      LegendreXv() { }

      virtual ~LegendreXv() { } // NOCOVER_LINE

    private:

      virtual e_float AtIdenticallyZero      (const e_float& v) const = 0;
      virtual e_float AtReflectNegativeDegree(const e_float& v, const e_float& x) const = 0;
      virtual e_float AtOnePlus              (const e_float& v, const e_float& x) const = 0;
      virtual e_float AtOneMinus             (const e_float& v, const e_float& x) const = 0;

      virtual e_float Legendre_n(const INT32 n, const e_float& x) const = 0;

    public:

      e_float MyLegendre(const e_float& v, const e_float& x) const;
    };
  }

#endif // _LEGENDRE_XV_2009_03_29_H_
