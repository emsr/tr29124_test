
#ifndef _LEGENDRE_XNM_2009_03_29_H_
  #define _LEGENDRE_XNM_2009_03_29_H_

  #include <functions/hypergeometric/legendre_x.h>

  namespace Legendre_Series
  {
    class LegendreXnm : public LegendreX
    {
    protected:

      LegendreXnm() { }

      virtual ~LegendreXnm() { } // NOCOVER_LINE

    private:

      virtual e_float AtIdenticallyZero        (const INT32 n, const INT32 m) const = 0;
      virtual e_float AtIdenticallyOne         (const INT32 n, const INT32 m) const = 0;
      virtual e_float AtReflectNegativeDegree  (const INT32 n, const INT32 m, const e_float& x) const = 0;
              e_float AtReflectNegativeOrder   (const INT32 n, const INT32 m, const e_float& x) const;
      virtual e_float AtReflectNegativeArgument(const INT32 n, const INT32 m, const e_float& x) const = 0;

      virtual bool AtResultForDegreeAndOrderIsZero(const INT32 n, const INT32 m) const = 0;

      virtual e_float L_00(const e_float& x) const = 0;
      virtual e_float L_10(const e_float& x) const = 0;
      virtual e_float L_01(const e_float& x) const = 0;
      virtual e_float L_11(const e_float& x) const = 0;

      virtual e_float Legendre_n(const INT32 n, const e_float& x) const = 0;

      virtual e_float MyRecursion(const INT32 n, const INT32 m, const e_float& x) const = 0;

    public:

      e_float MyLegendre(const INT32 n, const INT32 m, const e_float& x) const;
    };
  }

#endif // _LEGENDRE_XNM_2009_03_29_H_
