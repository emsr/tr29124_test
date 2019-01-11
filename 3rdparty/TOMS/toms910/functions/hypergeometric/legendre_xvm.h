
#ifndef _LEGENDRE_XVM_2009_03_29_H_
  #define _LEGENDRE_XVM_2009_03_29_H_

  #include <functions/hypergeometric/legendre_x.h>

  namespace Legendre_Series
  {
    class LegendreXvm : public LegendreX
    {
    protected:

      LegendreXvm() { }

      virtual ~LegendreXvm() { } // NOCOVER_LINE

    private:

      // NOCOVER_BLK_BEG
      virtual e_float AtIdenticallyZero        (const e_float& v, const INT32 m) const = 0;
      virtual e_float AtReflectNegativeDegree  (const e_float& v, const INT32 m, const e_float& x) const = 0;
      virtual e_float AtReflectNegativeOrder   (const e_float& v, const INT32 m, const e_float& x) const = 0;
      virtual e_float AtReflectNegativeArgument(const e_float& v, const INT32 m, const e_float& x) const = 0;
      virtual e_float AtOnePlus                (const e_float& v, const INT32 m, const e_float& x) const = 0;

      virtual bool NeedsRecurOfM(const e_float& x) const = 0;

      virtual e_float Legendre_nm(const INT32 n, const INT32 m, const e_float& x) const = 0;
      virtual e_float Legendre_v (const e_float& v, const e_float& x) const = 0;
      // NOCOVER_BLK_END

    public:

      e_float MyLegendre(const e_float& v, const INT32 m, const e_float& x) const;
    };
  }

#endif // _LEGENDRE_XVM_2009_03_29_H_
