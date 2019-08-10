
#ifndef _LEGENDRE_XNU_2009_03_29_H_
  #define _LEGENDRE_XNU_2009_03_29_H_

  #include <functions/hypergeometric/legendre_x.h>

  namespace Legendre_Series
  {
    class LegendreXnu : public LegendreX
    {
    protected:

      LegendreXnu() { }

      virtual ~LegendreXnu() { } // NOCOVER_LINE

    private:

      // NOCOVER_BLK_BEG
      virtual e_float AtReflectNegativeDegree  (const INT32 n, const e_float& u, const e_float& x) const = 0;
      virtual e_float AtReflectNegativeOrder   (const INT32 n, const e_float& u, const e_float& x) const { static_cast<void>(n); static_cast<void>(u); static_cast<void>(x); return ef::zero(); }
      virtual e_float AtReflectNegativeArgument(const INT32 n, const e_float& u, const e_float& x) const = 0;
      virtual e_float AtOnePlus                (const INT32 n, const e_float& u, const e_float& x) const = 0;

      virtual e_float Legendre_nm(const INT32 n, const INT32 m, const e_float& x) const = 0;
      // NOCOVER_BLK_END

      virtual bool NeedsReflectNegativeOrder(const e_float& u) const { return ef::isneg(u); }

    public:

      e_float MyLegendre(const INT32 n, const e_float& u, const e_float& x) const;
    };
  }

#endif // _LEGENDRE_XNU_2009_03_29_H_
