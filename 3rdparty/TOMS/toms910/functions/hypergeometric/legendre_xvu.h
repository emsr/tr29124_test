
#ifndef _LEGENDRE_XVU_2009_03_29_H_
  #define _LEGENDRE_XVU_2009_03_29_H_

  #include <functions/hypergeometric/legendre_x.h>

  namespace Legendre_Series
  {
    class LegendreXvu : public LegendreX
    {
    protected:

      LegendreXvu() { }

      virtual ~LegendreXvu() { } // NOCOVER_LINE

    private:

      // NOCOVER_BLK_BEG
      virtual e_float AtReflectNegativeDegree  (const e_float& v, const e_float& u, const e_float& x) const = 0;
      virtual e_float AtReflectNegativeOrder   (const e_float& v, const e_float& u, const e_float& x) const { static_cast<void>(v); static_cast<void>(u); static_cast<void>(x); return ef::zero(); }
      virtual e_float AtReflectNegativeArgument(const e_float& v, const e_float& u, const e_float& x) const = 0;
      virtual e_float AtOnePlus                (const e_float& v, const e_float& u, const e_float& x) const = 0;

      virtual e_float Legendre_vm(const e_float& v, const INT32    m, const e_float& x) const = 0;
      virtual e_float Legendre_nu(const INT32    n, const e_float& u, const e_float& x) const = 0;
      // NOCOVER_BLK_END

      virtual bool NeedsReflectNegativeOrder(const e_float& u) const { return ef::isneg(u); }

    public:

      e_float MyLegendre(const e_float& v, const e_float& u, const e_float& x) const;
    };
  }

#endif // _LEGENDRE_XVU_2009_03_29_H_
