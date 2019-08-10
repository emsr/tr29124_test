
#ifndef _LEGENDRE_X_2009_03_29_H_
  #define _LEGENDRE_X_2009_03_29_H_

  #include <utility/util_noncopyable.h>

  class e_float;

  namespace Legendre_Series
  {
    class LegendreX : private Util::noncopyable
    {
    protected:

      LegendreX() { }

      virtual ~LegendreX() { } // NOCOVER_LINE
    };
  }

#endif // _LEGENDRE_X_2009_03_29_H_
