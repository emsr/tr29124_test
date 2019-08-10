#ifndef _UTIL_TIMER_2010_01_26_H_
  #define _UTIL_TIMER_2010_01_26_H_

  #include <utility/util_noncopyable.h>

  namespace Util
  {
    struct timer : private noncopyable
    {
    private:

      const volatile double offset;
      const volatile double start;

      static double get_sec(void);

    public:

      timer(const double ofs = static_cast<double>(0.0)) : offset(ofs),
                                                           start (get_sec()) { }

      virtual ~timer() { }

      double elapsed(void) const;
    };
  }

#endif // _UTIL_TIMER_2010_01_26_H_
