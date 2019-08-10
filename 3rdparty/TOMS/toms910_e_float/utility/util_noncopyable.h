
#ifndef _UTIL_NONCOPYABLE_2009_03_30_H_
  #define _UTIL_NONCOPYABLE_2009_03_30_H_

  namespace Util
  {
    class noncopyable
    {
    protected:

      noncopyable() {}
      virtual ~noncopyable() {}

    private:  // emphasize the following members are private

      noncopyable(const noncopyable&);
      const noncopyable& operator=(const noncopyable&);
    };
  }

#endif // _UTIL_NONCOPYABLE_2009_03_30_H_
