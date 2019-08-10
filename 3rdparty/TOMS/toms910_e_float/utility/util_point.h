
#ifndef _UTIL_POINT_2009_10_27_H_
  #define _UTIL_POINT_2009_10_27_H_

  namespace Util
  {
    template<typename T1, typename T2 = T1> struct point
    {
      T1 x;
      T2 y;

      point(const T1& X = T1(), const T2& Y = T2()) : x(X), y(Y) { }
    };

    template<typename T1, typename T2> bool inline operator<(const point<T1, T2>& left, const point<T1, T2>& right)
    {
      return (left.x < right.x);
    }
  }

#endif // _UTIL_POINT_2009_10_27_H_
