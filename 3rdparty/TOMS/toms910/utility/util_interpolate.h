
#ifndef _UTIL_INTERPOLATE_2009_10_27_H_
  #define _UTIL_INTERPOLATE_2009_10_27_H_

  #include <vector>
  #include <algorithm>

  #include <utility/util_point.h>

  namespace Util
  {
    template<typename T1, typename T2 = T1> struct linear_interpolate
    {
      static T2 interpolate(const T1& x, const std::vector<Util::point<T1, T2> >& points)
      {
        if(points.empty())
        {
          return T2(0);
        }
        else if(x <= points.front().x || points.size() == static_cast<std::size_t>(1u))
        {
          return points.front().y;
        }
        else if(x >= points.back().x)
        {
          return points.back().y;
        }
        else
        {
          const Util::point<T1, T2> x_find(x);

          const typename std::vector<Util::point<T1, T2> >::const_iterator it =
            std::lower_bound(points.begin(), points.end(), x_find);

          const T1 xn            = (it - 1u)->x;
          const T1 xnp1_minus_xn = it->x - xn;
          const T1 delta_x       = x - xn;
          const T2 yn            = (it - 1u)->y;
          const T2 ynp1_minus_yn = it->y - yn;

          return T2(yn + T2((delta_x * ynp1_minus_yn) / xnp1_minus_xn));
        }
      }
    };
  }

#endif // _UTIL_INTERPOLATE_2009_10_27_H_
