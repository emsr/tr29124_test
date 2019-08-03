#ifndef MY_HYPOT_H
#define MY_HYPOT_H 1

#if __cplusplus > 201402L

namespace std _GLIBCXX_VISIBILITY(default)
{
_GLIBCXX_BEGIN_NAMESPACE_VERSION

  // [c.math.hypot3], three-dimensional hypotenuse
#define __cpp_lib_hypot 201603

#include <limits>

  // Avoid including all of <algorithm>
  template<typename _Tp>
    constexpr _Tp
    __fmax3(_Tp __x, _Tp __y, _Tp __z)
    { return std::fmax(std::fmax(__x, __y), std::fmax(__y, __z)); }

  template<typename _Tp>
    constexpr _Tp
    __hypot3(_Tp __x, _Tp __y, _Tp __z)
    {
      if (std::isnan(__x)
       || std::isnan(__y)
       || std::isnan(__z))
	return std::numeric_limits<_Tp>::quiet_NaN();
      else
	{
	  __x = std::abs(__x);
	  __y = std::abs(__y);
	  __z = std::abs(__z);
	  const auto __amax = __fmax3(__x, __y, __z);
	  if (__amax == _Tp{0})
	    return _Tp{0};
	  else if (std::isinf(__amax))
	    return std::numeric_limits<_Tp>::infinity();
	  else
	    {
	      __x /= __amax;
	      __y /= __amax;
	      __z /= __amax;
	      return __amax * std::sqrt(__x * __x + __y * __y + __z * __z);
            }
	}
    }

  /**
   * Return the three-dimensional hypoteneuse @f$ \sqrt{x^2 + y^2 + z^2} @f$
   * for @c float arguments x, y, and z.
   */
  constexpr inline float
  hypot(float __x, float __y, float __z)
  { return std::__hypot3<float>(__x, __y, __z); }

  /**
   * Return the three-dimensional hypoteneuse @f$ \sqrt{x^2 + y^2 + z^2} @f$
   * for <tt>long double</tt> arguments x, y, and z.
   */
  constexpr inline long double
  hypot(long double __x, long double __y, long double __z)
  { return std::__hypot3<long double>(__x, __y, __z); }

  /**
   * Return the three-dimensional hypoteneuse @f$ \sqrt{x^2 + y^2 + z^2} @f$
   * for real arguments x, y, and z.
   */
  template<typename _Tp, typename _Up, typename _Vp>
    constexpr typename __gnu_cxx::__promote_3<_Tp, _Up, _Vp>::__type
    hypot(_Tp __x, _Up __y, _Vp __z)
    {
      using __type = typename __gnu_cxx::__promote_3<_Tp, _Up, _Vp>::__type;
      return std::__hypot3<__type>(__x, __y, __z);
    }

_GLIBCXX_END_NAMESPACE_VERSION
} // namespace

#endif // C++17

#endif // MY_HYPOT_H
