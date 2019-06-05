/**
 *
 */

#include <limits>
#include <iostream>
#include <iomanip>

template<typename _Tp>
  constexpr _Tp
  __static_sqrt(_Tp __x)
  {
    if (__x < _Tp{0})
      return std::numeric_limits<_Tp>::quiet_NaN();
    else if (__x > _Tp{1})
      return __x * __static_sqrt(_Tp{1} / __x);
    else
      {
	auto __lo = __x;
	auto __hi = _Tp{1};
	auto __mid = (__lo + __hi) / _Tp{2};
	for (int __i = 0; __i < std::numeric_limits<_Tp>::digits; ++__i)
	  {
	    if (__mid * __mid > __x)
	      __hi = __mid;
	    else
	      __lo = __mid;
	    __mid = (__lo + __hi) / _Tp{2};
	  }
	return __mid;
      }
  }

template<typename _Tp>
  constexpr _Tp
  __static_cbrt(_Tp __x)
  {
    if (__x < _Tp{0})
      return std::numeric_limits<_Tp>::quiet_NaN();
    else if (__x > _Tp{1})
      return _Tp{1} / __static_cbrt(_Tp{1} / __x);
    else
      {
	auto __lo = __x;
	auto __hi = _Tp{1};
	auto __mid = (__lo + __hi) / _Tp{2};
	for (int __i = 0; __i < std::numeric_limits<_Tp>::digits; ++__i)
	  {
	    if (__mid * __mid * __mid > __x)
	      __hi = __mid;
	    else
	      __lo = __mid;
	    __mid = (__lo + __hi) / _Tp{2};
	  }
	return __mid;
      }
  }

template<typename _Tp>
  constexpr _Tp
  __static_powi(_Tp __x, unsigned int __r)
  {
    // There are better ways...
    _Tp __pp = _Tp{1};
    for (unsigned int __i = 0; __i < __r; ++__i)
      __pp *= __x;
    return __pp;
  }

template<typename _Tp>
  constexpr _Tp
  __static_root(unsigned int __r, _Tp __x)
  {
    if (__x < _Tp{0})
      return std::numeric_limits<_Tp>::quiet_NaN();
    else if (__x > _Tp{1})
      return _Tp{1} / __static_cbrt(_Tp{1} / __x);
    else
      {
	auto __lo = __x;
	auto __hi = _Tp{1};
	auto __mid = (__lo + __hi) / _Tp{2};
	for (int __i = 0; __i < std::numeric_limits<_Tp>::digits; ++__i)
	  {
	    if (__static_powi(__mid, __r) > __x)
	      __hi = __mid;
	    else
	      __lo = __mid;
	    __mid = (__lo + __hi) / _Tp{2};
	  }
	return __mid;
      }
  }

int
main()
{
  std::cout << __static_sqrt(std::numeric_limits<double>::epsilon()) << '\n';
  std::cout << __static_cbrt(std::numeric_limits<double>::epsilon()) << '\n';
  std::cout << __static_root(6, std::numeric_limits<double>::epsilon()) << '\n';
}
