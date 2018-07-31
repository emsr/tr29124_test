/*
$HOME/bin/bin/g++ -std=c++17 -o test_static_root test_static_root.cpp
*/

#include <limits>
#include <iostream>
#include <iomanip>

template<typename _FloatTp>
  constexpr _FloatTp
  __static_sqrt(_FloatTp __x)
  {
    if (__x > _FloatTp{1})
      return __x * __static_sqrt(_FloatTp{1} / __x);
    else
      {
	auto __lo = __x;
	auto __hi = _FloatTp{1};
	auto __mid = (__lo + __hi) / _FloatTp{2};
	for (int __i = 0; __i < std::numeric_limits<_FloatTp>::digits; ++__i)
	  {
	    if (__mid * __mid > __x)
	      __hi = __mid;
	    else
	      __lo = __mid;
	    __mid = (__lo + __hi) / _FloatTp{2};
	  }
	return __mid;
      }
  }

template<typename _FloatTp>
  constexpr _FloatTp
  __static_cbrt(_FloatTp __x)
  {
    if (__x > _FloatTp{1})
      return _FloatTp{1} / __static_cbrt(_FloatTp{1} / __x);
    else
      {
	auto __lo = __x;
	auto __hi = _FloatTp{1};
	auto __mid = (__lo + __hi) / _FloatTp{2};
	for (int __i = 0; __i < std::numeric_limits<_FloatTp>::digits; ++__i)
	  {
	    if (__mid * __mid * __mid > __x)
	      __hi = __mid;
	    else
	      __lo = __mid;
	    __mid = (__lo + __hi) / _FloatTp{2};
	  }
	return __mid;
      }
  }

int
main()
{
  std::cout << __static_sqrt(std::numeric_limits<double>::epsilon()) << '\n';
  std::cout << __static_cbrt(std::numeric_limits<double>::epsilon()) << '\n';
}
