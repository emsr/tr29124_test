/**
 *
 */

#include <limits>
#include <iostream>
#include <iomanip>

template<typename _Tp>
  constexpr _Tp
  static_sqrt(_Tp x)
  {
    if (x < _Tp{0})
      return std::numeric_limits<_Tp>::quiet_NaN();
    else if (x > _Tp{1})
      return x * static_sqrt(_Tp{1} / x);
    else
      {
	auto lo = x;
	auto hi = _Tp{1};
	auto mid = (lo + hi) / _Tp{2};
	for (int i = 0; i < std::numeric_limits<_Tp>::digits; ++i)
	  {
	    if (mid * mid > x)
	      hi = mid;
	    else
	      lo = mid;
	    mid = (lo + hi) / _Tp{2};
	  }
	return mid;
      }
  }

template<typename _Tp>
  constexpr _Tp
  static_cbrt(_Tp x)
  {
    if (x < _Tp{0})
      return std::numeric_limits<_Tp>::quiet_NaN();
    else if (x > _Tp{1})
      return _Tp{1} / static_cbrt(_Tp{1} / x);
    else
      {
	auto lo = x;
	auto hi = _Tp{1};
	auto mid = (lo + hi) / _Tp{2};
	for (int i = 0; i < std::numeric_limits<_Tp>::digits; ++i)
	  {
	    if (mid * mid * mid > x)
	      hi = mid;
	    else
	      lo = mid;
	    mid = (lo + hi) / _Tp{2};
	  }
	return mid;
      }
  }

template<typename _Tp>
  constexpr _Tp
  static_powi(_Tp x, unsigned int r)
  {
    // There are better ways...
    _Tp pp = _Tp{1};
    for (unsigned int i = 0; i < r; ++i)
      pp *= x;
    return pp;
  }

template<typename _Tp>
  constexpr _Tp
  static_root(unsigned int r, _Tp x)
  {
    if (x < _Tp{0})
      return std::numeric_limits<_Tp>::quiet_NaN();
    else if (x > _Tp{1})
      return _Tp{1} / static_cbrt(_Tp{1} / x);
    else
      {
	auto lo = x;
	auto hi = _Tp{1};
	auto mid = (lo + hi) / _Tp{2};
	for (int i = 0; i < std::numeric_limits<_Tp>::digits; ++i)
	  {
	    if (static_powi(mid, r) > x)
	      hi = mid;
	    else
	      lo = mid;
	    mid = (lo + hi) / _Tp{2};
	  }
	return mid;
      }
  }

int
main()
{
  std::cout << static_sqrt(std::numeric_limits<double>::epsilon()) << '\n';
  std::cout << static_cbrt(std::numeric_limits<double>::epsilon()) << '\n';
  std::cout << static_root(6, std::numeric_limits<double>::epsilon()) << '\n';
}
