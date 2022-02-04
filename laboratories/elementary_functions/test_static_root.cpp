/**
 *
 */

#include <limits>
#include <iostream>
#include <iomanip>

template<typename Tp>
  constexpr Tp
  static_sqrt(Tp x)
  {
    if (x < Tp{0})
      return std::numeric_limits<Tp>::quiet_NaN();
    else if (x > Tp{1})
      return x * static_sqrt(Tp{1} / x);
    else
      {
	auto lo = x;
	auto hi = Tp{1};
	auto mid = (lo + hi) / Tp{2};
	for (int i = 0; i < std::numeric_limits<Tp>::digits; ++i)
	  {
	    if (mid * mid > x)
	      hi = mid;
	    else
	      lo = mid;
	    mid = (lo + hi) / Tp{2};
	  }
	return mid;
      }
  }

template<typename Tp>
  constexpr Tp
  static_cbrt(Tp x)
  {
    if (x < Tp{0})
      return std::numeric_limits<Tp>::quiet_NaN();
    else if (x > Tp{1})
      return Tp{1} / static_cbrt(Tp{1} / x);
    else
      {
	auto lo = x;
	auto hi = Tp{1};
	auto mid = (lo + hi) / Tp{2};
	for (int i = 0; i < std::numeric_limits<Tp>::digits; ++i)
	  {
	    if (mid * mid * mid > x)
	      hi = mid;
	    else
	      lo = mid;
	    mid = (lo + hi) / Tp{2};
	  }
	return mid;
      }
  }

template<typename Tp>
  constexpr Tp
  static_powi(Tp x, unsigned int r)
  {
    // There are better ways...
    Tp pp = Tp{1};
    for (unsigned int i = 0; i < r; ++i)
      pp *= x;
    return pp;
  }

template<typename Tp>
  constexpr Tp
  static_root(unsigned int r, Tp x)
  {
    if (x < Tp{0})
      return std::numeric_limits<Tp>::quiet_NaN();
    else if (x > Tp{1})
      return Tp{1} / static_cbrt(Tp{1} / x);
    else
      {
	auto lo = x;
	auto hi = Tp{1};
	auto mid = (lo + hi) / Tp{2};
	for (int i = 0; i < std::numeric_limits<Tp>::digits; ++i)
	  {
	    if (static_powi(mid, r) > x)
	      hi = mid;
	    else
	      lo = mid;
	    mid = (lo + hi) / Tp{2};
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
