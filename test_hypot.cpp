/*
$HOME/bin_tr29124/bin/g++ -std=gnu++14 -g -Wall -Wextra -o test_hypot test_hypot.cpp -L$HOME/bin/lib64 -lquadmath
LD_LIBRARY_PATH=$HOME/bin_tr29124/lib64:$LD_LIBRARY_PATH ./test_hypot > test_hypot.txt

g++ -std=gnu++14 -Wall -Wextra -DNO_LOGBQ -I. -o test_hypot test_hypot.cpp -lquadmath
./test_hypot > test_hypot.txt
*/

#include <limits>
#include <iostream>
#include <iomanip>
#include <algorithm>
#include <cmath>

#define __cpp_lib_hypot 201603L

namespace std
{

  template<typename _Tp>
    constexpr _Tp
    hypot(_Tp __x, _Tp __y, _Tp __z)
    {
      auto __abs_max = [](_Tp __a, _Tp __b)
			-> bool
			{ return std::abs(__a) < std::abs(__b); };

      auto __amax = std::max({__x, __y, __z}, __abs_max);
      if (__amax == _Tp{0})
	return _Tp{0};
      else
	{
	  __x /= __amax;
	  __y /= __amax;
	  __z /= __amax;
	  return __amax * std::sqrt(__x * __x + __y * __y + __z * __z);
        }
    }

} // namespace std

template<typename _Tp>
  void
  test()
  {
    constexpr auto tiny = _Tp{8} * std::numeric_limits<double>::lowest();
    constexpr auto huge = _Tp{0.125L} * std::numeric_limits<double>::max();

    constexpr auto m123 = std::hypot(_Tp{1}, _Tp{2}, _Tp{3});
    constexpr auto m1huge = std::hypot(huge, _Tp{2}, _Tp{3});
    constexpr auto m2huge = std::hypot(huge, _Tp{2} * huge, _Tp{3});
    constexpr auto m3huge = std::hypot(huge, _Tp{2} * huge, _Tp{3} * huge);
    constexpr auto m1tiny = std::hypot(tiny, _Tp{2}, _Tp{3});
    constexpr auto m2tiny = std::hypot(tiny, _Tp{2} * tiny, _Tp{3});
    constexpr auto m3tiny = std::hypot(tiny, _Tp{2} * tiny, _Tp{3} * tiny);
  }

int
main()
{
  std::cout.precision(std::numeric_limits<double>::digits10);
  auto w = 8 + std::cout.precision();

  auto m123 = std::hypot(1.0, 2.0, 3.0);
  std::cout << "m123    = " << std::setw(w) << m123 << '\n';

  auto m1big = std::hypot(1.0e300, 2.0, 3.0);
  std::cout << "m1big   = " << std::setw(w) << m1big << '\n';

  auto m2big = std::hypot(1.0e300, 2.0e300, 3.0);
  std::cout << "m2big   = " << std::setw(w) << m2big << '\n';

  auto m3big = std::hypot(1.0e300, 2.0e300, 3.0e300);
  std::cout << "m3big   = " << std::setw(w) << m3big << '\n';

  auto m1small = std::hypot(1.0e-300, 2.0, 3.0);
  std::cout << "m1small = " << std::setw(w) << m1small << '\n';

  auto m2small = std::hypot(1.0e-300, 2.0e-300, 3.0);
  std::cout << "m2small = " << std::setw(w) << m2small << '\n';

  auto m3small = std::hypot(1.0e-300, 2.0e-300, 3.0e-300);
  std::cout << "m3small = " << std::setw(w) << m3small << '\n';

  auto m3zero = std::hypot(0.0, 0.0, 0.0);
  std::cout << "m3zero  = " << std::setw(w) << m3zero << '\n';
}

/*
m123    =        3.74165738677394
m1big   =                  1e+300
m2big   =   2.23606797749979e+300
m3big   =   3.74165738677394e+300
m1small =        3.60555127546399
m2small =                       3
m3small =   3.74165738677394e-300
m3zero  =                       0
*/
