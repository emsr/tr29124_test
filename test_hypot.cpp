// $HOME/bin_specfun/bin/g++ -std=gnu++14 -g -Wall -Wextra -o test_hypot test_hypot.cpp -L$HOME/bin/lib64 -lquadmath

// LD_LIBRARY_PATH=$HOME/bin_specfun/lib64:$LD_LIBRARY_PATH ./test_hypot > test_hypot.txt

// g++ -std=gnu++14 -Wall -Wextra -DNO_LOGBQ -I. -o test_hypot test_hypot.cpp -lquadmath

// ./test_hypot > test_hypot.txt

#include <limits>
#include <iostream>
#include <iomanip>
#include <algorithm>
#include <cmath>

namespace std
{

#define __cpp_lib_hypot 201603L

  template<typename _Tp>
    _Tp
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

int
main()
{
  auto m123 = std::hypot(1.0, 2.0, 3.0);
  std::cout << "m123 = " << m123 << '\n';
  auto m1big = std::hypot(1.0e300, 2.0, 3.0);
  std::cout << "m1big = " << m1big << '\n';
  auto m2big = std::hypot(1.0e300, 2.0e300, 3.0);
  std::cout << "m2big = " << m2big << '\n';
  auto m3big = std::hypot(1.0e300, 2.0e300, 3.0e300);
  std::cout << "m3big = " << m3big << '\n';
  auto m1small = std::hypot(1.0e-300, 2.0, 3.0);
  std::cout << "m1small = " << m1small << '\n';
  auto m2small = std::hypot(1.0e-300, 2.0e-300, 3.0);
  std::cout << "m2small = " << m2small << '\n';
  auto m3small = std::hypot(1.0e-300, 2.0e-300, 3.0e-300);
  std::cout << "m3small = " << m3small << '\n';
  auto m3zero = std::hypot(0.0, 0.0, 0.0);
  std::cout << "m3zero = " << m3zero << '\n';
}
