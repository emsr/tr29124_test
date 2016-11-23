/*
$HOME/bin_tr29124/bin/g++ -std=gnu++17 -g -Wall -Wextra -Wno-psabi -I. -o test_pochhammer test_pochhammer.cpp wrap_boost.cpp -lquadmath
./test_pochhammer > test_pochhammer.txt

$HOME/bin/bin/g++ -std=gnu++17 -g -I. -o test_pochhammer test_pochhammer.cpp wrap_boost.cpp
./test_pochhammer > test_pochhammer.txt
*/

#include <limits>
#include <vector>
#include <iostream>
#include <iomanip>
#include <ext/cmath>
#include "wrap_boost.h"

namespace __gnu_cxx
{

  template<typename _Tp>
    _Tp
    __pochhammer_upper_prod(int __a, int __n)
    {
      if (__a < __n)
	return _Tp{0};
      else
	{
	  auto __prod = 1;
	  for (int __k = 0; __k < __n; ++__k)
	    __prod *= __a--;
	  return __prod;
	}
    }

  template<typename _Tp>
    _Tp
    __pochhammer_prod(_Tp __a, int __n)
    {
      auto __prod = _Tp{1};
      for (int __k = 0; __k < __n; ++__k)
	__prod *= __a++;
      return __prod;
    }

  template<typename _Tp>
    _Tp
    __pochhammer_fake(_Tp __a, _Tp __x)
    {
      auto __n = int(std::nearbyint(__x));
      if (_Tp(__n) == __x)
	{
	  if (__n == 0)
	    return _Tp{1};
	  else
	    {
	      auto __m = int(std::nearbyint(__a));
	      if (int(__m) == __a)
		return __pochhammer_prod<_Tp>(__m, __n);
	      else
		return __pochhammer_prod(__a, __n);
	    }
	}
      else
	return std::__detail::__gamma(__a + __x)
	     / std::__detail::__gamma(__a);
    }

}

int
main()
{
  std::cout.precision(std::numeric_limits<double>::digits10);
  std::cout.flags(std::ios::showpoint);
  auto width = 8 + std::cout.precision();

  std::vector<double> xv{0.0, 1.0, 2.0, 5.0, 10.0, 20.0, 50.0, 100.0};

  // Test call (log form).
  std::cout << "\nLogarithmic form\n";
  for (int ia = 0; ia <= +500; ++ia)
    {
      auto a = ia * 0.01;
      std::cout << '\n';
      for (auto x : xv)
	{
	  auto pochg = __gnu_cxx::pochhammer(a, x);
	  auto pochb = beast::pochhammer(a, x);
	  std::cout << ' ' << std::setw(width) << a
		    << ' ' << std::setw(width) << x
		    << ' ' << std::setw(width) << pochg
		    << ' ' << std::setw(width) << pochb
		    << ' ' << std::setw(width) << pochg - pochb
		    << ' ' << std::setw(width) << (pochg - pochb) / pochb
		    << '\n';
	}
    }

  // Test naive product form.
  std::cout << "\nProduct form\n";
  for (int ia = 0; ia <= +500; ++ia)
    {
      auto a = ia * 0.01;
      std::cout << '\n';
      for (auto x : xv)
	{
	  auto pochg = __gnu_cxx::__pochhammer_fake(a, x);
	  auto pochb = beast::pochhammer(a, x);
	  std::cout << ' ' << std::setw(width) << a
		    << ' ' << std::setw(width) << x
		    << ' ' << std::setw(width) << pochg
		    << ' ' << std::setw(width) << pochb
		    << ' ' << std::setw(width) << pochg - pochb
		    << ' ' << std::setw(width) << (pochg - pochb) / pochb
		    << '\n';
	}
    }
}
