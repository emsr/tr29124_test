/*
$HOME/bin/bin/g++ -std=gnu++2a -g -Wall -Wextra -Wno-psabi -I. -o test_falling_factorial test_falling_factorial.cpp -lquadmath -Lwrappers/debug -lwrap_boost
LD_LIBRARY_PATH=wrappers/debug:$LD_LIBRARY_PATH ./test_falling_factorial > test_falling_factorial.txt

$HOME/bin/bin/g++ -std=gnu++2a -g -Wall -Wextra -I. -o test_falling_factorial test_falling_factorial.cpp -lquadmath -Lwrappers/debug -lwrap_boost -lwrap_gsl
PATH=wrappers/debug:$PATH ./test_falling_factorial > test_falling_factorial.txt
*/

#include <limits>
#include <vector>
#include <iostream>
#include <iomanip>
#include <cmath>

#include <wrap_boost.h>

namespace __gnu_cxx
{

  template<typename _Tp>
    _Tp
    __falling_factorial_prod(int __a, int __n)
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
    __falling_factorial_prod(_Tp __a, int __n)
    {
      auto __prod = _Tp{1};
      for (int __k = 0; __k < __n; ++__k)
	__prod *= __a--;
      return __prod;
    }

  template<typename _Tp>
    _Tp
    __falling_factorial_fake(_Tp __a, _Tp __x)
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
		return __falling_factorial_prod<_Tp>(__m, __n);
	      else
		return __falling_factorial_prod(__a, __n);
	    }
	}
      else
	return std::__detail::__gamma(__a + _Tp{1})
	     / std::__detail::__gamma(__a - __x + _Tp{1});
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
	  auto pochg = __gnu_cxx::falling_factorial(a, x);
	  auto pochb = beast::falling_factorial(a, x);
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
	  auto pochg = __gnu_cxx::__falling_factorial_fake(a, x);
	  auto pochb = beast::falling_factorial(a, x);
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
