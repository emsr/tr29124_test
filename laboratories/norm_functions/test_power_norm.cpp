/*
$HOME/bin_specfun/bin/g++ -std=c++2a -g -Wall -Wextra -Wno-psabi -I. -o test_power_norm test_power_norm.cpp
./test_power_norm > test_power_norm.txt

$HOME/bin/bin/g++ -std=c++2a -g -Wall -Wextra -Wno-psabi -I. -o test_power_norm test_power_norm.cpp
./test_power_norm > test_power_norm.txt
*/

#include <initializer_list>
#include <cmath>
#include <iostream>
#include <iomanip>
#include <limits>
#include <algorithm>

/**
 * 
 */
template<typename _Real, typename _Iter>
  _Real
  __impl_power_norm_0(_Iter __begin, _Iter __end)
  {
    std::size_t __n = 0;
    auto __prod = _Real{1};
    while (__begin != __end)
      {
	++__n;
	__prod *= *__begin;
	++__begin;
      }
    return std::pow(std::abs(__prod), _Real{1} / _Real(__n));
  }

/**
 * 
 */
template<typename _Real, typename _Iter>
  _Real
  __impl_power_norm(_Real __p, _Iter __begin, _Iter __end)
  {
    if (__begin == __end)
      return _Real{0};

    using _Tp = std::decay_t<decltype(*__begin)>;
    if (__p == _Real{0})
      return __impl_power_norm_0<_Tp>(__begin, __end);
    else if (std::isinf(__p))
      {
	if (__p < _Real{0})
	  return std::abs(*std::min_element(__begin, __end));
	else
	  return std::abs(*std::max_element(__begin, __end));
      }
    else
      {
	std::size_t __n = 0;
	auto __sum = _Real{0};
	while (__begin != __end)
	  {
	    ++__n;
	    __sum += std::pow(std::abs(*__begin), __p);
	    ++__begin;
	  }
	return std::pow(__sum / _Real(__n), _Real{1} / _Real(__p));
      }
  }

/**
 * 
 */
template<typename _Real>
  _Real
  power_norm(_Real __p, std::initializer_list<_Real> __val)
  {
    using _Tp = std::decay_t<decltype(*__val.begin())>;
    return __impl_power_norm<_Tp>(__p, __val.begin(), __val.end());
  }

template<typename _Tp>
  void
  test_power_norm()
  {
    std::cout.precision(std::numeric_limits<_Tp>::digits10);
    std::cout << std::showpoint << std::scientific;
    auto w = 8 + std::cout.precision();

    _Tp inf = std::numeric_limits<_Tp>::infinity();
    _Tp a = 1.5;
    _Tp b = 6.7;
    std::cout << '\n';
    std::cout << ' ' << std::setw(w) << -inf
	      << ' ' << std::setw(w) << power_norm(-inf, {a, b})
	      << '\n';
    for (_Tp p = -10.0; p <= 10.0; p += 0.015625)
      std::cout << ' ' << std::setw(w) << p
		<< ' ' << std::setw(w) << power_norm(p, {a, b})
		<< '\n';
    std::cout << ' ' << std::setw(w) << inf
	      << ' ' << std::setw(w) << power_norm(inf, {a, b})
	      << '\n';

    std::cout << '\n';
    std::cout << ' ' << std::setw(w) << -inf
	      << ' ' << std::setw(w) << power_norm(-inf, {-5.0, -3.0, -1.0, 1.0, 2.0, 3.0, 4.0, 5.0, 8.0})
	      << '\n';
    for (_Tp p = -10.0; p <= 10.0; p += 0.015625)
      std::cout << ' ' << std::setw(w) << p
		<< ' ' << std::setw(w) << power_norm(p, {-5.0, -3.0, -1.0, 1.0, 2.0, 3.0, 4.0, 5.0, 8.0})
		<< '\n';
    std::cout << ' ' << std::setw(w) << inf
	      << ' ' << std::setw(w) << power_norm(inf, {-5.0, -3.0, -1.0, 1.0, 2.0, 3.0, 4.0, 5.0, 8.0})
	      << '\n';
  }

int
main()
{
  test_power_norm<double>();
}
