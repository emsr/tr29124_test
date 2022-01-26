/**
 *
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
  impl_power_norm_0(_Iter begin, _Iter end)
  {
    std::size_t n = 0;
    auto prod = _Real{1};
    while (begin != end)
      {
	++n;
	prod *= *begin;
	++begin;
      }
    return std::pow(std::abs(prod), _Real{1} / _Real(n));
  }

/**
 * 
 */
template<typename _Real, typename _Iter>
  _Real
  impl_power_norm(_Real p, _Iter begin, _Iter end)
  {
    if (begin == end)
      return _Real{0};

    using _Tp = std::decay_t<decltype(*begin)>;
    if (p == _Real{0})
      return impl_power_norm_0<_Tp>(begin, end);
    else if (std::isinf(p))
      {
	if (p < _Real{0})
	  return std::abs(*std::min_element(begin, end));
	else
	  return std::abs(*std::max_element(begin, end));
      }
    else
      {
	std::size_t n = 0;
	auto sum = _Real{0};
	while (begin != end)
	  {
	    ++n;
	    sum += std::pow(std::abs(*begin), p);
	    ++begin;
	  }
	return std::pow(sum / _Real(n), _Real{1} / _Real(p));
      }
  }

/**
 * 
 */
template<typename _Real>
  _Real
  power_norm(_Real p, std::initializer_list<_Real> val)
  {
    using _Tp = std::decay_t<decltype(*val.begin())>;
    return impl_power_norm<_Tp>(p, val.begin(), val.end());
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
