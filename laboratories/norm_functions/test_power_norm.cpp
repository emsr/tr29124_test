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
template<typename Real, typename Iter>
  Real
  impl_power_norm_0(Iter begin, Iter end)
  {
    std::size_t n = 0;
    auto prod = Real{1};
    while (begin != end)
      {
	++n;
	prod *= *begin;
	++begin;
      }
    return std::pow(std::abs(prod), Real{1} / Real(n));
  }

/**
 * 
 */
template<typename Real, typename Iter>
  Real
  impl_power_norm(Real p, Iter begin, Iter end)
  {
    if (begin == end)
      return Real{0};

    using Tp = std::decay_t<decltype(*begin)>;
    if (p == Real{0})
      return impl_power_norm_0<Tp>(begin, end);
    else if (std::isinf(p))
      {
	if (p < Real{0})
	  return std::abs(*std::min_element(begin, end));
	else
	  return std::abs(*std::max_element(begin, end));
      }
    else
      {
	std::size_t n = 0;
	auto sum = Real{0};
	while (begin != end)
	  {
	    ++n;
	    sum += std::pow(std::abs(*begin), p);
	    ++begin;
	  }
	return std::pow(sum / Real(n), Real{1} / Real(p));
      }
  }

/**
 * 
 */
template<typename Real>
  Real
  power_norm(Real p, std::initializer_list<Real> val)
  {
    using Tp = std::decay_t<decltype(*val.begin())>;
    return impl_power_norm<Tp>(p, val.begin(), val.end());
  }

template<typename Tp>
  void
  test_power_norm()
  {
    std::cout.precision(std::numeric_limits<Tp>::digits10);
    std::cout << std::showpoint << std::scientific;
    auto w = 8 + std::cout.precision();

    Tp inf = std::numeric_limits<Tp>::infinity();
    Tp a = 1.5;
    Tp b = 6.7;
    std::cout << '\n';
    std::cout << ' ' << std::setw(w) << -inf
	      << ' ' << std::setw(w) << power_norm(-inf, {a, b})
	      << '\n';
    for (Tp p = -10.0; p <= 10.0; p += 0.015625)
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
    for (Tp p = -10.0; p <= 10.0; p += 0.015625)
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
