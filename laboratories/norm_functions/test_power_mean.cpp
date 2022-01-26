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
  impl_power_mean_0(_Iter begin, _Iter end)
  {
    std::size_t n = 0;
    auto prod = _Real{1};
    while (begin != end)
      {
	++n;
	if (*begin < _Real{0})
	  throw std::domain_error("impl_power_mean_0: numbers must be non-negative");
	prod *= *begin;
	++begin;
      }
    return std::pow(prod, _Real{1} / _Real(n));
  }

/**
 * 
 */
template<typename _Real, typename _Iter>
  _Real
  impl_power_mean(_Real p, _Iter begin, _Iter end)
  {
    if (begin == end)
      return _Real{0};

    using _Tp = std::decay_t<decltype(*begin)>;
    if (p == _Real{0})
      return impl_power_mean_0<_Tp>(begin, end);
    else if (std::isinf(p))
      {
	if (p < _Real{0})
	  return *std::min_element(begin, end);
	else
	  return *std::max_element(begin, end);
      }
    else
      {
	std::size_t n = 0;
	auto sum = _Real{0};
	while (begin != end)
	  {
	    ++n;
	    if (*begin < _Real{0})
	      throw std::domain_error("impl_power_mean: numbers must be non-negative");
	    sum += std::pow(*begin, p);
	    ++begin;
	  }
	return std::pow(sum / _Real(n), _Real{1} / _Real(p));
      }
  }

/**
 * 
 */
template<typename _Real, typename _Iter, typename _IterW>
  _Real
  impl_power_mean_0(_Iter vbegin, _Iter vend, _IterW wbegin)
  {
    auto w = _Real{0};
    auto prod = _Real{1};
    while (vbegin != vend)
      {
	if (*vbegin < _Real{0})
	  throw std::domain_error("impl_power_mean: numbers must be non-negative");
	w += *wbegin;
	prod *= std::pow(*vbegin, *wbegin);
	++vbegin;
	++wbegin;
      }
    return std::pow(prod, _Real{1} / w);
  }

/**
 * 
 */
template<typename _Real, typename _Iter, typename _IterW>
  _Real
  impl_power_mean(_Real p, _Iter vbegin, _Iter vend, _IterW wbegin)
  {
    if (vbegin == vend)
      return _Real{0};

    using _Tp = decltype(*vbegin * *wbegin);
    if (p == _Real{0})
      return impl_power_mean_0<_Tp>(vbegin, vend, wbegin);

    auto w = _Real{0};
    auto sum = _Real{0};
    while (vbegin != vend)
      {
	if (*vbegin < _Real{0})
	  throw std::domain_error("impl_power_mean: numbers must be non-negative");
	w += *wbegin;
	sum += *wbegin * std::pow(*vbegin, p);
	++vbegin;
	++wbegin;
      }
    return std::pow(sum / w, _Real{1} / _Real(p));
  }

/**
 * Compute the power mean of a sequence of non-negative numbers:
 * @f[
 *    M_p(x_1, x_2, ..., n_n) = \left( \prod_{k=1}^{n} x_k^p \right)^{1/p}
 * @f]
 */
template<typename _Real>
  _Real
  power_mean(_Real p, std::initializer_list<_Real> val)
  {
    using _Tp = std::decay_t<decltype(*val.begin())>;
    return impl_power_mean<_Tp>(p, val.begin(), val.end());
  }

/**
 * Compute the power mean of a sequence of non-negative numbers
 * with non-negative weights:
 * @f[
 *    M_p((x_1, w_1), (x_2, w_2), ..., (n_n, w_n)
 *        = \left( \prod_{k=1}^{n} x_k^p \right)^{1/p}
 * @f]
 */
template<typename _Real>
  _Real
  power_mean(_Real p, std::initializer_list<_Real> val,
			std::initializer_list<_Real> weight)
  {
    using _Tp = decltype(*val.begin() * *weight.begin());
    if (weight.size() != val.size())
      throw std::domain_error("power_mean: number of weights must equal number of values.");
    return impl_power_mean<_Tp>(p, val.begin(), val.end(),
				  weight.begin());
  }

template<typename _Tp>
  void
  test_power_mean()
  {
    std::cout.precision(std::numeric_limits<_Tp>::digits10);
    std::cout << std::showpoint << std::scientific;
    auto w = 8 + std::cout.precision();

    _Tp inf = std::numeric_limits<_Tp>::infinity();
    _Tp a = 1.5;
    _Tp b = 6.7;
    std::cout << ' ' << std::setw(w) << -inf
	      << ' ' << std::setw(w) << power_mean(-inf, {a, b})
	      << '\n';
    for (_Tp p = -10.0; p <= 10.0; p += 0.015625)
      std::cout << ' ' << std::setw(w) << p
		<< ' ' << std::setw(w) << power_mean(p, {a, b})
		<< '\n';
    std::cout << ' ' << std::setw(w) << inf
	      << ' ' << std::setw(w) << power_mean(inf, {a, b})
	      << '\n';
  }

int
main()
{
  test_power_mean<double>();
}
