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
  __impl_power_mean_0(_Iter __begin, _Iter __end)
  {
    std::size_t __n = 0;
    auto __prod = _Real{1};
    while (__begin != __end)
      {
	++__n;
	if (*__begin < _Real{0})
	  std::__throw_domain_error("__impl_power_mean_0: "
				    "numbers must be non-negative");
	__prod *= *__begin;
	++__begin;
      }
    return std::pow(__prod, _Real{1} / _Real(__n));
  }

/**
 * 
 */
template<typename _Real, typename _Iter>
  _Real
  __impl_power_mean(_Real __p, _Iter __begin, _Iter __end)
  {
    if (__begin == __end)
      return _Real{0};

    using _Tp = std::decay_t<decltype(*__begin)>;
    if (__p == _Real{0})
      return __impl_power_mean_0<_Tp>(__begin, __end);
    else if (std::isinf(__p))
      {
	if (__p < _Real{0})
	  return *std::min_element(__begin, __end);
	else
	  return *std::max_element(__begin, __end);
      }
    else
      {
	std::size_t __n = 0;
	auto __sum = _Real{0};
	while (__begin != __end)
	  {
	    ++__n;
	    if (*__begin < _Real{0})
	      std::__throw_domain_error("__impl_power_mean: "
				    "numbers must be non-negative");
	    __sum += std::pow(*__begin, __p);
	    ++__begin;
	  }
	return std::pow(__sum / _Real(__n), _Real{1} / _Real(__p));
      }
  }

/**
 * 
 */
template<typename _Real, typename _Iter, typename _IterW>
  _Real
  __impl_power_mean_0(_Iter __vbegin, _Iter __vend, _IterW __wbegin)
  {
    auto __w = _Real{0};
    auto __prod = _Real{1};
    while (__vbegin != __vend)
      {
	if (*__vbegin < _Real{0})
	  std::__throw_domain_error("__impl_power_mean: "
				    "numbers must be non-negative");
	__w += *__wbegin;
	__prod *= std::pow(*__vbegin, *__wbegin);
	++__vbegin;
	++__wbegin;
      }
    return std::pow(__prod, _Real{1} / __w);
  }

/**
 * 
 */
template<typename _Real, typename _Iter, typename _IterW>
  _Real
  __impl_power_mean(_Real __p, _Iter __vbegin, _Iter __vend, _IterW __wbegin)
  {
    if (__vbegin == __vend)
      return _Real{0};

    using _Tp = decltype(*__vbegin * *__wbegin);
    if (__p == _Real{0})
      return __impl_power_mean_0<_Tp>(__vbegin, __vend, __wbegin);

    auto __w = _Real{0};
    auto __sum = _Real{0};
    while (__vbegin != __vend)
      {
	if (*__vbegin < _Real{0})
	  std::__throw_domain_error("__impl_power_mean: "
				    "numbers must be non-negative");
	__w += *__wbegin;
	__sum += *__wbegin * std::pow(*__vbegin, __p);
	++__vbegin;
	++__wbegin;
      }
    return std::pow(__sum / __w, _Real{1} / _Real(__p));
  }

/**
 * Compute the power mean of a sequence of non-negative numbers:
 * @f[
 *    M_p(x_1, x_2, ..., n_n) = \left( \prod_{k=1}^{n} x_k^p \right)^{1/p}
 * @f]
 */
template<typename _Real>
  _Real
  power_mean(_Real __p, std::initializer_list<_Real> __val)
  {
    using _Tp = std::decay_t<decltype(*__val.begin())>;
    return __impl_power_mean<_Tp>(__p, __val.begin(), __val.end());
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
  power_mean(_Real __p, std::initializer_list<_Real> __val,
			std::initializer_list<_Real> __weight)
  {
    using _Tp = decltype(*__val.begin() * *__weight.begin());
    if (__weight.size() != __val.size())
      std::__throw_domain_error("power_mean: number of weights"
				" must equal number of values.");
    return __impl_power_mean<_Tp>(__p, __val.begin(), __val.end(),
				  __weight.begin());
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
