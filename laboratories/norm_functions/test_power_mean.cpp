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
  impl_power_mean_0(Iter begin, Iter end)
  {
    std::size_t n = 0;
    auto prod = Real{1};
    while (begin != end)
      {
	++n;
	if (*begin < Real{0})
	  throw std::domain_error("impl_power_mean_0: numbers must be non-negative");
	prod *= *begin;
	++begin;
      }
    return std::pow(prod, Real{1} / Real(n));
  }

/**
 * 
 */
template<typename Real, typename Iter>
  Real
  impl_power_mean(Real p, Iter begin, Iter end)
  {
    if (begin == end)
      return Real{0};

    using Tp = std::decay_t<decltype(*begin)>;
    if (p == Real{0})
      return impl_power_mean_0<Tp>(begin, end);
    else if (std::isinf(p))
      {
	if (p < Real{0})
	  return *std::min_element(begin, end);
	else
	  return *std::max_element(begin, end);
      }
    else
      {
	std::size_t n = 0;
	auto sum = Real{0};
	while (begin != end)
	  {
	    ++n;
	    if (*begin < Real{0})
	      throw std::domain_error("impl_power_mean: numbers must be non-negative");
	    sum += std::pow(*begin, p);
	    ++begin;
	  }
	return std::pow(sum / Real(n), Real{1} / Real(p));
      }
  }

/**
 * 
 */
template<typename Real, typename Iter, typename _IterW>
  Real
  impl_power_mean_0(Iter vbegin, Iter vend, _IterW wbegin)
  {
    auto w = Real{0};
    auto prod = Real{1};
    while (vbegin != vend)
      {
	if (*vbegin < Real{0})
	  throw std::domain_error("impl_power_mean: numbers must be non-negative");
	w += *wbegin;
	prod *= std::pow(*vbegin, *wbegin);
	++vbegin;
	++wbegin;
      }
    return std::pow(prod, Real{1} / w);
  }

/**
 * 
 */
template<typename Real, typename Iter, typename _IterW>
  Real
  impl_power_mean(Real p, Iter vbegin, Iter vend, _IterW wbegin)
  {
    if (vbegin == vend)
      return Real{0};

    using Tp = decltype(*vbegin * *wbegin);
    if (p == Real{0})
      return impl_power_mean_0<Tp>(vbegin, vend, wbegin);

    auto w = Real{0};
    auto sum = Real{0};
    while (vbegin != vend)
      {
	if (*vbegin < Real{0})
	  throw std::domain_error("impl_power_mean: numbers must be non-negative");
	w += *wbegin;
	sum += *wbegin * std::pow(*vbegin, p);
	++vbegin;
	++wbegin;
      }
    return std::pow(sum / w, Real{1} / Real(p));
  }

/**
 * Compute the power mean of a sequence of non-negative numbers:
 * @f[
 *    M_p(x_1, x_2, ..., n_n) = \left( \prod_{k=1}^{n} x_k^p \right)^{1/p}
 * @f]
 */
template<typename Real>
  Real
  power_mean(Real p, std::initializer_list<Real> val)
  {
    using Tp = std::decay_t<decltype(*val.begin())>;
    return impl_power_mean<Tp>(p, val.begin(), val.end());
  }

/**
 * Compute the power mean of a sequence of non-negative numbers
 * with non-negative weights:
 * @f[
 *    M_p((x_1, w_1), (x_2, w_2), ..., (n_n, w_n)
 *        = \left( \prod_{k=1}^{n} x_k^p \right)^{1/p}
 * @f]
 */
template<typename Real>
  Real
  power_mean(Real p, std::initializer_list<Real> val,
			std::initializer_list<Real> weight)
  {
    using Tp = decltype(*val.begin() * *weight.begin());
    if (weight.size() != val.size())
      throw std::domain_error("power_mean: number of weights must equal number of values.");
    return impl_power_mean<Tp>(p, val.begin(), val.end(),
				  weight.begin());
  }

template<typename Tp>
  void
  test_power_mean()
  {
    std::cout.precision(std::numeric_limits<Tp>::digits10);
    std::cout << std::showpoint << std::scientific;
    auto w = 8 + std::cout.precision();

    Tp inf = std::numeric_limits<Tp>::infinity();
    Tp a = 1.5;
    Tp b = 6.7;
    std::cout << ' ' << std::setw(w) << -inf
	      << ' ' << std::setw(w) << power_mean(-inf, {a, b})
	      << '\n';
    for (Tp p = -10.0; p <= 10.0; p += 0.015625)
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
