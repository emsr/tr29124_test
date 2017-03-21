/*
$HOME/bin_tr29124/bin/g++ -std=gnu++17 -g -Wall -Wextra -Wno-psabi -I. -o test_debye test_debye.cpp -lquadmath -Lwrappers/debug -lwrap_boost
LD_LIBRARY_PATH=wrappers/debug:$LD_LIBRARY_PATH ./test_debye > test_debye.txt

PATH=wrappers/debug:$PATH $HOME/bin/bin/g++ -std=gnu++17 -g -Wall -Wextra -Wno-psabi -I. -o test_debye test_debye.cpp -lquadmath -Lwrappers/debug -lwrap_boost
./test_debye > test_debye.txt
*/

#include <iostream>
#include <iomanip>
#include <limits>
#include <ext/cmath>

  /**
   * Return the Debye integral or the incomplete Riemann zeta function:
   * @todo: We should return both the integral and it's complement.
   * @f[
   *    \zeta_x(s) = \frac{1}{\Gamma(s)}\int_{0}^{x}\frac{t^{s-1}}{e^t-1}dt
   *          = \sum{k=1}{\infty}k^{-s}P(s,kx)
   * @f]
   * @f[
   *    \Zeta_x(s) = \frac{1}{\Gamma(s)}\int_{x}^{\infty}\frac{t^{s-1}}{e^t-1}dt
   *          = \sum{k=1}{\infty}k^{-s}Q(s,kx)
   * @f]
   * where @f$ P(a,x), Q(a,x) @f$ is the incomplete gamma function ratios.
   * The Debye integrals are:
   * @f[
   *    D_n(x) = \frac{n}{x^n}\int_{0}^{x}\frac{t^n}{e^t-1}dt
   *           = \Gamma(n+1)[\zeta(n+1)-\zeta_x(n+1)]
   * @f]
   * and
   * @f[
   *    \int_{0}^{x}\frac{t^n}{e^t-1}dt = \Gamma(n+1)\zeta_x(n+1)
   * @f]
   */
  template<typename _Tp>
    _Tp
    __debye(unsigned int __n, _Tp __x)
    {
      if (std::isnan(__x))
	return std::numeric_limits<_Tp>::quiet_NaN();
      else if (__n < 1)
	std::__throw_domain_error("__debye: Degree n must be positive.");
      else if (__x >= _Tp{3})
	{
	  // For values up to 4.80 the list of zeta functions
	  // and the sum up to k < K are huge enough to gain
	  // numeric stability in the sum

	  const auto _S_eps = __gnu_cxx::__epsilon(__x);
	  const std::size_t _S_max_iter = 100;

	  // n!zeta(n) is the integral for x=inf, Abramowitz & Stegun 27.1.3
	  auto __sum = _Tp{0};
	  if (std::__detail::_S_num_factorials<_Tp>)
	    __sum += std::__detail::__factorial<_Tp>(__n)
		   * std::__detail::__riemann_zeta<_Tp>(__n + 1);
	  else
	    return __gnu_cxx::__infinity(__x);

	  /**
	   * Compute the Debye integral:
	   * @f[
	   *    D_n(x) = 1 - \sum_{k = 1}^{\infty} e^{-kx}
	   *       \frac{n}{k}\sum_{m=0}^{n}\frac{n!}{(n-m)!}frac{1}{(kx)^m}
	   * @f]
	   * Abramowitz & Stegun 27.1.2
	   */
	  auto __term = _Tp{0};
	  for(unsigned int __k = 1; __k < _S_max_iter; ++__k)
	    {
	      const auto __xk = __x * __k;
	      auto __ksum = _Tp{1} / __xk;
	      auto __kterm = _Tp(__n) * __ksum / __xk;  // n / (xk)^2
	      for (unsigned int __s = 1; __s <= __n; ++__s)
		__ksum += std::exchange(__kterm,
					_Tp(__n - __s) * __kterm / __xk);

	      __term -= std::exp(-__xk) * __ksum * std::pow(__x, _Tp(__n + 1));
	    }
	  __sum += __term;
	  if (std::abs(__term) < _S_eps * std::abs(__sum))
	    return __sum;
	  return __sum;
	}
      else
	{
	  /**
	   * Compute the Debye integral:
	   * @f[
	   *    D_n(x) = 1 - \frac{n x}{2(n+1)}
	   *       + n \sum_{k = 1}^{\infty} \frac{B_{2k} x^{2k}}{(2k + n)(2k)!}
	   * @f]
           * for @f$ |x| < 2\pi @f$.
	   * Abramowitz-Stegun 27.1.1
	   */
	  const auto _S_eps = __gnu_cxx::__epsilon(__x);
	  const std::size_t _S_max_iter = 200;
	  const auto _S_1_2pi = __gnu_cxx::__const_one_div_2_pi(__x);
	  const auto __x2pi = __x * _S_1_2pi;
	  const auto __x2pi2 = __x2pi * __x2pi;
	  auto __x2pi2k = __x2pi2;
	  auto __sum = _Tp{0};
	  for(unsigned int __k = 1; __k < _S_max_iter; ++__k)
	    {
	      const auto __term = _Tp{2}
				* std::__detail::__riemann_zeta<_Tp>(2 * __k)
				* __x2pi2k / _Tp(2 * __k + __n);
	      __sum += __term;
	      if (std::abs(__term) < _S_eps * std::abs(__sum))
        	break;
	      __x2pi2k *= -__x2pi2;
	    }
	  __sum += _Tp{1} / _Tp(__n) - __x / _Tp(2 * (1 + __n));
	  return __sum * std::pow(__x, _Tp(__n));
	}
    }

template<typename _Tp>
  void
  test_debye(_Tp __proto = _Tp{})
  {
    std::cout.precision(__gnu_cxx::__max_digits10(__proto));
    std::cout << std::showpoint << std::scientific;
    auto width = 8 + std::cout.precision();

    std::cout << '\n' << '\n';
    std::cout << ' ' << std::setw(width) << "x"
	      << ' ' << std::setw(width) << "Debye_1(x)"
	      << ' ' << std::setw(width) << "Debye_2(x)"
	      << ' ' << std::setw(width) << "Debye_3(x)"
	      << ' ' << std::setw(width) << "Debye_4(x)"
	      << ' ' << std::setw(width) << "Debye_5(x)"
	      << '\n';
    for (int i = -50; i <= +200; ++i)
      {
	auto x = _Tp{0.1L} * i;
	std::cout << ' ' << std::setw(width) << x;
	for (int n = 1; n <= 5; ++n)
	  std::cout << ' ' << std::setw(width) << __debye(n, x);
	std::cout << '\n';
      }
  }

int
main()
{
  test_debye(1.0);

  return 0;
}

