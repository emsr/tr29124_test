/*
$HOME/bin/bin/g++ -std=gnu++2a -g -Wall -Wextra -Wno-psabi -I. -o test_debye test_debye.cpp -lquadmath -Lwrappers/debug -lwrap_gsl
LD_LIBRARY_PATH=wrappers/debug:$LD_LIBRARY_PATH ./test_debye > test_debye.txt

$HOME/bin/bin/g++ -std=gnu++2a -g -Wall -Wextra -Wno-psabi -I. -o test_debye test_debye.cpp -lquadmath -Lwrappers/debug -lwrap_gsl
LD_LIBRARY_PATH=wrappers/debug:$LD_LIBRARY_PATH ./test_debye > test_debye.txt
*/

#include <iostream>
#include <iomanip>
#include <limits>
#include <ext/cmath>
#include <wrap_gsl.h>

  /**
   * Return the Debye function.
   * The Debye functions are related to the incomplete Riemann zeta function:
   * @f[
   *    \zeta_x(s) = \frac{1}{\Gamma(s)}\int_{0}^{x}\frac{t^{s-1}}{e^t-1}dt
   *               = \sum_{k=1}^{\infty}\frac{P(s,kx)}{k^s}
   * @f]
   * @f[
   *    Z_x(s) = \frac{1}{\Gamma(s)}\int_{x}^{\infty}\frac{t^{s-1}}{e^t-1}dt
   *           = \sum_{k=1}^{\infty}\frac{Q(s,kx)}{k^s}
   * @f]
   * where @f$ P(a,x), Q(a,x) @f$ is the incomplete gamma function ratios.
   * The Debye function is:
   * @f[
   *    D_n(x) = \frac{n}{x^n}\int_{0}^{x}\frac{t^n}{e^t-1}dt
   *           = \Gamma(n+1)\zeta_x(n+1)
   * @f]
   * Note the infinite limit:
   * @f[
   *    D_n(\infty) = \int_{0}^{\infty}\frac{t^n}{e^t-1}dt = n!\zeta(n+1)
   * @f]
   *
   * @todo: We should return both the Debye function and it's complement.
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

	  // n!zeta(n) is the integral for x=inf, Abramowitz & Stegun 27.1.3
	  auto __sum = _Tp{0};
	  if (__n < std::__detail::_S_num_factorials<_Tp>)
	    __sum += std::__detail::__factorial<_Tp>(__n)
		   * std::__detail::__riemann_zeta<_Tp>(__n + 1);
	  else
	    return __gnu_cxx::__infinity(__x);

	  /**
	   * Compute the Debye function:
	   * @f[
	   *    D_n(x) = 1 - \sum_{k = 1}^{\infty} e^{-kx}
	   *       \frac{n}{k}\sum_{m=0}^{n}\frac{n!}{(n-m)!}frac{1}{(kx)^m}
	   * @f]
	   * Abramowitz & Stegun 27.1.2
	   */
	  const std::size_t _S_max_iter = 100;
	  auto __term = _Tp{0};
	  const auto __expmx = std::exp(-__x);
	  auto __expmkx = _Tp{1};
	  const auto __xn = std::pow(__x, _Tp(__n));
	  for(unsigned int __k = 1; __k < _S_max_iter; ++__k)
	    {
	      const auto __kx = __k * __x;
	      __expmkx *= __expmx;
	      auto __ksum = _Tp{1};
	      auto __kterm = _Tp(__n) * __ksum / __kx;  // n / (xk)^2
	      for (unsigned int __m = 1; __m <= __n; ++__m)
		__ksum += std::exchange(__kterm,
					_Tp(__n - __m) * __kterm / __kx);

	      __term -= __expmkx * __ksum * __xn / _Tp(__k);
	    }
	  __sum += __term;
	  return _Tp(__n) * __sum / __xn;
	}
      else if (std::abs(__x) < _Tp{2} * __gnu_cxx::__const_pi(__x))
	{
	  /**
	   * Compute the Debye function:
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
	  __sum *= _Tp(__n);
	  __sum += _Tp{1} - _Tp(__n) * __x / _Tp(2 * (__n + 1));
	  return __sum;
	}
      else
	return _Tp{0}; /// @todo Find Debye for x < -2pi!
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
	      << ' ' << std::setw(width) << "Debye_6(x)"
	      << '\n';
    for (int i = -50; i <= +200; ++i)
      {
	auto x = _Tp{0.1L} * i;
	std::cout << ' ' << std::setw(width) << x;
	for (int n = 1; n <= 6; ++n)
	  std::cout << ' ' << std::setw(width) << __debye(n, x);
	std::cout << '\n';
      }
  }

template<typename _Tp>
  void
  test_debye_gsl(_Tp __proto = _Tp{})
  {
    std::cout.precision(__gnu_cxx::__max_digits10(__proto));
    std::cout << std::showpoint << std::scientific;
    auto width = 8 + std::cout.precision();

    std::cout << '\n' << '\n';
    std::cout << ' ' << std::setw(width) << "x"
	      << ' ' << std::setw(width) << "Debye_1(x) GSL"
	      << ' ' << std::setw(width) << "Debye_2(x) GSL"
	      << ' ' << std::setw(width) << "Debye_3(x) GSL"
	      << ' ' << std::setw(width) << "Debye_4(x) GSL"
	      << ' ' << std::setw(width) << "Debye_5(x) GSL"
	      << ' ' << std::setw(width) << "Debye_6(x) GSL"
	      << '\n';
    for (int i = 0; i <= +200; ++i)
      {
	auto x = _Tp{0.1L} * i;
	std::cout << ' ' << std::setw(width) << x;
	for (int n = 1; n <= 6; ++n)
	  std::cout << ' ' << std::setw(width) << gsl::debye(n, x);
	std::cout << '\n';
      }
  }

int
main()
{
  test_debye(1.0);

  test_debye_gsl(1.0);

  return 0;
}

