/**
 *
 */

#include <iostream>
#include <iomanip>
#include <limits>
#include <cmath>

#include <emsr/numeric_limits.h>
#include <emsr/sf_gamma.h> // factorial
#include <emsr/sf_zeta.h>

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
  template<typename Tp>
    Tp
    debye(unsigned int n, Tp x)
    {
      if (std::isnan(x))
	return std::numeric_limits<Tp>::quiet_NaN();
      else if (n < 1)
	throw std::domain_error("debye: Degree n must be positive.");
      else if (x >= Tp{3})
	{
	  // For values up to 4.80 the list of zeta functions
	  // and the sum up to k < K are huge enough to gain
	  // numeric stability in the sum

	  // n!zeta(n) is the integral for x=inf, Abramowitz & Stegun 27.1.3
	  auto sum = Tp{0};
	  if (n < emsr::detail::s_num_factorials<Tp>)
	    sum += emsr::detail::factorial<Tp>(n)
		   * emsr::detail::riemann_zeta<Tp>(n + 1);
	  else
	    return emsr::infinity(x);

	  /**
	   * Compute the Debye function:
	   * @f[
	   *    D_n(x) = 1 - \sum_{k = 1}^{\infty} e^{-kx}
	   *       \frac{n}{k}\sum_{m=0}^{n}\frac{n!}{(n-m)!}frac{1}{(kx)^m}
	   * @f]
	   * Abramowitz & Stegun 27.1.2
	   */
	  const std::size_t s_max_iter = 100;
	  auto term = Tp{0};
	  const auto expmx = std::exp(-x);
	  auto expmkx = Tp{1};
	  const auto xn = std::pow(x, Tp(n));
	  for(unsigned int k = 1; k < s_max_iter; ++k)
	    {
	      const auto kx = k * x;
	      expmkx *= expmx;
	      auto ksum = Tp{1};
	      auto kterm = Tp(n) * ksum / kx;  // n / (xk)^2
	      for (unsigned int m = 1; m <= n; ++m)
		ksum += std::exchange(kterm,
					Tp(n - m) * kterm / kx);

	      term -= expmkx * ksum * xn / Tp(k);
	    }
	  sum += term;
	  return Tp(n) * sum / xn;
	}
      else if (std::abs(x) < Tp{2} * emsr::pi_v<Tp>)
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
	  const auto s_eps = emsr::epsilon(x);
	  const std::size_t s_max_iter = 200;
	  const auto s_1_2pi = emsr::inv_tau_v<Tp>;
	  const auto x2pi = x * s_1_2pi;
	  const auto x2pi2 = x2pi * x2pi;
	  auto x2pi2k = x2pi2;
	  auto sum = Tp{0};
	  for(unsigned int k = 1; k < s_max_iter; ++k)
	    {
	      const auto term = Tp{2}
				* emsr::detail::riemann_zeta<Tp>(2 * k)
				* x2pi2k / Tp(2 * k + n);
	      sum += term;
	      if (std::abs(term) < s_eps * std::abs(sum))
        	break;
	      x2pi2k *= -x2pi2;
	    }
	  sum *= Tp(n);
	  sum += Tp{1} - Tp(n) * x / Tp(2 * (n + 1));
	  return sum;
	}
      else
	return Tp{0}; /// @todo Find Debye for x < -2pi!
    }

template<typename Tp>
  void
  test_debye(Tp proto = Tp{})
  {
    std::cout.precision(emsr::max_digits10(proto));
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
	auto x = Tp{0.1L} * i;
	std::cout << ' ' << std::setw(width) << x;
	for (int n = 1; n <= 6; ++n)
	  std::cout << ' ' << std::setw(width) << debye(n, x);
	std::cout << '\n';
      }
  }

template<typename Tp>
  void
  test_debye_gsl(Tp proto = Tp{})
  {
    std::cout.precision(emsr::max_digits10(proto));
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
	auto x = Tp{0.1L} * i;
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

