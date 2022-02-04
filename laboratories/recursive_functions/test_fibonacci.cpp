/**
 *
 */

#include <cmath>
#include <limits>
#include <iostream>
#include <iomanip>

#include <emsr/math_constants.h>
#include <emsr/math_util.h>
#include <emsr/sf_trig.h>

/**
 * @todo Test these relations between the Lucas and Fibonacci numbers
 *     and the Chebyshev polynomials of the first and second kinds respectively.
 * @f[
 *   2i^{-n} T_n(i/2) = L_n
 * @f]
 * and
 * @f[
 *   i^n U_n(i/2) = F_{n+1}
 * @f]
 */

/**
 * Return the Fibonacci number for unsigned integers by recursion:
 * @f[
 *    F_n = F_{n-1} + F_{n-2}, F_0 = 0, F_1 = 1
 * @f]
 */
template<typename _UIntTp>
  _UIntTp
  fibonacci_recur(_UIntTp n)
  {
    _UIntTp Fnm2 = 0;
    if (n == 0)
      return Fnm2;
    _UIntTp Fnm1 = 1;
    if (n == 1)
      return Fnm1;
    _UIntTp Fn = Fnm1 + Fnm2;
    for (_UIntTp k = 3; k <= n; ++k)
      {
	Fnm2 = Fnm1;
	Fnm1 = Fn;
	if (__builtin_add_overflow(Fnm1, Fnm2, &Fn))
	    throw std::runtime_error("fibonacci: integer overflow");
      }
    return Fn;
  }

/**
 * Return the Fibonacci number for integral or real degree.
 * The real-valued Fibonacci function matches the Fibonacci numbers
 * for integral arguments:
 * @f[
 *    F_\nu = \frac{1}{\sqrt{5}}\left[\phi^\nu - \cos(\nu\pi)\phi^{-\nu}\right]
 * @f]
 * where @f$ \phi = (1+\sqrt{5})/2 @g$ is the Golden ratio.
 * Note that the Fibonacci numbers for negative integers are given
 * by @f$ F_{-n} = (-1)^{n+1} F_n @f$.
 */
template<typename Tp>
  Tp
  fibonacci(Tp nu)
  {
    if constexpr (std::is_integral_v<Tp>)
      {
	if (std::is_unsigned_v<Tp>)
	  return fibonacci_recur(nu);
	else
	  {
	    if (nu < 0)
	      return emsr::parity<Tp>(-nu + 1)
		   * fibonacci_recur(-nu);
	    else
	      return fibonacci_recur(nu);
	  }
      }
    else if constexpr (std::is_floating_point_v<Tp>)
      {
	const auto _S_phi = emsr::phi_v<Tp>;
	const auto _S_sqrt5 = emsr::sqrt5_v<Tp>;
	const auto phinu = std::pow(_S_phi, nu);
	return (phinu - emsr::cos_pi(nu) / phinu) / _S_sqrt5;
      }
  }

/**
 * Return the Fibonacci polynomial of order @c n and argument @c x by recursion:
 * @f[
 *    F_n(x) = xF_{n-1}(x) + F_{n-2}(x), F_0 = 2, F_1 = 1
 * @f]
 */
template<typename _UIntTp, typename _RealTp>
  _RealTp
  fibonacci(_UIntTp n, _RealTp x)
  {
    //const auto _S_log_phi
    //  = _RealTp{4.812118250596034474977589134243684231358e-1L};
    auto Fnm2 = _RealTp{0};
    if (n == 0)
      return Fnm2;
    auto Fnm1 = _RealTp{1};
    if (n == 1)
      return Fnm1;
    auto Fn = x * Fnm1 + Fnm2;
    for (_UIntTp k = 3; k <= n; ++k)
      {
	Fnm2 = Fnm1;
	Fnm1 = Fn;
	Fn = x * Fnm1 + Fnm2;
      }
    return Fn;
  }

/**
 * Return the Lucas number for unsigned integers by recursion:
 * @f[
 *    L_n = L_{n-1} + L_{n-2}, L_0 = 2, L_1 = 1
 * @f]
 */
template<typename _UIntTp>
  _UIntTp
  lucas_recur(_UIntTp n)
  {
    _UIntTp Lnm2 = 2;
    if (n == 0)
      return Lnm2;
    _UIntTp Lnm1 = 1;
    if (n == 1)
      return Lnm1;
    _UIntTp Ln = Lnm1 + Lnm2;
    for (_UIntTp k = 3; k <= n; ++k)
      {
	Lnm2 = Lnm1;
	Lnm1 = Ln;
	if (__builtin_add_overflow(Lnm1, Lnm2, &Ln))
	    throw std::runtime_error("lucas: integer overflow");
      }
    return Ln;
  }

/**
 * Return the Lucas number for integral or real degree.
 * The real-valued Lucas function matches the Lucas numbers
 * for integral arguments:
 * @f[
 *    L_\nu = \phi^\nu + \cos(\nu\pi) \phi^{-\nu}
 * @f]
 * where @f$ \phi = (1+\sqrt{5})/2 @g$ is the Golden ratio.
 * Note that the Lucas numbers for negative integers are given
 * by @f$ L_{-n} = (-1)^n L_n @f$.
 */
template<typename Tp>
  Tp
  lucas(Tp nu)
  {
    if constexpr (std::is_integral_v<Tp>)
      {
	if (std::is_unsigned_v<Tp>)
	  return lucas_recur(nu);
	else
	  {
	    if (nu < 0)
	      return emsr::parity<Tp>(nu) * lucas_recur(-nu);
	    else
	      return lucas_recur(nu);
	  }
      }
    else if constexpr (std::is_floating_point_v<Tp>)
      {
	const auto _S_phi = emsr::phi_v<Tp>;
	const auto phinu = std::pow(_S_phi, nu);
	return phinu + emsr::cos_pi(nu) / phinu;
      }
  }

/**
 * Return the Lucas polynomial of order @c n and argument @c x by recursion:
 * @f[
 *    L_n(x) = xL_{n-1}(x) + L_{n-2}(x), L_0 = 2, L_1 = x
 * @f]
 */
template<typename _UIntTp, typename _RealTp>
  _RealTp
  lucas(_UIntTp n, _RealTp x)
  {
    //const auto _S_log_phi
    //  = _RealTp{4.812118250596034474977589134243684231358e-1L};
    auto Lnm2 = _RealTp{2};
    if (n == 0)
      return Lnm2;
    auto Lnm1 = x;
    if (n == 1)
      return Lnm1;
    auto Ln = x * Lnm1 + Lnm2;
    for (_UIntTp k = 3; k <= n; ++k)
      {
	Lnm2 = Lnm1;
	Lnm1 = Ln;
	Ln = x * Lnm1 + Lnm2;
      }
    return Ln;
  }

template<typename Tp>
  void
  test_fibonacci()
  {
    std::cout.precision(std::numeric_limits<Tp>::digits10);
    std::cout << std::showpoint << std::scientific;
    auto width = 8 + std::cout.precision();

    const auto max_number = 50ll;
    const auto delnu = Tp{1} / Tp{50};
    const auto max_order = 50ll;

    std::cout << "\n\n Fibonacci numbers\n";
    for (auto n = -max_number; n <= max_number; ++n)
      {
	auto _F_n = fibonacci(n);
	std::cout << ' ' << std::setw(4) << n
		  << ' ' << std::setw(width) << _F_n
		  << '\n';
      }

    std::cout << "\n\n Fibonacci function\n";
    for (int n = -500; n <= 500; ++n)
      {
	auto nu = n * delnu;
	auto F_nu = fibonacci(nu);
	std::cout << ' ' << std::setw(width) << nu
		  << ' ' << std::setw(width) << F_nu
		  << '\n';
      }

    std::cout << "\n\n Fibonacci polynomials\n";
    for (auto n = 0; n <= max_order; ++n)
      {
	std::cout << '\n' << '\n' << ' ' << std::setw(4) << n << '\n';
	const auto del = Tp{1} / Tp{10};
	for (int i = -50; i <= 50; ++i)
	  {
	    auto x = del * i;
	    auto _F_n = fibonacci(n, x);
	    std::cout << ' ' << std::setw(width) << x
		      << ' ' << std::setw(width) << _F_n
		      << '\n';
	  }
      }

    std::cout << "\n\n Lucas numbers\n";
    for (auto n = -max_number; n <= max_number; ++n)
      {
	auto _L_n = lucas(n);
	std::cout << ' ' << std::setw(4) << n
		  << ' ' << std::setw(width) << _L_n
		  << '\n';
      }

    std::cout << "\n\n Lucas function\n";
    for (int n = -500; n <= 500; ++n)
      {
	auto nu = n * delnu;
	auto L_nu = lucas(nu);
	std::cout << ' ' << std::setw(width) << nu
		  << ' ' << std::setw(width) << L_nu
		  << '\n';
      }

    std::cout << "\n\n Lucas polynomials\n";
    for (auto n = 0; n <= max_order; ++n)
      {
	std::cout << '\n' << '\n' << ' ' << std::setw(4) << n << '\n';
	const auto del = Tp{1} / Tp{10};
	for (int i = -50; i <= 50; ++i)
	  {
	    auto x = del * i;
	    auto _L_n = lucas(n, x);
	    std::cout << ' ' << std::setw(width) << x
		      << ' ' << std::setw(width) << _L_n
		      << '\n';
	  }
      }
  }

int
main()
{
  test_fibonacci<double>();
}
