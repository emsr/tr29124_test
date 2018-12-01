/*
$HOME/bin/bin/g++ -std=gnu++2a -g -Wall -Wextra -Wno-psabi -I. -o test_fibonacci test_fibonacci.cpp -Lwrappers/debug -lquadmath -Lwrappers/debug -lwrap_gsl -lwrap_burkhardt -lgfortran
LD_LIBRARY_PATH=wrappers/debug:$HOME/bin/lib64:$LD_LIBRARY_PATH ./test_fibonacci > test_fibonacci.txt

$HOME/bin/bin/g++ -std=gnu++2a -g -Wall -Wextra -I. -o test_fibonacci test_fibonacci.cpp -lquadmath -Lwrappers/debug -lwrap_gsl -lwrap_burkhardt
./test_fibonacci > test_fibonacci.txt
*/

#include <cmath>
#include <limits>
#include <iostream>
#include <iomanip>

/**
 * @todo Test these relations between the Lucas and Fibonacci numbers
 *     and the Chebyshev polynomials of the first and second kinds respectively.
 * @f[
 *   2i^{-n} T_n(i/2) = L_n
 * @f]
 * and
 * @f[
 *   i^n U(i/2) = F_{n+1}
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
  __fibonacci_recur(_UIntTp __n)
  {
    _UIntTp __Fnm2 = 0;
    if (__n == 0)
      return __Fnm2;
    _UIntTp __Fnm1 = 1;
    if (__n == 1)
      return __Fnm1;
    _UIntTp __Fn = __Fnm1 + __Fnm2;
    for (_UIntTp __k = 3; __k <= __n; ++__k)
      {
	__Fnm2 = __Fnm1;
	__Fnm1 = __Fn;
	if (__builtin_add_overflow(__Fnm1, __Fnm2, &__Fn))
	    std::__throw_runtime_error(__N("__fibonacci: "
					   "integer overflow"));
      }
    return __Fn;
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
template<typename _Tp>
  _Tp
  __fibonacci(_Tp __nu)
  {
    if constexpr (std::is_integral_v<_Tp>)
      {
	if (std::is_unsigned_v<_Tp>)
	  return __fibonacci_recur(__nu);
	else
	  {
	    if (__nu < 0)
	      return __gnu_cxx::__parity<_Tp>(-__nu + 1)
		   * __fibonacci_recur(-__nu);
	    else
	      return __fibonacci_recur(__nu);
	  }
      }
    else if constexpr (std::is_floating_point_v<_Tp>)
      {
	const auto _S_phi = __gnu_cxx::__const_phi(__nu);
	const auto _S_sqrt5 = __gnu_cxx::__const_root_5(__nu);
	const auto __phinu = std::pow(_S_phi, __nu);
	return (__phinu - __gnu_cxx::cos_pi(__nu) / __phinu) / _S_sqrt5;
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
  __fibonacci(_UIntTp __n, _RealTp __x)
  {
    //const auto _S_log_phi
    //  = _RealTp{4.812118250596034474977589134243684231358e-1L};
    auto __Fnm2 = _RealTp{0};
    if (__n == 0)
      return __Fnm2;
    auto __Fnm1 = _RealTp{1};
    if (__n == 1)
      return __Fnm1;
    auto __Fn = __x * __Fnm1 + __Fnm2;
    for (_UIntTp __k = 3; __k <= __n; ++__k)
      {
	__Fnm2 = __Fnm1;
	__Fnm1 = __Fn;
	__Fn = __x * __Fnm1 + __Fnm2;
      }
    return __Fn;
  }

/**
 * Return the Lucas number for unsigned integers by recursion:
 * @f[
 *    L_n = L_{n-1} + L_{n-2}, L_0 = 2, L_1 = 1
 * @f]
 */
template<typename _UIntTp>
  _UIntTp
  __lucas_recur(_UIntTp __n)
  {
    _UIntTp __Lnm2 = 2;
    if (__n == 0)
      return __Lnm2;
    _UIntTp __Lnm1 = 1;
    if (__n == 1)
      return __Lnm1;
    _UIntTp __Ln = __Lnm1 + __Lnm2;
    for (_UIntTp __k = 3; __k <= __n; ++__k)
      {
	__Lnm2 = __Lnm1;
	__Lnm1 = __Ln;
	if (__builtin_add_overflow(__Lnm1, __Lnm2, &__Ln))
	    std::__throw_runtime_error(__N("__lucas: "
					   "integer overflow"));
      }
    return __Ln;
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
template<typename _Tp>
  _Tp
  __lucas(_Tp __nu)
  {
    if constexpr (std::is_integral_v<_Tp>)
      {
	if (std::is_unsigned_v<_Tp>)
	  return __lucas_recur(__nu);
	else
	  {
	    if (__nu < 0)
	      return __gnu_cxx::__parity<_Tp>(__nu) * __lucas_recur(-__nu);
	    else
	      return __lucas_recur(__nu);
	  }
      }
    else if constexpr (std::is_floating_point_v<_Tp>)
      {
	const auto _S_phi = __gnu_cxx::__const_phi(__nu);
	const auto __phinu = std::pow(_S_phi, __nu);
	return __phinu + __gnu_cxx::cos_pi(__nu) / __phinu;
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
  __lucas(_UIntTp __n, _RealTp __x)
  {
    //const auto _S_log_phi
    //  = _RealTp{4.812118250596034474977589134243684231358e-1L};
    auto __Lnm2 = _RealTp{2};
    if (__n == 0)
      return __Lnm2;
    auto __Lnm1 = __x;
    if (__n == 1)
      return __Lnm1;
    auto __Ln = __x * __Lnm1 + __Lnm2;
    for (_UIntTp __k = 3; __k <= __n; ++__k)
      {
	__Lnm2 = __Lnm1;
	__Lnm1 = __Ln;
	__Ln = __x * __Lnm1 + __Lnm2;
      }
    return __Ln;
  }

template<typename _Tp>
  void
  test_fibonacci()
  {
    std::cout.precision(__gnu_cxx::__digits10<_Tp>());
    std::cout << std::showpoint << std::scientific;
    auto width = 8 + std::cout.precision();
    const auto max_number = 50ll;
    const auto delnu = _Tp{1} / _Tp{50};
    const auto max_order = 50ll;

    std::cout << "\n\n Fibonacci numbers\n";
    for (auto n = -max_number; n <= max_number; ++n)
      {
	auto _F_n = __fibonacci(n);
	std::cout << ' ' << std::setw(4) << n
		  << ' ' << std::setw(width) << _F_n
		  << '\n';
      }

    std::cout << "\n\n Fibonacci function\n";
    for (int n = -500; n <= 500; ++n)
      {
	auto nu = n * delnu;
	auto F_nu = __fibonacci(nu);
	std::cout << ' ' << std::setw(width) << nu
		  << ' ' << std::setw(width) << F_nu
		  << '\n';
      }

    //std::cout << "\n\n Fibonacci polynomials\n";
    for (auto n = 0; n <= max_order; ++n)
      {
	std::cout << '\n' << '\n' << ' ' << std::setw(4) << n << '\n';
	const auto del = _Tp{1} / _Tp{10};
	for (int i = -50; i <= 50; ++i)
	  {
	    auto x = del * i;
	    auto _F_n = __fibonacci(n, x);
	    std::cout << ' ' << std::setw(width) << x
		      << ' ' << std::setw(width) << _F_n
		      << '\n';
	  }
      }

    std::cout << "\n\n Lucas numbers\n";
    for (auto n = -max_number; n <= max_number; ++n)
      {
	auto _L_n = __lucas(n);
	std::cout << ' ' << std::setw(4) << n
		  << ' ' << std::setw(width) << _L_n
		  << '\n';
      }

    std::cout << "\n\n Lucas function\n";
    for (int n = -500; n <= 500; ++n)
      {
	auto nu = n * delnu;
	auto L_nu = __lucas(nu);
	std::cout << ' ' << std::setw(width) << nu
		  << ' ' << std::setw(width) << L_nu
		  << '\n';
      }

    //std::cout << "\n\n Lucas polynomials\n";
    for (auto n = 0; n <= max_order; ++n)
      {
	std::cout << '\n' << '\n' << ' ' << std::setw(4) << n << '\n';
	const auto del = _Tp{1} / _Tp{10};
	for (int i = -50; i <= 50; ++i)
	  {
	    auto x = del * i;
	    auto _L_n = __lucas(n, x);
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
