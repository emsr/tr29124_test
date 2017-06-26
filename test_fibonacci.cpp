/*
$HOME/bin_tr29124/bin/g++ -std=gnu++17 -g -Wall -Wextra -Wno-psabi -I. -o test_fibonacci test_fibonacci.cpp -Lwrappers/debug -lquadmath -Lwrappers/debug -lwrap_gsl -lwrap_burkhardt -lgfortran
LD_LIBRARY_PATH=wrappers/debug:$HOME/bin_tr29124/lib64:$LD_LIBRARY_PATH ./test_fibonacci > test_fibonacci.txt

$HOME/bin/bin/g++ -std=gnu++17 -g -Wall -Wextra -I. -o test_fibonacci test_fibonacci.cpp -lquadmath -Lwrappers/debug -lwrap_gsl -lwrap_burkhardt
./test_fibonacci > test_fibonacci.txt
*/

#include <cmath>
#include <limits>
#include <iostream>
#include <iomanip>

/**
 * 
 */
template<typename _UIntTp>
  _UIntTp
  __fibonacci(_UIntTp __n)
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
	__Fn = __Fnm1 + __Fnm2;
	if (__builtin_add_overflow(__Fnm1, __Fnm2, &__Fn))
	    std::__throw_runtime_error(__N("__fibonacci: "
					   "integer overflow"));	  
      }
    return __Fn;
  }

/**
 * 
 */
template<typename _UIntTp, typename _RealTp>
  _RealTp
  __fibonacci(_UIntTp __n, _RealTp __x)
  {
    _UIntTp __Fnm2 = 0;
    if (__n == 0)
      return __Fnm2;
    _UIntTp __Fnm1 = 1;
    if (__n == 1)
      return __Fnm1;
    _UIntTp __Fn = __x * __Fnm1 + __Fnm2;
    for (_UIntTp __k = 3; __k <= __n; ++__k)
      {
	__Fnm2 = __Fnm1;
	__Fnm1 = __Fn;
	__Fn = __x * __Fnm1 + __Fnm2;
      }
    return __Fn;
  }

/**
 * 
 */
template<typename _UIntTp>
  _UIntTp
  __lucas(_UIntTp __n)
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
 * 
 */
template<typename _UIntTp, typename _RealTp>
  _RealTp
  __lucas(_UIntTp __n, _RealTp __x)
  {
    _UIntTp __Lnm2 = 2;
    if (__n == 0)
      return __Lnm2;
    _UIntTp __Lnm1 = __x;
    if (__n == 1)
      return __Lnm1;
    _UIntTp __Ln = __x * __Lnm1 + __Lnm2;
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

    std::cout << "\n Fibonacci numbers\n";
    for (auto n = 0ull; n <= 200; ++n)
      {
	auto _F_n = __fibonacci(n);
	std::cout << ' ' << std::setw(4) << n
		  << ' ' << std::setw(width) << _F_n
		  << '\n';
      }

    std::cout << "\n Fibonacci polynomials\n";
    for (auto n = 0ull; n <= 50; ++n)
      {
	std::cout << '\n' << ' ' << std::setw(4) << n << '\n';
	const auto del = _Tp{1} / _Tp{10};
	for (auto i = 0ull; i <= 50; ++i)
	  {
	    auto x = del * i;
	    auto _F_n = __fibonacci(n, x);
	    std::cout << ' ' << std::setw(width) << x
		      << ' ' << std::setw(width) << _F_n
		      << '\n';
	  }
      }

    std::cout << "\n Lucas numbers\n";
    for (auto n = 0ull; n <= 200; ++n)
      {
	auto _L_n = __lucas(n);
	std::cout << ' ' << std::setw(4) << n
		  << ' ' << std::setw(width) << _L_n
		  << '\n';
      }

    std::cout << "\n Lucas polynomials\n";
    for (auto n = 0ull; n <= 50; ++n)
      {
	std::cout << '\n' << ' ' << std::setw(4) << n << '\n';
	const auto del = _Tp{1} / _Tp{10};
	for (auto i = 0ull; i <= 50; ++i)
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
