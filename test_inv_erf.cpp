/*
$HOME/bin_tr29124/bin/g++ -std=gnu++17 -g -Wall -Wextra -Wno-psabi -I. -o test_inv_erf test_inv_erf.cpp -lquadmath
LD_LIBRARY_PATH=$HOME/bin_tr29124/lib64:$LD_LIBRARY_PATH ./test_inv_erf > test_inv_erf.txt

$HOME/bin/bin/g++ -std=gnu++17 -g -Wall -Wextra -Wno-psabi -I. -o test_inv_erf test_inv_erf.cpp -lquadmath
./test_inv_erf > test_inv_erf.txt
*/

// AAOF pp. 408-409.

#include <vector>
#include <iostream>
#include <iomanip>
#include <limits>
#include <cmath>
#include <complex>
#include <bits/float128_io.h>

  template<typename _Tp>
    _Tp
    __experfc(_Tp __x)
    { return std::exp(__x * __x) * std::erfc(__x); }

  /**
   * Return the inverse error function.
   */
  template<typename _Tp>
    _Tp
    __erfi(_Tp __p)
    {
//std::cout.precision(std::numeric_limits<_Tp>::digits10);
//auto w = 8 + std::cout.precision();
      constexpr auto _S_eps = _Tp{10} * std::numeric_limits<_Tp>::epsilon();
      // Iterate experfc(x^2).
      if (__p < _Tp{0})
	return -__erfi(-__p);
      else
	{
	  auto __x2 = _Tp{25};
	  auto __x2prev2 = _Tp{0}, __x2prev = _Tp{0};
	  const auto _S_max_iter = 500;
	  auto __iter = 0;
	  while (++__iter < _S_max_iter)
	    {
		    
	      __x2prev2 = __x2prev;
	      __x2prev = __x2;
	      __x2 = std::log(__experfc(__x2) / (_Tp{1} - __p));
//std::cout
// << ' ' << std::setw(w) << __x2
// << ' ' << std::setw(w) << __x2 - __x2prev
// << ' ' << std::setw(w) << __x2 - __x2prev2
// << '\n';
	      // If the fraction jumps < 0 just bop it back.
	      if (__x2 < _Tp{0})
		__x2 = -__x2;
	      if (_S_eps * std::abs(__x2) > std::abs(__x2 - __x2prev))
		break;
	      if (_S_eps * std::abs(__x2) > std::abs(__x2 - __x2prev2))
		break;
	    }
	  return std::sqrt(__x2);
	}
    }

/**
 * Test the inverse error function.
 */
template<typename _Tp>
  void
  test_inv_erf()
  {
    //  Build the series coefficients.
    const auto _S_eps = std::numeric_limits<_Tp>::epsilon();
    const auto _S_sqrt_pi = _Tp{1.772453850905516027298167483341145182797Q};

    std::cout.precision(std::numeric_limits<_Tp>::digits10);
    decltype(std::cout.precision()) xw = 22;
    auto w = std::max(xw, 8 + std::cout.precision());

    const int n_max = 250;
    std::vector<_Tp> a;
    a.push_back(1);
    for (int __n = 1; __n < n_max; ++__n)
      {
	auto __atemp = _Tp{0};
	for (int __k = 1; __k <= __n; ++__k)
	  __atemp += _Tp(2 * (__k - 1) + 1) * a[__k - 1]
		 * _Tp(2 * (__n - __k) + 1) * a[__n - __k]
		 / _Tp(__k * (2 * __k - 1));
	__atemp /= _Tp(2 * __n + 1);
	a.push_back(__atemp);
      }

    std::cout << "\n\n" << std::setw(w) << " a_k" << '\n';
    for (auto __aa : a)
      std::cout << ' ' << std::setw(w) << __aa << '\n';

    std::vector<_Tp> c;
    c.push_back(1);
    for (int __n = 1; __n < n_max; ++__n)
      {
	auto __ctemp = _Tp{0.0Q};
	for (int __k = 1; __k <= __n; ++__k)
	  __ctemp += c[__k - 1] * c[__n - __k] / _Tp(__k * (2 * __k - 1));
	c.push_back(__ctemp);
      }
    for (int __n = 1; __n < n_max; ++__n)
      c[__n] /= _Tp(2 * __n + 1);

    std::cout << "\n\n " << std::setw(w) << "c_k" << '\n';
    for (auto __cc : c)
      std::cout << ' ' << std::setw(w) << __cc << '\n';

    std::cout << "\n\n"
	      << ' ' << std::setw(w) << "\"p\""
	      << ' ' << std::setw(w) << "\"inv_erf(p)\""
	      << ' ' << std::setw(w) << "\"erf(inv_erf(p))\""
	      << ' ' << std::setw(w) << "\"erf(inv_erf(p)) - p\""
	      << '\n';
    for (int __i = -100; __i <= 100; ++__i)
      {
	auto __x = __i * _Tp{0.01Q};
	auto __chi = _S_sqrt_pi * __x / _Tp{2};
	auto __chi2 = __chi * __chi;
	auto __chip = __chi;
	auto __inverf = _Tp{0.0Q};
	for (int __k = 0; __k < n_max; ++__k)
	  {
	    auto __term = a[__k] * __chip;
	    __inverf += __term;
	    if (std::abs(__term) < _S_eps * std::abs(__inverf))
	      break;
	    __chip *= __chi2;
	  }
	std::cout << ' ' << std::setw(w) << __x
		  << ' ' << std::setw(w) << __inverf
		  << ' ' << std::setw(w) << erf(__inverf)
		  << ' ' << std::setw(w) << erf(__inverf) - __x
		  << '\n';
      }

    std::cout << "\n\n"
	      << ' ' << std::setw(w) << "\"x\""
	      << ' ' << std::setw(w) << "\"erf(x)\""
	      << ' ' << std::setw(w) << "\"inv_erf(erf(x))\""
	      << ' ' << std::setw(w) << "\"inv_erf(erf(x)) - x\""
	      << ' ' << std::setw(w) << "\"inv_erf(erf(x))\""
	      << ' ' << std::setw(w) << "\"inv_erf(erf(x)) - x\""
	      << '\n';
    for (int __i = -200; __i <= 200; ++__i)
      {
	auto __x = __i * _Tp{0.01Q};
	auto __erfx = erf(__x);
	auto __chi = _S_sqrt_pi * __erfx / _Tp{2};
	auto __chi2 = __chi * __chi;
	auto __chip = __chi;
	auto __inverf = _Tp{0.0Q};
	for (int __k = 0; __k < n_max; ++__k)
	  {
	    auto __term = a[__k] * __chip;
	    __inverf += __term;
	    if (std::abs(__term) < _S_eps * std::abs(__inverf))
	      break;
	    __chip *= __chi2;
	  }
	std::cout << ' ' << std::setw(w) << __x
		  << ' ' << std::setw(w) << __erfx
		  << ' ' << std::setw(w) << __inverf
		  << ' ' << std::setw(w) << __inverf - __x
		  << ' ' << std::setw(w) << __erfi(__erfx)
		  << ' ' << std::setw(w) << __erfi(__erfx) - __x
		  << '\n';
      }
  }

int
main()
{
  std::cout << "\n\n  float\n";
  std::cout << "  =====\n";
  test_inv_erf<float>();

  std::cout << "\n\n  double\n";
  std::cout << "  ======\n";
  test_inv_erf<double>();

  std::cout << "\n\n  long double\n";
  std::cout << "  ===========\n";
  test_inv_erf<long double>();

  std::cout << "\n\n  __float128\n";
  std::cout << "  ==========\n";
  test_inv_erf<__float128>();
}
